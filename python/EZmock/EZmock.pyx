#
# EZmock.pyx: this file is part of the EZmock python package.
#
# EZmock: Effective Zel'dovich approximation mock generator.
#
# Github repository:
#       https://github.com/cheng-zhao/EZmock
#
# Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from .cEZmock cimport *
from libc.stdlib cimport free
from cython.view cimport array as cvarray
import numpy as np


cdef class EZmock:
  cdef EZMOCK *_ez
  cdef int _err
  cdef int _Ng


  def __cinit__(self, Lbox, Ngrid, rng = 1, seed = 1, nthread = 0):
    '''
    Initialize the EZmock instance.

    Parameters
    ----------
    Lbox: float
        Side length of the cubic simulation box.
    Ngrid: int
        Number of grid cells per box side for the density field.
    rng: Literal[0, 1], optional
        Random number generation algorithm (0: MRG32K3A, 1: MT19937).
        Default: 1.
    seed: int, optional
        Random seed.
        Default: 1.
    nthread: int, optional
        Number of OpenMP threads to be used.
        Default: 0 (maximum number of threads available).
    '''

    # Call the C initializer.
    self._err = 0
    self._Ng = Ngrid
    self._ez = EZmock_init(<double> Lbox, <int> Ngrid, <int> rng,    \
        <uint64_t> seed, <int> nthread, &self._err)
    if self._ez is NULL:
      raise RuntimeError(EZmock_errmsg(self._err).decode('utf-8'))


  def __dealloc__(self):
    '''
    Deconstruct the EZmock instance.
    '''
    EZmock_destroy(self._ez)


  cdef _assert_return_value(self):
    '''
    Raise a runtime error with the error message from the EZmock library
    if an exception occurs.
    '''
    if self._err != 0:
      raise RuntimeError(EZmock_errmsg(self._err).decode('utf-8'))


  cpdef set_growth_params(self, const double growth_pk, const double vel_fac):
    '''
    Set the structure growth parameters directly.

    Parameters
    ----------
    growth_pk: float
        Factor that normalizes the input linear power spectrum at the
        desired redshift of the output catalog: ( D(z_out) / D(z_pk) )**2,
        where D indicates the linear growth factor.
    vel_fac: float
        Factor for converting Lagrangian displacements into physical peculiar
        velocities: f(z_out) * H(z_out) * a(z_out) / h, where f is the linear
        growth rate, H the Hubble parameter, a the scale factor, and h the
        dimensionless Hubble parameter.
    '''
    self._err = EZmock_set_cosmology(self._ez, growth_pk, vel_fac,   \
        False, <double> 0, <double> 0, <double> 0, <double> 0, <double> 0,   \
        &self._err)
    self._assert_return_value()


  cpdef eval_growth_params(self, const double z_out, const double Omega_m,   \
      double z_pk = 0, double Omega_nu = 0, double w = -1):
    '''
    Evaluate the structure growth parameters with cosmological parameters.

    Parameters
    ----------
    z_out: float
        Redshift of the output catalog.
    Omega_m: float
        Matter (without neutrino) density parameter at z = 0.
    z_pk: float, optional
        Redshift of the input linear power spectrum.
        Default: 0.
    Omega_nu: float, optional
        Neutrino density parameter at z = 0.
        Default: 0.
    w: float, optional
        Dark energy equation of state.
        Default: -1.
    '''
    self._err = EZmock_set_cosmology(self._ez, <double> 0,           \
        <double> 0, True, z_out, z_pk, Omega_m, Omega_nu, w, &self._err)
    self._assert_return_value()


  cpdef setup_linear_pk(self, k, Plin, Pnw = None, double BAO_enhance = 0,   \
      bint interp_log = False):
    '''
    Setup the linear power spectrum for generating the initial condition.

    Parameters
    ----------
    k: float numpy array
        An ascending array of wavenumbers for the input power spectra.
    Plin: float numpy array
        The linear power spectrum.
    Pnw: float numpy array, optional
        The linear power spectrum without BAO wiggles (non-wiggle).
        It is only necessary if `BAO_enhance` != 0.
        Default: None.
    BAO_enhance: float, optional
        The BAO enhancement parameter (positive: enhance BAO;
        negative: damp BAO; 0: no effect).
        Default: 0.
    interp_log: boolean, optional
        True for interpolating the power spectra in log scale; False for
        interpolating in linear scale.
        Default: False.

    Prerequisite
    ------------
    Function `set_growth_params` or `eval_growth_params`.
    '''
    try:
      if not isinstance(k, np.ndarray):    k = np.array(k)
      if not isinstance(Plin, np.ndarray): Plin = np.array(Plin)
      k = k.astype(np.float64)
      Plin = Plin.astype(np.float64)
    except:
      raise TypeError('`k` and `Plin` are expected to be float numpy arrays')

    assert k.ndim == 1 and Plin.ndim == 1 and k.shape[0] == Plin.shape[0],   \
        '`k` and `Plin` must be 1d arrays with the same length'

    cdef int nbin = k.shape[0]
    # Create cython memoryviews of the arrays
    if not k.flags['C_CONTIGUOUS']:    k = np.ascontiguousarray(k)
    if not Plin.flags['C_CONTIGUOUS']: Plin = np.ascontiguousarray(Plin)

    if BAO_enhance != 0:
      try:
        if not isinstance(Pnw, np.ndarray): Pnw = np.array(Pnw)
        Pnw = Pnw.astype(np.float64)
      except:
        raise TypeError('`Pnw` is expected to be a float numpy array')

      assert Pnw.ndim == 1 and nbin == Pnw.shape[0],                         \
        '`Pnw` must be an 1d array with the same length as k'

      if not Pnw.flags['C_CONTIGUOUS']: Pnw = np.ascontiguousarray(Pnw)
    else:
      Pnw = None

    cdef double[::1] k_mv = k
    cdef double[::1] Plin_mv = Plin
    cdef double[::1] Pnw_mv = Pnw

    if Pnw_mv is None:  # BAO_enhance = 0: no BAO enhancement
      self._err = EZmock_setup_linear_pk(self._ez, &k_mv[0], nbin,   \
          &Plin_mv[0], NULL, BAO_enhance, interp_log, &self._err)
    else:
      self._err = EZmock_setup_linear_pk(self._ez, &k_mv[0], nbin,   \
          &Plin_mv[0], &Pnw_mv[0], BAO_enhance, interp_log, &self._err)
    self._assert_return_value()


  cpdef create_dens_field_from_ic(self, bint fixamp = False,                 \
      bint iphase = False):
    '''
    Generate the EZmock density field from a fresh initial condition based on
    the input linear power spectra.

    Parameters
    ----------
    fixamp: boolean, optional
        True for fixing the amplitude of the initial condition.
        Default: False.
    iphase: boolean, optional
        True for inverting the phases of the initial condition.
        Default: False.

    Prerequisite
    ------------
    Function `setup_linear_pk`.
    '''
    self._err = EZmock_create_dens_field(self._ez, NULL,            \
        <bint> False, NULL, fixamp, iphase, &self._err)
    self._assert_return_value()


  cpdef create_dens_field_from_wn(self, delta):
    '''
    Generate the EZmock density field given a white noise field and the
    input linear power spectra.

    Parameters
    ----------
    delta: float numpy array
        The configuration-space white noise field, with dimension `Ngrid`**3,
        and C-order indices.

    Prerequisite
    ------------
    Function `setup_linear_pk`.
    '''
    try:
      if not isinstance(delta, np.ndarray): delta = np.array(delta)
      if sizeof(real) == 8: rtype = np.float64
      else: rtype = np.float32
      delta = delta.astype(rtype)
    except:
      raise TypeError('`delta` is expected to be a float numpy array')

    assert delta.ndim == 1 and delta.shape[0] == int(self._Ng) * self._Ng *  \
        self._Ng, f'`delta` must be a 1d array with {self._Ng:d}**3 elements'

    # Make sure that the memory is contiguous.
    if not delta.flags['C_CONTIGUOUS']: delta = np.ascontiguousarray(delta)
    cdef real[::1] delta_view = delta

    self._err = EZmock_create_dens_field(self._ez, NULL,             \
        <bint> False, &delta_view[0], <bint> False, <bint> False, &self._err)
    self._assert_return_value()


  cpdef create_dens_field_from_disp(self, dx, dy, dz, bint deepcopy = False):
    '''
    Generate the EZmock density field given the displacement fields on
    different directions.

    Parameters
    ----------
    dx, dy, dz: float numpy array
        The displacement fields on the x, y, and z directions, respectively,
        with dimensions being `Ngrid`**3, and indices in C order.

    deepcopy: boolean, optional
        True for maintaining a copy of the displacements by the EZmock library,
        otherwise only the references are passed.
        In any case it is the user's responsibility to release memory for the
        displacement field arrays.
        Default: False.
    '''
    try:
      if not isinstance(dx, np.ndarray): dx = np.array(dx)
      if not isinstance(dy, np.ndarray): dy = np.array(dy)
      if not isinstance(dz, np.ndarray): dz = np.array(dz)
      if sizeof(real) == 8: rtype = np.float64
      else: rtype = np.float32
      dx = dx.astype(rtype)
      dy = dy.astype(rtype)
      dz = dz.astype(rtype)
    except:
      raise TypeError('`dx`, `dy`, and `dz` must be float numpy arrays')

    Ntot = int(self._Ng) * self._Ng * self._Ng
    assert dx.ndim == 1 and dy.ndim == 1 and dz.ndim == 1 and                \
        dx.shape[0] == Ntot and dy.shape[0] == Ntot and dz.shape[0] == Ntot, \
        f'`dx`, `dy`, `dz` must be 1d arrays with {self._Ng:d}**3 elements'

    # Make sure that the memory is contiguous.
    if not dx.flags['C_CONTIGUOUS']: dx = np.ascontiguousarray(dx)
    if not dy.flags['C_CONTIGUOUS']: dy = np.ascontiguousarray(dy)
    if not dz.flags['C_CONTIGUOUS']: dz = np.ascontiguousarray(dz)
    cdef real[::1] dx_view = dx
    cdef real[::1] dy_view = dy
    cdef real[::1] dz_view = dz

    cdef real *disp[3]
    disp[0] = &dx_view[0]
    disp[1] = &dy_view[0]
    disp[2] = &dz_view[0]
    self._err = EZmock_create_dens_field(self._ez, disp, deepcopy,   \
        NULL, <bint> False, <bint> False, &self._err)
    self._assert_return_value()


  cpdef populate_tracer(self, const double rho_c, const double rho_exp,      \
      const double pdf_base, const double sigma_v, const size_t ntracer,     \
      bint att_part = False):
    '''
    Generate the EZmock tracer catalog based on the density field and
    effective tracer bias model.

    Parameters
    ----------
    rho_c: float
        Critical dark matter density as the threshold of tracer formation.
        See Eq. (16) of the reference. Allowed range: [0, infty).
    rho_exp: float
        Exponential cut-off of the effective bias model.
        See Eq. (16) of the reference. Allowed range: (0, infty).
    pdf_base: float
        Base number of the power-law tracer probability distribution function.
        See Eq. (18) of the reference. Allowed range: (0, 1).
    sigma_v: float
        Standard deviation for random local peculiar motions.
        See Eq. (24) of the reference. Allowed range: [0, infty).
    ntracer: int
        The expected number of tracers to be generated.
    att_part: boolean, optional
        True for attaching tracers to dark matter partciles whenever possible.
        Default: False.

    Return
    ------
    x, y, z, vx, vy, vz: float numpy arrays
        Arrays for the coordinates and velocities of the tracers.

    Prerequisite
    ------------
    Function `create_dens_field_from_ic` or `create_dens_field_from_wn`
    or `create_dens_field_from_disp`.

    Reference
    ---------
    https://arxiv.org/abs/2007.08997
    '''
    cdef real params[4]
    params[0] = rho_c
    params[1] = rho_exp
    params[2] = pdf_base
    params[3] = sigma_v

    cdef real *x
    cdef real *y
    cdef real *z
    cdef real *vx
    cdef real *vy
    cdef real *vz
    cdef size_t num
    self._err = EZmock_populate_tracer(self._ez, params, ntracer,    \
        att_part, &num, &x, &y, &z, &vx, &vy, &vz, &self._err)
    self._assert_return_value()

    # Convert the cython arrays into numpy arrays.
    cdef cvarray ax  = <real[:num]> x
    cdef cvarray ay  = <real[:num]> y
    cdef cvarray az  = <real[:num]> z
    cdef cvarray avx = <real[:num]> vx
    cdef cvarray avy = <real[:num]> vy
    cdef cvarray avz = <real[:num]> vz

    ax.free_data = True
    ay.free_data = True
    az.free_data = True
    avx.free_data = True
    avy.free_data = True
    avz.free_data = True

    return np.asarray(ax), np.asarray(ay), np.asarray(az),                   \
        np.asarray(avx), np.asarray(avy), np.asarray(avz)


  cpdef populate_tracer_to_file(self, const double rho_c,                    \
      const double rho_exp, const double pdf_base, const double sigma_v,     \
      const size_t ntracer, str fname, bint att_part = False,                \
      double rsd_fac = 0, bint header = True):
    '''
    Generate the EZmock tracer catalog based on the density field and
    effective tracer bias model, then save the catalog to an ASCII file.

    Parameters
    ----------
    rho_c: float
        Critical dark matter density as the threshold of tracer formation.
        See Eq. (16) of the reference. Allowed range: [0, infty).
    rho_exp: float
        Exponential cut-off of the effective bias model.
        See Eq. (16) of the reference. Allowed range: (0, infty).
    pdf_base: float
        Base number of the power-law tracer probability distribution function.
        See Eq. (18) of the reference. Allowed range: (0, 1).
    sigma_v: float
        Standard deviation for random local peculiar motions.
        See Eq. (24) of the reference. Allowed range: [0, infty).
    ntracer: int
        The expected number of tracers to be generated.
    fname: str
        Filename of the output tracer catalog.
    att_part: boolean, optional
        True for attaching tracers to dark matter partciles whenever possible.
        Default: False.
    rsd_fac: float, optional
        Positive for multiplying this factor to the z-velocity for
        redshift-space z coordinate, then write the redshift-space coordinates;
        otherwise, write the real-space coordinates and velocities.
        Default: 0.
    header: boolean, optional
        True for writing the header to the output file.
        Default: True.

    Prerequisite
    ------------
    Function `create_dens_field_from_ic` or `create_dens_field_from_wn`
    or `create_dens_field_from_disp`.

    Reference
    ---------
    https://arxiv.org/abs/2007.08997
    '''
    cdef real params[4]
    params[0] = rho_c
    params[1] = rho_exp
    params[2] = pdf_base
    params[3] = sigma_v

    cdef real *x
    cdef real *y
    cdef real *z
    cdef real *vx
    cdef real *vy
    cdef real *vz
    cdef size_t num
    self._err = EZmock_populate_tracer(self._ez, params, ntracer,    \
        att_part, &num, &x, &y, &z, &vx, &vy, &vz, &self._err)
    self._assert_return_value()

    cdef bytes bfname = fname.encode('utf-8')
    cdef const char *cfname = bfname
    self._err = EZmock_write_ascii(self._ez, x, y, z, vx, vy, vz,    \
        num, rsd_fac, header, cfname, &self._err)
    self._assert_return_value()

    free(x)
    free(y)
    free(z)
    free(vx)
    free(vy)
    free(vz)
