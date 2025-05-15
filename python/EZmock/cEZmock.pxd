#
# cEZmock.pxd: this file is part of the EZmock python package.
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

from libc.stdint cimport uint64_t

cdef extern from "EZmock.h":
  ctypedef float real
  ctypedef struct EZMOCK:
    pass

  EZMOCK *EZmock_init(const double Lbox, const int Ngrid, const int randgen, \
      const uint64_t seed, const int nthread, int *err)

  int EZmock_set_cosmology(EZMOCK *ez, const double pk_norm,                 \
      const double fHa, const bint eval, const double z, const double z_pk,  \
      const double Omega_m, const double Omega_nu, const double w, int *err)

  int EZmock_setup_linear_pk(EZMOCK *ez, const double *k, const int n,       \
      const double *Pk, const double *Pnw, const double mBAO,                \
      const bint logint, int *err)

  int EZmock_create_dens_field(EZMOCK *ez, real *psi[3], bint deepcopy,      \
      real *delta, const bint fixamp, const bint iphase, int *err)

  int EZmock_populate_tracer(EZMOCK *ez, const real *params,                 \
      const size_t nexp, const bint att_part, size_t *ntracer,               \
      real **x, real **y, real **z, real **vx, real **vy, real **vz, int *err)

  int EZmock_write_ascii(EZMOCK *ez, const real *x, const real *y,           \
      const real *z, const real *vx, const real *vy, const real *vz,         \
      const size_t ntracer, const double rsd_fac, const bint header,         \
      const char *fname, int *err)

  void EZmock_destroy(EZMOCK *ez)

  const char *EZmock_errmsg(const int err)

