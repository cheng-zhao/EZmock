/*******************************************************************************
* EZmock.h: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __EZMOCK_H__
#define __EZMOCK_H__

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef SINGLE_PREC
typedef float   real;
#else
typedef double  real;
#endif

/*============================================================================*\
                         Definition of data structures
\*============================================================================*/

/* Structure for the EZmock instance. */
typedef struct {
  void *conf;           /* simulation settings (box size, grid size)    */
  void *rng;            /* random number generator                      */
  void *cosmo;          /* structure growth parameters                  */
  void *pk;             /* linear power spectrum                        */
  void *mesh;           /* density and displacement fields              */
} EZMOCK;


/*============================================================================*\
                       Interfaces of the EZmock generator
\*============================================================================*/

/******************************************************************************
Function `EZmock_init`:
  Initialise the EZmock instance.
Arguments:
  * `Lbox`:     side length of the simulation box;
  * `Ngrid`:    number of grids per side for the density field;
  * `randgen`:  random number generation algorithm;
  * `seed`:     random seed;
  * `nthread`:  number of OpenMP threads, 0 for omp_get_max_threads();
  * `err`:      integer storing the error code.
Return:
  Interface of EZmock generator.
******************************************************************************/
EZMOCK *EZmock_init(const double Lbox, const int Ngrid, const int randgen,
    const uint64_t seed, const int nthread, int *err);

/******************************************************************************
Function `EZmock_set_cosmology`:
  Set or compute structure growth parameters in a flat wCDM cosmology.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `pk_norm`:  linear P(k) renormalization parameter, i.e., (D(z)/D(z_pk))^2;
  * `fHa`:      the factor for computing peculiar velocity, i.e., f*H(a)*a/h;
  * `eval`:     if true, set the above 2 parameters directly, otherwise
                evaluating them in a flat-wCDM cosmology with the following
                parameters;
  * `z`:        redshift of the EZmock snapshot to be produced;
  * `z_pk`:     redshift at which the input linear power spectrum is normalized;
  * `Omega_m`:  matter (without neutrino) density parameter at present (z=0);
  * `Omega_nu`: neutrino density parameter at present (z=0);
  * `w`:        dark energy equation of state;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on err.
******************************************************************************/
int EZmock_set_cosmology(EZMOCK *ez, const double pk_norm, const double fHa,
    const bool eval, const double z, const double z_pk, const double Omega_m,
    const double Omega_nu, const double w, int *err);

/******************************************************************************
Function `EZmock_setup_linear_pk`:
  Setup the EZmock linear power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `k`:        ascending array for wavenumbers;
  * `n`:        number of `k` bins;
  * `Pk`:       array for the power spectrum at `k`;
  * `Pnw`:      array for the non-wiggle power spectrum at `k`,
                not used if `mBAO` = 0;
  * `mBAO`:     positive: enhance BAO, negative: damp BAO, zero: no effect;
  * `logint`:   indicate if the power spectrum interpolation is going to be
                performed in log scale (log(k) vs. log(pk));
  * `err`:      integer storing the error message.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int EZmock_setup_linear_pk(EZMOCK *ez, const double *k, const int n,
    const double *Pk, const double *Pnw, const double mBAO, const bool logint,
    int *err);

/******************************************************************************
Function `EZmock_create_dens_field`:
  Generate the EZmock density field in three possible ways:
  * using input displacement fields;
  * using an input white noise field (in configuration space);
  * from the input linear power spectra.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `psi`:      if not NULL, set the displacements directly;
  * `deepcopy`: indicate if saving a copy of `psi` or using only references,
                it is the user's responsibility to free `psi`.
  * `delta`:    if not NULL and `psi` is NULL, set a white noise field for
                generating the density field;
  * `fixamp`:   if `psi` and `delta` are NULL, indicate whether the
                initial amplitudes are fixed;
  * `iphase`:   if `psi` and `delta` are NULL, indicate whether the
                initial phases are inverted;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int EZmock_create_dens_field(EZMOCK *ez, real *psi[3], bool deepcopy,
    real *delta, const bool fixamp, const bool iphase, int *err);

/******************************************************************************
Function `EZmock_populate_tracer`:
  Create the tracer catalogue based on the effective bias model.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `params`:   free parameters of the effective bias model:
                * rho_cut: critical dark matter density for structure formation;
                * rho_exp: exponential cut-off of the effective bias model;
                * pdf_base: base number of the power-law tracer PDF;
                * sigma_v: standard deviation of random local peculiar motions;
  * `nexp`:     expected number of tracers to be generated;
  * `att_part`: indicate whether to attach tracers to particles;
  * `ntracer`:  actual number of tracers generated;
  * `x`:        array to be filled x coordinates of the tracers;
  * `y`:        array to be filled y coordinates of the tracers;
  * `z`:        array to be filled z coordinates of the tracers;
  * `vx`:       array to be filled velocities of the tracers along x;
  * `vy`:       array to be filled velocities of the tracers along y;
  * `vz`:       array to be filled velocities of the tracers along z;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int EZmock_populate_tracer(EZMOCK *ez, const real *params, const size_t nexp,
    const bool att_part, size_t *ntracer, real **x, real **y, real **z,
    real **vx, real **vy, real **vz, int *err);

/******************************************************************************
Function `EZmock_destroy`:
  Release memory allocated for the EZmock generator interface.
Arguments:
  * `ez`:       instance of the EZmock generator.
******************************************************************************/
void EZmock_destroy(EZMOCK *ez);

/******************************************************************************
Function `EZmock_errmsg`:
  Return the error message given the error code.
Arguments:
  * `err`:      the error code.
Return:
  Error message corresponding to the input error code.
******************************************************************************/
const char *EZmock_errmsg(const int err);

#endif

