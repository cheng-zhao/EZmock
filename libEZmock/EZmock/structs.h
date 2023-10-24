/*******************************************************************************
* structs.h: this file is part of the EZmock library.

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

#ifndef __STRUCTS_H__
#define __STRUCTS_H__

#include "prand.h"
#include "fftw_def.h"
#include <stdbool.h>

/*============================================================================*\
                         Definition of data structures
\*============================================================================*/

/* Structure for the simulation settings. */
typedef struct {
  int nthread;          /* number of OpenMP threads                        */
  int Ng;               /* number of grids per side for the density field  */
  double Lbox;          /* box size                                        */
} EZMOCK_CONF;

/* Structure for the random number generator. */
typedef struct {
  prand_t *rng;         /* interface of the random number generator        */
  prand_rng_enum ran;   /* random number generation algorithm              */
  uint64_t seed;        /* random seed                                     */
} EZMOCK_RNG;

/* Structure for cosmological parameters. */
typedef struct {
  double growth2;       /* squared growth factor for P(k) renormalization  */
  double vfac;          /* factor of linear peculiar velocity: f*H(a)*a/h  */
} EZMOCK_COSMO;

/* Structure for the linear power spectrum. */
typedef struct pk_struct {
  int n;                /* number of k bins                                */
  bool log;             /* indicate if doing log-scale interpolation       */
  double kmin;          /* minimum k                                       */
  double kmax;          /* maximum k                                       */
  double *k;            /* ascending array for wavenumbers                 */
  double *P;            /* array for P(k)                                  */
  double *P2;           /* second derivatives of P(k)                      */
} EZMOCK_PK;

/* Structure for EZmock free parameters. */
typedef struct {
  double bao_enhance;   /* BAO enhancement parameter                       */
  double rho_c;         /* critical density                                */
  double rho_exp;       /* exponential density cutoff                      */
  double pdf_base;      /* base number of the power-law tracer PDF         */
  double sigma_v;       /* standard deviation of the random local motion   */
} EZMOCK_PAR;

/* Structure for the density and displacement fields. */
typedef struct {
  bool rho_replaced;    /* indicate if the DM density field is overwritten */
  bool psi_ref;         /* indicate if the displacements are references    */
  FFT_REAL *psi[3];     /* Lagrangian displacement field                   */
  FFT_CMPLX *rhok;      /* Fourier-space dark matter density field         */
  FFT_CMPLX *rhok2;     /* Fourier-space field for computing displacements */
  FFT_REAL *rho;        /* configuration-space dark matter density field   */
  FFT_REAL *rhot;       /* number density field of tracers                 */
  uint16_t *np;         /* number of particles at each cell                */
  int16_t *nt;          /* number of tracers to be assigned to each cell   */
} EZMOCK_MESH;

/* Structure for the probability distribution function (PDF) of tracers. */
typedef struct {
  size_t ntot;          /* total number of tracers to be generated         */
  int nmax;             /* maximum number of tracers in a given grid cell  */
  size_t ncell;         /* total number of cells that host tracers         */
  size_t *pdf;          /* PDF: number of cells hosting (nmax-i) tracers   */
  size_t *hist;         /* histogram of tracer number density in log2 bins */
  int ilogb_max;        /* maximum ilogb of tracer number density bins     */
} EZMOCK_PDF;


/*============================================================================*\
                 Functions for (de)constructing data structures
\*============================================================================*/

/******************************************************************************
Function `EZmock_rng_destroy`:
  Release memory allocated for the EZmock random number generator.
Arguments:
  * `rng`:      instance of the random number generator.
******************************************************************************/
void EZmock_rng_destroy(void *rng);

/******************************************************************************
Function `EZmock_pk_destroy`:
  Release memory allocated for the EZmock linear power spectrum.
Arguments:
  * `pk`:       instance of the linear power spectrum.
******************************************************************************/
void EZmock_pk_destroy(void *pk);

/******************************************************************************
Function `EZmock_mesh_destroy`:
  Release memory allocated for the EZmock density and displacement fields.
Arguments:
  * `mesh`:     instance of the EZmock fields.
******************************************************************************/
void EZmock_mesh_destroy(void *mesh);

/******************************************************************************
Function `EZmock_pdf_init`:
  Initialise the data structure for the EZmock tracer PDF.
return:
  Instance of the EZmock tracer PDF.
******************************************************************************/
EZMOCK_PDF *EZmock_pdf_init(void);

/******************************************************************************
Function `EZmock_pdf_destroy`:
  Release memory allocated for the EZmock tracer PDF.
Arguments:
  * `pdf`:      instance of the EZmock tracer PDF.
******************************************************************************/
void EZmock_pdf_destroy(EZMOCK_PDF *pdf);

#endif
