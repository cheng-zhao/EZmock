/*******************************************************************************
* dens_field.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "EZmock.h"
#include "linear_pk.h"
#include "config.h"
#include "errmsg.h"

#ifdef OMP
#include <omp.h>
#define OMP_ATOMIC      _Pragma("omp atomic")
#else
#define OMP_ATOMIC
#endif

/*============================================================================*\
                Functions for computing Zel'dovich displacement
\*============================================================================*/

/* Clean all the relevant macros first */
#ifdef EZMOCK_PT_WHITENOISE
  #undef EZMOCK_PT_WHITENOISE
#endif
#ifdef EZMOCK_PT_LOGPK
  #undef EZMOCK_PT_LOGPK
#endif
#ifdef EZMOCK_PT_FIXAMP
  #undef EZMOCK_PT_FIXAMP
#endif
#ifdef EZMOCK_PT_IPHASE
  #undef EZMOCK_PT_IPHASE
#endif

/******************************************************************************
Function `EZmock_ZA_disp<EZMOCK_LOGPK_NAME><EZMOCK_FIXAMP_NAME>
    <EZMOCK_PT_IPHASE>`:
  Generate the Zel'dovich displacement field given the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/

#define EZMOCK_PT_LOGPK         0
#define EZMOCK_PT_FIXAMP        0
#define EZMOCK_PT_IPHASE        0
#include "perturb.c"

#define EZMOCK_PT_LOGPK         0
#define EZMOCK_PT_FIXAMP        0
#define EZMOCK_PT_IPHASE        1
#include "perturb.c"

#define EZMOCK_PT_LOGPK         0
#define EZMOCK_PT_FIXAMP        1
#define EZMOCK_PT_IPHASE        0
#include "perturb.c"

#define EZMOCK_PT_LOGPK         0
#define EZMOCK_PT_FIXAMP        1
#define EZMOCK_PT_IPHASE        1
#include "perturb.c"

#define EZMOCK_PT_LOGPK         1
#define EZMOCK_PT_FIXAMP        0
#define EZMOCK_PT_IPHASE        0
#include "perturb.c"

#define EZMOCK_PT_LOGPK         1
#define EZMOCK_PT_FIXAMP        0
#define EZMOCK_PT_IPHASE        1
#include "perturb.c"

#define EZMOCK_PT_LOGPK         1
#define EZMOCK_PT_FIXAMP        1
#define EZMOCK_PT_IPHASE        0
#include "perturb.c"

#define EZMOCK_PT_LOGPK         1
#define EZMOCK_PT_FIXAMP        1
#define EZMOCK_PT_IPHASE        1
#include "perturb.c"


/******************************************************************************
Function `EZmock_ZA_disp_wn<EZMOCK_LOGPK_NAME>`:
  Generate the Zel'dovich displacement field given the white noise field and
  the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/

#define EZMOCK_PT_WHITENOISE    1
#define EZMOCK_PT_LOGPK         0
#include "perturb.c"

#define EZMOCK_PT_WHITENOISE    1
#define EZMOCK_PT_LOGPK         1
#include "perturb.c"


/*============================================================================*\
                Functions for generating the displacement field
\*============================================================================*/

/******************************************************************************
Function `EZmock_disp_with_rng`:
  Generate the EZmock displacement field with a random number generator.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `fixamp`:   indicate whether the initial amplitudes are fixed;
  * `iphase`:   indicate whether the initial phases are inverted;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int EZmock_disp_with_rng(EZMOCK *ez, const bool fixamp,
    const bool iphase, int *err) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  /* Create FFTW plan. */
#ifdef OMP
  if (FFT_INIT_OMP() == 0) return (*err = EZMOCK_ERR_FFT_INIT);
  FFT_PLAN_OMP(conf->nthread);
#endif
  FFT_PLAN plan = FFT_PLAN_C2R(conf->Ng, conf->Ng, conf->Ng,
     mesh->rhok, mesh->psi[0], FFT_FLAG);

  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  if (pk->log) {
    if (fixamp) {
      if (iphase) EZmock_ZA_disp_logpk_fixamp_iphase(ez, plan, err);
      else EZmock_ZA_disp_logpk_fixamp(ez, plan, err);
    }
    else {
      if (iphase) EZmock_ZA_disp_logpk_iphase(ez, plan, err);
      else EZmock_ZA_disp_logpk(ez, plan, err);
    }
  }
  else {
    if (fixamp) {
      if (iphase) EZmock_ZA_disp_fixamp_iphase(ez, plan, err);
      else EZmock_ZA_disp_fixamp(ez, plan, err);
    }
    else {
      if (iphase) EZmock_ZA_disp_iphase(ez, plan, err);
      else EZmock_ZA_disp(ez, plan, err);
    }
  }

  /* Cleanup FFTW. */
  FFT_DESTROY(plan);
#ifdef OMP
  FFT_CLEAN_OMP();
#else
  FFT_CLEAN();
#endif

  return *err;
}


/******************************************************************************
Function `EZmock_disp_with_whitenoise`:
  Generate the EZmock displacement field with the given white noise field.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `delta`:    the white noise field;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int EZmock_disp_with_whitenoise(EZMOCK *ez, real *delta, int *err) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  /* Create FFTW plans. */
#ifdef OMP
  if (FFT_INIT_OMP() == 0) return (*err = EZMOCK_ERR_FFT_INIT);
  FFT_PLAN_OMP(conf->nthread);
#endif
  FFT_PLAN plan = FFT_PLAN_R2C(conf->Ng, conf->Ng, conf->Ng,
      mesh->psi[0], mesh->rhok2, FFT_FLAG);

  /* Generate the Fourier-space white noise field via inverse FFT. */
  FFT_EXEC_R2C(plan, delta, mesh->rhok2);
  FFT_DESTROY(plan);

  plan = FFT_PLAN_C2R(conf->Ng, conf->Ng, conf->Ng,
      mesh->rhok, mesh->psi[0], FFT_FLAG);

  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  if (pk->log) EZmock_ZA_disp_wn_logpk(ez, plan, err);
  else EZmock_ZA_disp_wn(ez, plan, err);

  /* Cleanup FFTW. */
  FFT_DESTROY(plan);
#ifdef OMP
  FFT_CLEAN_OMP();
#else
  FFT_CLEAN();
#endif

  return *err;
}


/*============================================================================*\
                    Function for computing the density field
\*============================================================================*/

/******************************************************************************
Function `density_field_cic`:
  Compute the density field with cloud-in-cell (CIC) from the displacements.
Arguments:
  * `ez`:       instance of the EZmock generator.
******************************************************************************/
static void density_field_cic(EZMOCK *ez) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  const real igs = conf->Ng / conf->Lbox;       /* inverse grid size */

#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int i = idx_ij / conf->Ng;
      int j = idx_ij % conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      size_t idx0 = idx_ij * conf->Ng;
      for (int k = 0; k < conf->Ng; k++) {
        size_t idx = idx0 + k;
        real x = mesh->psi[0][idx] * igs + i;
        real y = mesh->psi[1][idx] * igs + j;
        real z = mesh->psi[2][idx] * igs + k;

        int x0 = (int) floor(x);
        int y0 = (int) floor(y);
        int z0 = (int) floor(z);

        /* Weights for neighbours. */
        real wx1 = x - x0;
        real wx0 = 1 - wx1;
        real wy1 = y - y0;
        real wy0 = 1 - wy1;
        real wz1 = z - z0;
        real wz0 = 1 - wz1;

        /* Periodic boundary condition. */
        if (x0 < 0) x0 += conf->Ng;
        else if (x0 >= conf->Ng) x0 -= conf->Ng;
        if (y0 < 0) y0 += conf->Ng;
        else if (y0 >= conf->Ng) y0 -= conf->Ng;
        if (z0 < 0) z0 += conf->Ng;
        else if (z0 >= conf->Ng) z0 -= conf->Ng;

        int x1 = (x0 == conf->Ng - 1) ? 0 : x0 + 1;
        int y1 = (y0 == conf->Ng - 1) ? 0 : y0 + 1;
        int z1 = (z0 == conf->Ng - 1) ? 0 : z0 + 1;

OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x0,y0,z0)] += wx0 * wy0 * wz0;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x0,y0,z1)] += wx0 * wy0 * wz1;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x0,y1,z0)] += wx0 * wy1 * wz0;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x0,y1,z1)] += wx0 * wy1 * wz1;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x1,y0,z0)] += wx1 * wy0 * wz0;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x1,y0,z1)] += wx1 * wy0 * wz1;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x1,y1,z0)] += wx1 * wy1 * wz0;
OMP_ATOMIC
        mesh->rho[IDX(conf->Ng,x1,y1,z1)] += wx1 * wy1 * wz1;
      }         /* for k */
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */
}


/*============================================================================*\
                   Interface for generating the density field
\*============================================================================*/

/******************************************************************************
Function `EZmock_create_dens_field`:
  Generate the EZmock density field in three possible ways:
  * using input displacement fields;
  * using an input white noise field (in configuration space);
  * from the input linear power spectra.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `psi`:      if not NULL, set the displacements directly;
  * `deepcopy`: indicate if saving a copy of `psi` or using only references;
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
    real *delta, const bool fixamp, const bool iphase, int *err) {
  /* Validate the arguments. */
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return (*err = EZMOCK_ERR_ARG_EZ);

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  size_t size = (size_t) conf->Ng * conf->Ng * conf->Ng * sizeof(FFT_REAL);
  /* Dereference existing displacements if applicable. */
  if (mesh->psi_ref) mesh->psi[0] = mesh->psi[1] = mesh->psi[2] = NULL;

  /* Create the displacement field. */
  if (psi && psi[0] && psi[1] && psi[2]) {      /* set the displacements */
    if (deepcopy) {
      if (!(mesh->psi[0] || (mesh->psi[0] = FFT_MALLOC(size))) ||
          !(mesh->psi[1] || (mesh->psi[1] = FFT_MALLOC(size))) ||
          !(mesh->psi[2] || (mesh->psi[2] = FFT_MALLOC(size))) ||
          !(mesh->rho || (mesh->rho = FFT_MALLOC(size)))) {
        return (*err = EZMOCK_ERR_MEMORY);
      }
      memcpy(mesh->psi[0], psi[0], size);
      memcpy(mesh->psi[1], psi[1], size);
      memcpy(mesh->psi[2], psi[2], size);
      mesh->psi_ref = false;
    }
    else {
      if (!(mesh->rho || (mesh->rho = FFT_MALLOC(size))))
        return (*err = EZMOCK_ERR_MEMORY);
      /* Override existing displacements if applicable. */
      if (mesh->psi[0]) FFT_FREE(mesh->psi[0]);
      if (mesh->psi[1]) FFT_FREE(mesh->psi[1]);
      if (mesh->psi[2]) FFT_FREE(mesh->psi[2]);
      for (int i = 0; i < 3; i++) mesh->psi[i] = psi[i];
      mesh->psi_ref = true;
    }
  }
  else {                                /* compute the displacement field */
    /* Validate the k range of input linear power spectrum. */
    EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
    if (!pk) return (*err = EZMOCK_ERR_ARG_PK);

    const double kmin = 2 * M_PI / conf->Lbox;
    const double kmax = M_PI * conf->Ng / conf->Lbox * sqrt(3);
    if (pk->kmin > kmin - DOUBLE_TOL || pk->kmax < kmax + DOUBLE_TOL)
      return (*err = EZMOCK_ERR_PK_RANGE);

    /* Allocate memory. */
    int Ngk = (conf->Ng >> 1) + 1;
    size_t csize = (size_t) conf->Ng * conf->Ng * Ngk * sizeof(FFT_CMPLX);
    if (!(mesh->psi[0] || (mesh->psi[0] = FFT_MALLOC(size))) ||
        !(mesh->psi[1] || (mesh->psi[1] = FFT_MALLOC(size))) ||
        !(mesh->psi[2] || (mesh->psi[2] = FFT_MALLOC(size))) ||
        !(mesh->rhok || (mesh->rhok = FFT_MALLOC(csize))) ||
        !(mesh->rhok2 || (mesh->rhok2 = FFT_MALLOC(csize)))) {
      return (*err = EZMOCK_ERR_MEMORY);
    }
    mesh->psi_ref = false;

    if (delta) EZmock_disp_with_whitenoise(ez, delta, err);
    else EZmock_disp_with_rng(ez, fixamp, iphase, err);

    if (*err != EZMOCK_SUCCESS) return *err;
    mesh->rho = (FFT_REAL *) mesh->rhok;
  }

  /* Compute the density field with CIC. */
  memset(mesh->rho, 0, size);
  density_field_cic(ez);

  mesh->fixamp = fixamp;
  mesh->iphase = iphase;

#ifdef EZMOCK_DEBUG
  FILE *fp = fopen("densfield.dat", "w");
  if (!fp ||
      fwrite(mesh->rho, size, 1, fp) != 1) {
    fprintf(stderr, "\nError: cannot write the density field\n");
    exit(1);
  }
  fclose(fp);
/*
  if (!(fp = fopen("psix.dat","w")) ||
      fwrite(mesh->psi[0], size, 1, fp) != 1 || fclose(fp) ||
      !(fp = fopen("psiy.dat","w")) ||
      fwrite(mesh->psi[1], size, 1, fp) != 1 || fclose(fp) ||
      !(fp = fopen("psiz.dat","w")) ||
      fwrite(mesh->psi[2], size, 1, fp) != 1 || fclose(fp)) {
    printf("\nerror saving psi\n");
    exit(1);
  }
*/
#endif
  return EZMOCK_SUCCESS;
}
