/*******************************************************************************
* convert_z.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "define.h"
#include "convert_z.h"
#include "legauss.h"
#include "cspline.h"

/*============================================================================*\
                       Functions for redshift integration
\*============================================================================*/

/******************************************************************************
Function `zcnvt_integrand`:
  Integrand for the redshift to radial comoving distance conversion.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `widx`:     power index with the dark energy equation of state;
  * `z`:        the redshift to be converted to comoving distance;
Return:
  The integrand for comoving distance integration.
******************************************************************************/
static inline double zcnvt_integrand(const double Omega_m, const double Omega_L,
    const double widx, const double z) {
  double z1 = z + 1;
  double d = Omega_m * pow(z1, 3);
  d += (widx) ? Omega_L * pow(z1, widx) : Omega_L;
  d = SPEED_OF_LIGHT * 0.01 / sqrt(d);
  return d;
}

/******************************************************************************
Function `zcnvt_legauss`:
  Convert redshift to comoving distance using the Legendre-Gauss integration.
Arguments:
  * `Omega_m`:  density parameter of matter at z = 0;
  * `Omega_L`:  density parameter of Lambda at z = 0;
  * `widx`:     power index with the dark energy equation of state;
  * `order`:    order of the integration;
  * `z`:        redshift to be converted to comoving distance;
Return:
  The radial comoving distance on success; HUGE_VAL on error.
******************************************************************************/
static inline double zcnvt_legauss(const double Omega_m, const double Omega_L,
    const double widx, const int order, const double z) {
  /* Variable transformation for integration from 0 to z. */
  double zp = z * 0.5;
  double sum = 0;
  int i;
  for (i = LEGAUSS_IDX(order);
      i < LEGAUSS_IDX(order) + LEGAUSS_LEN_NONZERO(order); i++) {
    /* Look up the abscissas and weights. */
    double x = legauss_x[i];
    double w = legauss_w[i];

    /* Integrate for both positive and negative abscissas. */
    double z1 = zp * (1 + x);
    double z2 = zp * (1 - x);
    sum += w * (zcnvt_integrand(Omega_m, Omega_L, widx, z1) +
        zcnvt_integrand(Omega_m, Omega_L, widx, z2));
  }
  /* For odd orders, there is also the abscissas x = 0. */
  if (order & 1)
    sum += legauss_w[i] * zcnvt_integrand(Omega_m, Omega_L, widx, zp);

  return sum * zp;
}


/*============================================================================*\
                          Functions for interpolation
\*============================================================================*/

/******************************************************************************
Function `create_sample`:
  Create redshift and comoving distance samples for interpolation.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cvt`:      structure for redshift conversion;
  * `zmin`:     minimum redshift of the sample;
  * `dz`:       width of redshift bins for the sample.
******************************************************************************/
static void create_zsample(const CONF *conf, ZCVT *cvt, const double zmin,
    const double dz) {
  /* Pre-compute the power index for Lambda. */
  const double widx = (conf->eos_w == -1) ? 0 : 3 * (1 + conf->eos_w);
  const double omega_m = conf->omega_m + conf->omega_nu;
  const double omega_l = 1 - omega_m;
  const int order = EZMOCK_ZCNVT_ORDER;

  /* Create the samples through Gauss-Legendre integration. */
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
  for (int i = 0; i < cvt->n; i++) {
    cvt->z[i] = zmin + dz * i;
    cvt->d2[i] = zcnvt_legauss(omega_m, omega_l, widx, order, cvt->z[i]);
    cvt->d2[i] *= cvt->d2[i];
  }

#ifdef EZMOCK_DEBUG
  /* Sanity check: squared distances should be in ascending order. */
  for (int i = 1; i < cvt->n; i++) {
    if (cvt->d2[i - 1] >= cvt->d2[i]) {
      P_ERR("the comoving distance sample is not strictly increasing\n");
      exit(EZMOCK_ERR_UNKNOWN);
    }
  }
#endif

  /* Compute the second derivative of z. */
  cspline_ypp(cvt->d2, cvt->z, cvt->n, cvt->zpp);

  /* Compute the minimum and maximum squared distances of interest,
     with tolerance for RSD. */
  double min = conf->zmin -
      EZMOCK_ZCNVT_MAX_V * (1 + conf->zmin) / SPEED_OF_LIGHT;
  if (min < 0) min = 0;
  double max = conf->zmax +
      EZMOCK_ZCNVT_MAX_V * (1 + conf->zmax) / SPEED_OF_LIGHT;
  double dmin = zcnvt_legauss(omega_m, omega_l, widx, order, min);
  double dmax = zcnvt_legauss(omega_m, omega_l, widx, order, max);
  cvt->d2min = dmin * dmin - DOUBLE_TOL;
  cvt->d2max = dmax * dmax + DOUBLE_TOL;

  /* Make sure that the samples are within `d2min` and `d2max`. */
  if (cvt->d2min < cvt->d2[0]) cvt->d2min = cvt->d2[0];
  if (cvt->d2max > cvt->d2[cvt->n - 1]) cvt->d2max = cvt->d2[cvt->n - 1];
}

/******************************************************************************
Function `bin_search`:
  Binary search the x coordinate for interpolation.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `n`:        number of the sample points;
  * `xv`:       x coordinate of the value to be evaluated;
Return:
  Index of the value in the sample to be evaluated.
******************************************************************************/
static inline int bin_search(const double *x, const int n, const double xv) {
  int l = 0;
  int u = n - 1;
  while (l <= u) {
    int i = ((unsigned int) l + (unsigned int) u) >> 1;
    if (i >= n - 1) {
      if (x[n - 1] == xv) return n - 1;
      else return INT_MAX;
    }
    if (x[i + 1] <= xv) l = i + 1;
    else if (x[i] > xv) u = i - 1;
    else return i;
  }
  return INT_MAX;
}


/*============================================================================*\
                      Interface for coordinate conversion
\*============================================================================*/

/******************************************************************************
Function `zcnvt_init`:
  Initialise cubic spline interpolation for converting (squared) comoving
  distances to redshifts.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Data structure for redshift conversion on success; NULL on error.
******************************************************************************/
ZCVT *zcnvt_init(const CONF *conf) {
  printf("Setting up interpolation for redshift conversion ...");
  if (!conf) {
    P_ERR("configurations are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Determine the range and interval of the redshift sample.
     Note that the redshift range should be sufficient for RSD. */
  double zmin = conf->zmin -
      EZMOCK_ZCNVT_MAX_V * (1 + conf->zmin) / SPEED_OF_LIGHT;
  double zmax = conf->zmax +
      EZMOCK_ZCNVT_MAX_V * (1 + conf->zmax) / SPEED_OF_LIGHT;

  /* Further extend the range on both sides for safe interpolation. */
  double dz = EZMOCK_ZCNVT_DZ;
  zmin -= EZMOCK_ZCNVT_EXT * dz;
  zmax += EZMOCK_ZCNVT_EXT * dz;
  if (zmin < 0) zmin = 0;

  size_t nsp = ceil((zmax - zmin) / dz);
  if (nsp < EZMOCK_ZCNVT_MIN_NSP || nsp > EZMOCK_ZCNVT_MAX_NSP) {
    nsp = (nsp < EZMOCK_ZCNVT_MIN_NSP) ? EZMOCK_ZCNVT_MIN_NSP :
        EZMOCK_ZCNVT_MAX_NSP;
    dz = (zmax - zmin) / nsp;
  }

  /* Allocate memory. */
  ZCVT *cvt = malloc(sizeof(ZCVT));
  if (!cvt) {
    P_ERR("failed to allocate memory for redshift conversion\n");
    return NULL;
  }
  cvt->z = cvt->d2 = cvt->zpp = NULL;
  if (!(cvt->z = malloc(nsp * sizeof(double))) ||
      !(cvt->d2 = malloc(nsp * sizeof(double))) ||
      !(cvt->zpp = malloc(nsp * 2 * sizeof(double)))) {
    P_ERR("failed to allocate memory for redshift conversion\n");
    zcnvt_destroy(cvt);
    return NULL;
  }
  cvt->n = nsp;

  create_zsample(conf, cvt, zmin, dz);
  if (conf->verbose)
    printf("  Interpolation sample created with %d points\n", cvt->n);

  printf(FMT_DONE);
  return cvt;
}

/******************************************************************************
Function `convert_z`:
  Convert a squared radial comoving distance to redshift.
Arguments:
  * `cvt`:      structure for redshift conversion;
  * `dist2`:    the input squared radial comoving distance.
Return:
  The corresponding redshift on success; HUGE_VAL on error.
******************************************************************************/
double convert_z(const ZCVT *cvt, const double dist2) {
  int idx = bin_search(cvt->d2, cvt->n, dist2);
  if (idx == INT_MAX) return HUGE_VAL;

  double z = (idx == cvt->n - 1) ? cvt->z[idx] :
      cspline_eval(cvt->d2, cvt->z, cvt->zpp, dist2, idx);
  return z;
}

/******************************************************************************
Function `zcnvt_destroy`:
  Initialise cubic spline interpolation for converting (squared) comoving
  distances to redshifts.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Data structure for coordinate conversion on success; NULL on error.
******************************************************************************/
void zcnvt_destroy(ZCVT *cvt) {
  if (!cvt) return;
  if (cvt->z) free(cvt->z);
  if (cvt->d2) free(cvt->d2);
  if (cvt->zpp) free(cvt->zpp);
  free(cvt);
}
