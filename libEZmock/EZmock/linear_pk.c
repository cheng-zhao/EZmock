/*******************************************************************************
* linear_pk.c: this file is part of the EZmock library.

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
#include "cspline.h"

/*============================================================================*\
                   Function for power spectrum interpolation
\*============================================================================*/

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
static inline int bin_search(const double *x, const int n,
    const double xv) {
  int l = 0;
  int u = n - 1;
  while (l <= u) {
    int i = (l + u) >> 1;
    if (i >= n - 1) {
      if (x[n - 1] == xv) return n - 1;
      else return -1;
    }
    if (x[i + 1] <= xv) l = i + 1;
    else if (x[i] > xv) u = i - 1;
    else return i;
  }
  return -1;
}


/*============================================================================*\
                    Interfaces for the linear power spectrum
\*============================================================================*/

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
    int *err) {
  /* Validate arguments. */
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return (*err = EZMOCK_ERR_ARG_EZ);
  if (!k || !Pk || (mBAO != 0 && !Pnw)) return (*err = EZMOCK_ERR_ARG_NULL);
  if (n <= 3 || n > EZMOCK_MAX_LINEAR_PK_VALUE)
    return(*err = EZMOCK_ERR_ARG_NBIN);

  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  if (cosmo->growth2 == HUGE_VAL) return (*err = EZMOCK_ERR_ARG_COSMO);

  /* Validate ranges of k and P(k). */
  for (int i = 0; i < n; i++) {
    if (!isfinite(k[i]) || !isfinite(Pk[i])) {
      return (*err = EZMOCK_ERR_PK_FINITE);
    }
    if (k[i] <= 0 || Pk[i] <= 0) {
      return (*err = EZMOCK_ERR_PK_NONPOS);
    }
    if (i && k[i] <= k[i - 1]) {
      return (*err = EZMOCK_ERR_PK_NONASC);
    }
  }
  if (mBAO != 0) {
    for (int i = 0; i < n; i++) {
      if (!isfinite(Pnw[i])) return (*err = EZMOCK_ERR_PK_FINITE);
      if (Pnw[i] <= 0) return (*err = EZMOCK_ERR_PK_NONPOS);
    }
  }

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  const double kmin = 2 * M_PI / conf->Lbox;
  const double kmax = M_PI * conf->Ng / conf->Lbox * sqrt(3);
  if (k[0] > kmin - DOUBLE_TOL || k[n - 1] < kmax + DOUBLE_TOL)
    return (*err = EZMOCK_ERR_PK_RANGE);

  /* Retrieve the effective k range. */
  int imin, imax;
  for (imin = 0; imin < n; imin++) {
    if (k[imin] > kmin) break;
  }
  for (imax = imin; imax < n; imax++) {
    if (k[imax] > kmax) break;
  }
  /* Extend the effective k range by 3 bins on both sides for cubic splines. */
  imin = (imin < 3) ? 0 : imin - 3;
  imax = (imax + 3 >= n) ? n : imax + 3;
  int nbin = imax - imin;
  if (nbin < 3) return (*err = EZMOCK_ERR_PK_NBIN);

  /* Allocate memory. */
  EZMOCK_PK *pk = calloc(1, sizeof(EZMOCK_PK));
  if (!pk) return (*err = EZMOCK_ERR_MEMORY);

  pk->k = pk->P = pk->P2 = NULL;
  if (!(pk->k = malloc(nbin * sizeof(double))) ||
      !(pk->P = malloc(nbin * sizeof(double))) ||
      !(pk->P2 = malloc(nbin * 2 * sizeof(double)))) {
    EZmock_pk_destroy(pk);
    return (*err = EZMOCK_ERR_MEMORY);
  }

  /* Modify the BAO strength. */
  if (mBAO == 0) {
    for (int i = 0; i < nbin; i++) pk->P[i] = Pk[i + imin] * cosmo->growth2;
  }
  else {
    for (int i = 0; i < nbin; i++) {
      int j = i + imin;
      pk->P[i] = (k[j] <= EZMOCK_MAX_K_FOR_BAO_ENHANCE) ?
          (Pk[j] - Pnw[j]) * exp(k[j] * k[j] * mBAO) + Pnw[j] : Pk[j];
      if (!isfinite(pk->P[i]) || pk->P[i] <= 0) {
        EZmock_pk_destroy(pk);
        return (*err = EZMOCK_ERR_PK_MODBAO);
      }
      pk->P[i] *= cosmo->growth2;
    }
  }

  if (logint) {
    for (int i = 0; i < nbin; i++) {
      pk->k[i] = log(k[i + imin]);
      pk->P[i] = log(pk->P[i]);
    }
  }
  else {
    for (int i = 0; i < nbin; i++) pk->k[i] = k[i + imin];
  }

  pk->n = nbin;
  pk->log = logint;
  pk->kmin = k[imin];
  pk->kmax = k[imax - 1];

  /* Setup the cubic spline interpolation. */
  cspline_ypp(pk->k, pk->P, pk->n, pk->P2);

  /* Override any existing pk. */
  if (ez->pk) EZmock_pk_destroy(ez->pk);
  ez->pk = pk;

  return EZMOCK_SUCCESS;
}

/******************************************************************************
Function `pk_interp`:
  Evaluate the power spectrum at given k with cubic spline interpolation.
Arguments:
  * `pk`:       instance of the power spectrum;
  * `k`:        wavenumber of the power spectrum to be evaluated;
Return:
  The interpolated power spectrum value.
******************************************************************************/
double pk_interp(const EZMOCK_PK *pk, const double k) {
  const int idx = bin_search(pk->k, pk->n, k);
  if (idx == pk->n - 1) return pk->P[idx];
  else if (idx < 0) {
    if (k >= pk->k[pk->n - 1]) return pk->P[pk->n - 1];
    if (k <= pk->k[0]) return pk->P[0];
  }
  return cspline_eval(pk->k, pk->P, pk->P2, k, idx);
}
