/*******************************************************************************
* read_pk.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_pk.h"
#include "read_file.h"
#include "cspline.h"
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/*============================================================================*\
                   Functions for power spectrum interpolation
\*============================================================================*/

/******************************************************************************
Function `bin_search`:
  Binary search the x coordinate for interpolation.
Arguments:
  * `x`:        x coordinates of the sample points;
  * `n`:        number of the sample points;
  * `xv`:       x coordinate of the value to be evaluated;
  * `istart`:   starting index for the search;
  * `iend`:     ending index for the search.
Return:
  Index of the value in the sample to be evaluated.
******************************************************************************/
static inline int bin_search(const double *x, const double xv,
    const int istart, const int iend) {
  int l = istart;
  int u = iend;
  while (l <= u) {
    int i = ((unsigned int) l + (unsigned int) u) >> 1;
    if (i >= iend) {
      if (x[iend] == xv) return iend;
      else return INT_MAX;
    }
    if (x[i + 1] <= xv) l = i + 1;
    else if (x[i] > xv) u = i - 1;
    else return i;
  }
  return INT_MAX;
}

/******************************************************************************
Function `pk_interp`:
  Interpolate the power spectrum in log scale, given reference wavenumbers.
Arguments:
  * `n`:        number of bins for the power spectrum to be interpolated;
  * `k`:        wavenumbers of the power spectrum to be interpolated;
  * `P`:        the power spectrum to be interpolated;
  * `nref`:     number of bins for the reference wavenumbers;
  * `kref`:     reference wavenumbers to be interpolated the power spectrum on.
Return:
  The interpolated power spectrum on success; NULL on error.
******************************************************************************/
static double *pk_interp(const int n, double *k, double *P, const int nref,
    const double *kref) {
  if (nref < 1) return NULL;

  /* Interpolate in log scale. */
  for (int i = 0; i < n; i++) {
    k[i] = log(k[i]);
    P[i] = log(P[i]);
  }

  /* Allocate memory for the interpolated array. */
  double *res = malloc(sizeof(double) * nref);
  if (!res) return NULL;

  /* Compute the second derivative of log(P) for cubic spline interpolation. */
  double *ypp = malloc(sizeof(double) * n * 2);
  if (!ypp) {
    free(res); return NULL;
  }
  cspline_ypp(k, P, n, ypp);

  /* Process the largest k value first, to reduce the searching range. */
  double kv = log(kref[nref - 1]);
  int end = bin_search(k, kv, 0, n - 1);
  if (end >= n - 1) {
    if (fabs(k[n - 1] - kv) < DOUBLE_EPSILON)
      res[nref - 1] = P[n - 1];
    else {
      free(res); free(ypp); return NULL;
    }
  }
  else res[nref - 1] = exp(cspline_eval(k, P, ypp, kv, end));

  if (nref == 1) {
    free(ypp); return res;
  }

  /* Process the rest k values in ascending order. */
  if (end < n - 1) end++;
  int pos = 0;
  for (int i = 0; i < nref - 1; i++) {
    kv = log(kref[i]);
    pos = bin_search(k, kv, pos, end);
    if (pos >= n - 1) {
      if (fabs(k[n - 1] - kv) < DOUBLE_EPSILON) res[i] = P[n - 1];
      else {
        free(res); free(ypp); return NULL;
      }
    }
    else res[i] = exp(cspline_eval(k, P, ypp, kv, pos));
  }

  free(ypp);
  return res;
}

/******************************************************************************
Function `pk_share_k`:
  Given two power spectrum, interpolating the one with fewer bins in the
  common wavenumber range, for them to share the same wavenumbers.
Arguments:
  * `pk`:       structure storing the resulting linear power spectra;
  * `n1`:       number of bins for the linear power spectrum;
  * `k1`:       wavenumbers of the linear power spectrum;
  * `P1`:       the linear power spectrum;
  * `n2`:       number of bins for the non-wiggle power spectrum;
  * `k2`:       wavenumbers of the non-wiggle power spectrum;
  * `P2`:       the non-wiggle power spectrum;
  * `verb`:     indicate whether to print detailed information.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int pk_share_k(PK *pk, const int n1, double *k1, double *P1,
    const int n2, double *k2, double *P2, const int verb) {
  /* Get the common k range assuming the wavenumbers are in ascending order. */
  const double kmin = (k1[0] < k2[0]) ? k2[0] : k1[0];
  const double kmax = (k1[n1 - 1] < k2[n2 - 1]) ? k1[n1 - 1] : k2[n2 - 1];
  int i1 = n1;          /* starting index of the common range in `k1` */
  int i2 = n2;          /* starting index of the common range in `k2` */
  int cnt1, cnt2;
  cnt1 = cnt2 = 0;

  /* Find the common k range in `k1` and `k2`. */
  for (int i = 0; i < n1; i++) {
    if (k1[i] >= kmin) {
      i1 = i;
      break;
    }
  }
  for (int i = i1; i < n1; i++) {
    if (k1[i] <= kmax) cnt1++;
    else break;
  }

  for (int i = 0; i < n2; i++) {
    if (k2[i] >= kmin) {
      i2 = i;
      break;
    }
  }
  for (int i = i2; i < n2; i++) {
    if (k2[i] <= kmax) cnt2++;
    else break;
  }

  if (cnt1 == 0 || cnt2 == 0) {
    P_ERR("no common k range for the linear and non-wiggle power spectra\n");
    free(k1); free(P1); free(k2); free(P2);
    return EZMOCK_ERR_PK;
  }

  /* Check if interpolation is necessary. */
  if (cnt1 == cnt2) {
    bool same = true;
    for (int i = 0; i < cnt1; i++) {
      if (fabs(k1[i + i1] - k2[i + i2]) > DOUBLE_EPSILON) {
        same = false;
        break;
      }
    }
    if (same) {         /* the two k arrays are identical in the common range */
      if (i1 != 0) {
        for (int i = 0; i < cnt1; i++) {
          k1[i] = k1[i + i1];
          P1[i] = P1[i + i1];
        }
      }
      if (cnt1 != n1) {         /* reduce memory cost if applicable */
        double *tmp = realloc(k1, cnt1 * sizeof(double));
        if (tmp) k1 = tmp;
        tmp = realloc(P1, cnt1 * sizeof(double));
        if (tmp) P1 = tmp;
      }
      pk->n = cnt1;
      pk->k = k1;
      pk->Plin = P1;

      if (i2 != 0) {
        for (int i = 0; i < cnt2; i++) P2[i] = P2[i + i2];
      }
      if (cnt2 != n2) {
        double *tmp = realloc(P2, cnt2 * sizeof(double));
        if (tmp) P2 = tmp;
      }
      pk->Pnw = P2;
      free(k2);

      if (verb)
        printf("  The two power spectra already share the same wavenumbers\n");
      return 0;
    }
  }

  /* Set the common wavenumbers by the array with more bins. */
  if (cnt1 >= cnt2) {           /* set common wavenumbers from k1 */
    if (i1 != 0) {
      for (int i = 0; i < cnt1; i++) {
        k1[i] = k1[i + i1];
        P1[i] = P1[i + i1];
      }
    }
    if (cnt1 != n1) {           /* reduce memory cost if applicable */
      double *tmp = realloc(k1, cnt1 * sizeof(double));
      if (tmp) k1 = tmp;
      tmp = realloc(P1, cnt1 * sizeof(double));
      if (tmp) P1 = tmp;
    }

    pk->n = cnt1;
    pk->k = k1;
    pk->Plin = P1;

    /* Interpolate (k2, P2). */
    if (!(pk->Pnw = pk_interp(n2, k2, P2, pk->n, pk->k))) {
      P_ERR("failed to interpolate the non-wiggle power spectrum\n");
      free(k1); free(P1); free(k2); free(P2);
      return EZMOCK_ERR_PK;
    }
    free(k2);
    free(P2);

    if (verb) printf("  The non-wiggle power spectrum is interpolated\n");
  }
  else {                        /* set common wavenumbers from k2 */
    if (i2 != 0) {
      for (int i = 0; i < cnt2; i++) {
        k2[i] = k2[i + i2];
        P2[i] = P2[i + i2];
      }
    }
    if (cnt2 != n2) {           /* reduce memory cost if applicable */
      double *tmp = realloc(k2, cnt2 * sizeof(double));
      if (tmp) k2 = tmp;
      tmp = realloc(P2, cnt2 * sizeof(double));
      if (tmp) P2 = tmp;
    }

    pk->n = cnt2;
    pk->k = k2;
    pk->Pnw = P2;

    /* Interpolate (k1, P1). */
    if (!(pk->Plin = pk_interp(n1, k1, P1, pk->n, pk->k))) {
      P_ERR("failed to interpolate the linear power spectrum\n");
      free(k1); free(P1); free(k2); free(P2);
      return EZMOCK_ERR_PK;
    }
    free(k1);
    free(P1);

    if (verb) printf("  The linear power spectrum is interpolated\n");
  }

  return 0;
}

/*============================================================================*\
             Interface for reading the linear matter power spectrum
\*============================================================================*/

/******************************************************************************
Function `read_pk`:
  Read the linear ( and non-wiggle if needed) matter power spectrum from file.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Structure storing the linear power spectrum on success; NULL on error.
******************************************************************************/
PK *read_pk(const CONF *conf) {
  printf("Reading the input matter power spectrum ...");
  if (!conf) {
    P_ERR("configurations are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Read the linear matter power spectrum. */
  size_t nlin;
  double *klin, *Plin;
  if (read_ascii_table(conf->plin, &klin, &Plin, &nlin, conf->verbose))
    return NULL;
  if (nlin >= INT_MAX) {
    P_ERR("too many records in the file\n");
    free(klin); free(Plin); return NULL;
  }
  /* The wavenumbers should be in ascending order. */
  for (size_t i = 1; i < nlin; i++) {
    if (klin[i - 1] >= klin[i]) {
      P_ERR("k values are not in ascending order: `%s'\n", conf->plin);
      free(klin); free(Plin); return NULL;
    }
  }
  /* All values should be positive. */
  for (size_t i = 0; i < nlin; i++) {
    if (klin[i] <= 0 || Plin[i] <= 0) {
      P_ERR("non-positive value found in file\n");
      free(klin); free(Plin); return NULL;
    }
  }
  if (conf->verbose)
    printf("  Linear matter power spectrum loaded\n");

  /* Allocate memory. */
  PK *pk = malloc(sizeof(PK));
  if (!pk) {
    P_ERR("failed to allocate memory for storing the power spectrum\n");
    free(klin); free(Plin); return NULL;
  }
  pk->k = klin;
  pk->Plin = Plin;
  pk->Pnw = NULL;

  if (conf->bao_mod == 0) {
    pk->n = nlin;
    printf(FMT_DONE);
    return pk;
  }

  /* Read the linear non-wiggle matter power spectrum. */
  size_t nnw;
  double *knw, *Pnw;
  if (read_ascii_table(conf->pnw, &knw, &Pnw, &nnw, conf->verbose)) {
    pk_destroy(pk); return NULL;
  }
  if (nnw >= INT_MAX) {
    P_ERR("too many records in the file\n");
    free(knw); free(Pnw); pk_destroy(pk); return NULL;
  }
  /* The wavenumbers should be in ascending order. */
  for (size_t i = 1; i < nnw; i++) {
    if (knw[i - 1] >= knw[i]) {
      P_ERR("k values are not in ascending order: `%s'\n", conf->pnw);
      free(knw); free(Pnw); pk_destroy(pk); return NULL;
    }
  }
  /* All values should be positive. */
  for (size_t i = 0; i < nnw; i++) {
    if (knw[i] <= 0 || Pnw[i] <= 0) {
      P_ERR("non-positive value found in file\n");
      free(knw); free(Pnw); pk_destroy(pk); return NULL;
    }
  }
  if (conf->verbose)
    printf("  Linear non-wiggle power spectrum loaded\n");

  /* Interpolate the power spectra if necessary, so they share the same k. */
  if (pk_share_k(pk, nlin, klin, Plin, nnw, knw, Pnw, conf->verbose)) {
    free(pk); return NULL;
  }

  printf(FMT_DONE);
  return pk;
}

/******************************************************************************
Function `pk_destroy`:
  Deconstruct the structure for the linear power spectrum.
Arguments:
  * `pk`:       structure for the linear power spectrum.
******************************************************************************/
void pk_destroy(PK *pk) {
  if (!pk) return;
  if (pk->k) free(pk->k);
  if (pk->Plin) free(pk->Plin);
  if (pk->Pnw) free(pk->Pnw);
  free(pk);
}
