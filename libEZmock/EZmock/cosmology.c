/*******************************************************************************
* cosmology.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "EZmock.h"
#include "structs.h"
#include "hypergeom.h"
#include "config.h"
#include "errmsg.h"
#include <math.h>

/*============================================================================*\
             Function for computing cosmological growth parameters
\*============================================================================*/

/*******************************************************************************
  Numerical evaluations of exact solutions of linear perturbations,
  in a flat-wCDM cosmology.
  In particular, the following parameters (solutions) are concerned:
     -  delta_m = a * hyp2f1((w-1)/2w, -1/(3w), 1-5/(6w), 1-1/Om(a))
     -  f = 1 + 3 * (w-1) / (6w-5) * (1-1/Om(a))
              * hyp2f1((3w-1)/2w, (3w-1)/(3w), 2-5/(6w), 1-1/Om(a))
              / hyp2f1((w-1)/2w, -1/(3w), 1-5/(6w), 1-1/Om(a))
  The allowed ranges of the relevant cosmological parameters are:
     -  0 < a <= 1
     -  w < -1/3
     -  0 < Om < 1
  Ref: https://doi.org/10.1088/1475-7516/2011/10/010
*******************************************************************************/

/******************************************************************************
Function `cosmo_growth`:
  Compute the linear growth factor and growth rate in a flat wCDM cosmology.
Arguments:
  * `a`:        the scale factor;
  * `ainit`:    reference scale factor for the normalization of growth factor;
  * `Om0`:      matter density parameter at present (z=0);
  * `w`:        dark energy equation of state;
  * `tol`:      the desired tolerance of the numerical evaluation;
  * `maxiter`:  the maximum allowed number of interations;
  * `D`:        the evaluated linear growth factor;
  * `f`:        the evaluated linear growth rate;
  * `err`:      integer storing the error code.
******************************************************************************/
static void cosmo_growth(const double a, const double ainit,
    const double Om0, const double w, const double tol, const int maxiter,
    double *D, double *f, int *err) {
  /* Validate the ranges of arguments. */
  if (a <= 0 || a > 1 || ainit <= 0 || ainit > 1 || Om0 <= 0 || Om0 > 1 ||
      w >= -1.0 / 3.0) {
    *D = *f = HUGE_VAL;
    *err = EZMOCK_ERR_PAR_COSMO; return;
  }

  double tl = (tol < DOUBLE_EPSILON) ? HYPERGEOM_DEFAULT_TOL : tol;
  int mx = (maxiter <= 0) ? HYPERGEOM_DEFAULT_MAXITER : maxiter;

  /* Compute the growth rate. */
  const double iw = 1 / w;
  const double f1 = 0.5 * (1 - iw);                     /* (w - 1) / (2w) */
  const double f2 = -0x1.5555555555555p-2 * iw;         /* -1 / (3w) */
  const double f3 = 1 - 0x1.aaaaaaaaaaaabp-1 * iw;      /* 1 - 5 / (6w) */
  const double f4 = (1 - 1 / Om0) * pow(a, -3 * w);

  double hg1, hg2;
  if (hyp2f1(f1, f2, f3, f4, tl, mx, &hg1) ||
      hyp2f1(f1 + 1, f2 + 1, f3 + 1, f4, tl, mx, &hg2)) {
    *D = *f = HUGE_VAL;
    *err = EZMOCK_ERR_COS_GROWTH; return;
  }

  *f = 1 + 3 * (w - 1) * f4 / (6 * w - 5) * hg2 / hg1;

  /* Compute the growth factor. */
  if (a == ainit) {
    *D = 1;
    return;
  }

  const double f5 = (1 - 1 / Om0) * pow(ainit, -3 * w);
  if (hyp2f1(f1, f2, f3, f5, tl, mx, &hg2)) {
    *D = HUGE_VAL;
    *err = EZMOCK_ERR_COS_GROWTH; return;
  }

  *D = (a * hg1) / (ainit * hg2);
}


/*============================================================================*\
          Interface for setting cosmological parameters used by EZmock
\*============================================================================*/

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
    const double Omega_nu, const double w, int *err) {
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return EZMOCK_ERR_ARG_EZ;
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;

  if (!eval) {  /* set the two parameters for EZmock generation directly */
    if (pk_norm <= 0 || fHa <= 0) return (*err = EZMOCK_ERR_PAR_COSMO);
    cosmo->growth2 = pk_norm;
    cosmo->vfac = fHa;
  }
  else {        /* evaluate the parameters for EZmock generation */
    if (z < 0 || z_pk < 0 || Omega_m <= 0 || Omega_m > 1 ||
        Omega_nu < 0 || Omega_nu >= 1 || w >= -1.0 / 3.0)
      return (*err = EZMOCK_ERR_PAR_COSMO);

    double a = 1 / (1 + z);
    double a_pk = 1 / (1 + z_pk);
    double D, f;
    cosmo_growth(a, a_pk, Omega_m, w, HYPERGEOM_DEFAULT_TOL,
        HYPERGEOM_DEFAULT_MAXITER, &D, &f, err);
    if (*err != EZMOCK_SUCCESS) return *err;

    const double Omn = Omega_m + Omega_nu;
    const double Hubble = 100 * sqrt((1 - Omn) * pow(a, -3 * (1 + w))
        + Omn * pow(a, -3));

    cosmo->growth2 = D * D;
    cosmo->vfac = f * Hubble * a;
  }

#ifdef EZMOCK_DEBUG
  printf("\ngrowth2 = %lf    vfac = %lf\n", cosmo->growth2, cosmo->vfac);
#endif
  return EZMOCK_SUCCESS;
}
