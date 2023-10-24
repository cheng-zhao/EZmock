/*******************************************************************************
* hypergeom.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "hypergeom.h"
#include <math.h>

/*******************************************************************************
  Implementation of the Gauss hypergeometric function calculator under
  certain conditions.
  Ref: J. Pearson, Computation of Hypergeometric Functions
       (https://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
*******************************************************************************/

/*============================================================================*\
           Functions for evaluating the Gauss hypergeometric function
\*============================================================================*/

/******************************************************************************
Function `hyp2f1_taylor`:
  Evaluate hyp2f1(a,b,c,z) using Taylor series in the following cases:
     -  |a|, |b| < 50
     -  |c| > 1
     -  |z| < 0.9
  (see Section 4.2 and Table 17 of the reference).
Arguments:
  * `a`,`b`,`c`,`z`:    arguments of the Gauss hypergeometric function;
  * `tol`:      the desired tolerance of the numerical evaluation;
  * `maxiter`:  the maximum allowed number of interations;
  * `res`:      the evaluated function value.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int hyp2f1_taylor(const double a, const double b, const double c,
    const double z, const double tol, const int maxiter, double *res) {
  /* Validate arguments. */
  if (fabs(a) >= 50 || fabs(b) >= 50 || fabs(c) <= 1 || fabs(z) >= 0.9) {
    *res = HUGE_VAL; return 1;
  }

  /* Iterations with Taylor series. */
  int stop_cnt = 0;
  double ci, si;
  ci = si = 1;

  for (int i = 0; i < maxiter; i++) {
    const double j = (double) i;
    /* The (i+1)-th term of the Taylor series. */
    ci *= (a + j) * (b + j) * z / ((c + j) * (j + 1));

    /* Check the termination criteria. */
    if (fabs(ci) < fabs(si) * tol) {
      if (++stop_cnt >= 3) {
        *res = si;
        return 0;
      }
    }
    else stop_cnt = 0;

    /* The sum of the first (i+1) terms. */
    si += ci;
  }

  return 1;
}

/******************************************************************************
Function `hyp2f1_frac`:
  Evaluate hyp2f1(a,b,c,z) as a single fraction by using recurrence relations,
  in the following cases:
     -  |a|, |b| < 50
     -  |c| < 20
     -  |z| < 0.9
  (see Section 4.3 and Table 17 of the reference).
Arguments:
  * `a`,`b`,`c`,`z`:    arguments of the Gauss hypergeometric function;
  * `tol`:      the desired tolerance of the numerical evaluation;
  * `maxiter`:  the maximum allowed number of interations;
  * `res`:      the evaluated function value.
Return:
  Status of the evaluation (zero on success; non-zero on error).
******************************************************************************/
static int hyp2f1_frac(const double a, const double b, const double c,
    const double z, const double tol, const int maxiter, double *res) {
  /* Validate arguments. */
  if (fabs(a) >= 50 || fabs(b) >= 50 || fabs(c) >= 20 || fabs(z) >= 0.9) {
    *res = HUGE_VAL; return 1;
  }

  /* Iterations with Taylor series. */
  int stop_cnt = 0;
  double alpha, beta, gamma, zeta;
  alpha = 0;
  beta = gamma = zeta = 1;

  for (int i = 0; i < maxiter; i++) {
    const double j = (double) i;
    /* Record the function value from the previous iteration. */
    const double prev = zeta;

    /* Recurrence relation. */
    const double fac = (j + 1) * (c + j);
    alpha = (alpha + beta) * fac;
    beta *= (a + j) * (b + j) * z;
    gamma *= fac;
    zeta = (alpha + beta) / gamma;

    /* Check the termination criteria. */
    if (fabs(zeta - prev) < fabs(prev) * tol) {
      if (++stop_cnt >= 3) {
        *res = zeta;
        return 0;
      }
    }
    else stop_cnt = 0;
  }

  return 1;
}


/*============================================================================*\
                    Interface of the hypergeometric function
\*============================================================================*/

/******************************************************************************
Function `hyp2f1`:
  Evaluate hyp2f1(a,b,c,z) in the following cases:
     -  |a|, |b| < 50
     -  z <= 0 or (after transformation) 0 < z <= 1/2
     -  (a - b) and (c - a - b) are not integers
     -  a, b, c, (c - a), (c - b), (a - b) are positive
     -  gamma(b - a) is negative as -1 < b - a < -1/2
  (the last two conditions are for the signs of the gamma functions).
Arguments:
  * `a`,`b`,`c`,`z`:    arguments of the Gauss hypergeometric function;
  * `tol`:      tolerance of the numerical evaluation (default: 1e-12);
  * `maxiter`:  maximum number of interations (default: 10000);
  * `res`:      the evaluated function value.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int hyp2f1(const double a, const double b, const double c, const double z,
    const double tol, const int maxiter, double *res) {
  /* Validate arguments. */
  if (fabs(a) >= 50 || fabs(b) >= 50) {
    *res = HUGE_VAL; return 1;
  }

  /* Apply the transformation formulae (see Table 13 of the reference). */
  if (z < -1) {         /* Eq. (4.16) of the reference */
    double w = 1 / (1 - z);
    double fac1, fac2;
    if (hyp2f1(a, c - b, a - b + 1, w, tol, maxiter, &fac1) ||
        hyp2f1(b, c - a, b - a + 1, w, tol, maxiter, &fac2)) {
      *res = HUGE_VAL; return 1;
    }

    const double gc = lgamma(c);
    const double gfac1 = gc + lgamma(b - a) - lgamma(b) - lgamma(c - a);
    const double gfac2 = gc + lgamma(a - b) - lgamma(a) - lgamma(c - b);

    /* Note that tgamma(b - a) is negative. */
    *res = -pow(w, a) * exp(gfac1) * fac1 + pow(w, b) * exp(gfac2) * fac2;
    return 0;
  }
  else if (z < 0) {     /* Eq. (4.17) of the reference */
    double w = z / (z - 1);
    double fac;
    if (hyp2f1(a, c - b, c, w, tol, maxiter, &fac)) {
      *res = HUGE_VAL; return 1;
    }

    *res = pow(1 - z, -a) * fac;
    return 0;
  }

  /* Now the z value shoule be between 0 and 1/2. */
  if (z < 0 || z > 0.5) {
    *res = HUGE_VAL; return 1;
  }

  /* Use either the Taylor series or single fraction method. */
  if (fabs(c) > 1) return hyp2f1_taylor(a, b, c, z, tol, maxiter, res);
  else return hyp2f1_frac(a, b, c, z, tol, maxiter, res);
}
