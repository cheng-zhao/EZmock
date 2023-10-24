/*******************************************************************************
* hypergeom.h: this file is part of the EZmock library.

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

#ifndef __HYPERGEOM_H__
#define __HYPERGEOM_H__

/*******************************************************************************
  Compute the Gauss hypergeometric function hyp2f1(a,b,c,z) with conditions:
    -  |a|, |b| < 5
    -  1 < c < 5
    -  z < 0
    -  1/2 < a - b < 1
    -  |c - a - b| = 1/2
  For simplicity, possible huge numerical errors when (a - b) approaches 1 is
  not dealt with, which rarely occurs in the cosmic structure growth context.
  Ref: J. Pearson, Computation of Hypergeometric Functions
       (https://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
*******************************************************************************/

#define HYPERGEOM_DEFAULT_TOL           1e-12
#define HYPERGEOM_DEFAULT_MAXITER       10000

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
    const double tol, const int maxiter, double *res);

#endif
