/*******************************************************************************
* convert_z.h: this file is part of the EZmock program.

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

#ifndef __CONVERT_Z_H__
#define __CONVERT_Z_H__

#include "load_conf.h"

/*============================================================================*\
                     Data structure for redshift conversion
\*============================================================================*/

typedef struct {
  int n;                /* number of distance conversion samples */
  double d2min;         /* minimum squared distance of interest  */
  double d2max;         /* maximum squared distance of interest  */
  double *z;            /* array for redshift values             */
  double *d2;           /* array for squared comoving distances  */
  double *zpp;          /* second derivative of redshift         */
} ZCVT;


/*============================================================================*\
                       Interfaces for redshift conversion
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
ZCVT *zcnvt_init(const CONF *conf);

/******************************************************************************
Function `convert_z`:
  Convert a squared radial comoving distance to redshift.
Arguments:
  * `cvt`:      structure for redshift conversion;
  * `dist2`:    the input squared radial comoving distance.
Return:
  The corresponding redshift on success; HUGE_VAL on error.
******************************************************************************/
double convert_z(const ZCVT *cvt, const double dist2);

/******************************************************************************
Function `zcnvt_destroy`:
  Initialise cubic spline interpolation for converting (squared) comoving
  distances to redshifts.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Data structure for coordinate conversion on success; NULL on error.
******************************************************************************/
void zcnvt_destroy(ZCVT *cvt);

#endif
