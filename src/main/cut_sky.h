/*******************************************************************************
* cut_sky.h: this file is part of the EZmock program.

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

#ifndef __CUT_SKY_H__
#define __CUT_SKY_H__

#include "EZmock.h"
#include "load_conf.h"
#include "convert_z.h"

/*============================================================================*\
                    Data structure for the cut-sky catalogue
\*============================================================================*/
typedef struct {
  size_t n;             /* number of tracers in the cut-sky catalogue */
  size_t cap;           /* capacity of the cut-sky catalogue          */
  real *x[4];           /* RA, Dec, redshift, real-space redshift     */
  uint16_t *status;     /* bit-code indicating the footprint          */
  real *rand;           /* random number for post-processing          */
  double num_dens;      /* number density of the catalog              */
} CDATA;

/*============================================================================*\
                 Interface for the cutsky catalogue generation
\*============================================================================*/

/******************************************************************************
Function `cutsky`:
  Contruct a cut-sky catalogue from a periodic cubic box.
Arguments:
  * `conf`:     structure for storing configurations;
  * `cvt`:      structure for redshift conversion;
  * `x`:        array for the x coordinates of the input cubic mock catalogue;
  * `y`:        array for the y coordinates of the input cubic mock catalogue;
  * `z`:        array for the z coordinates of the input cubic mock catalogue;
  * `vx`:       array for the peculiar velocities along the x direction;
  * `vy`:       array for the peculiar velocities along the y direction;
  * `vz`:       array for the peculiar velocities along the z direction;
  * `ndata`:    number of tracers in the input cubic mock catalogue.
Return:
  The cut-sky catalogue on success; NULL on error.
******************************************************************************/
CDATA *cutsky(const CONF *conf, const ZCVT *cvt, real *x, real *y, real *z,
    real *vx, real *vy, real *vz, const size_t ndata);

/******************************************************************************
Function `cutsky_init`:
  Deconstruct the cut-sky catalogue;
Arguments:
  * `data`:     instance of the cut-sky catalogue.
******************************************************************************/
void cutsky_destroy(CDATA *data);

#endif
