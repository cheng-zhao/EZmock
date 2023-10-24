/*******************************************************************************
* save_res.h: this file is part of the EZmock program.

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

#ifndef __SAVE_RES_H__
#define __SAVE_RES_H__

#include "load_conf.h"
#include "EZmock.h"
#include <stdio.h>

typedef enum {
  EZMOCK_OFMT_ASCII = 0,
  EZMOCK_OFMT_FITS
} EZMOCK_OFMT;

/*============================================================================*\
                    Function for saving the tracer catalogue
\*============================================================================*/

/******************************************************************************
Function `save_box`:
  Write the periodic tracer catalogue to file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `x`:        array for the x coordinates;
  * `y`:        array for the y coordinates;
  * `z`:        array for the z coordinates;
  * `vx`:       array for the peculiar velocities along the x direction;
  * `vy`:       array for the peculiar velocities along the y direction;
  * `vz`:       array for the peculiar velocities along the z direction;
  * `ndata`:    number of tracers in the tracer catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_box(CONF *conf, real *x, real *y, real *z,
    real *vx, real *vy, real *vz, const size_t ndata);

#endif
