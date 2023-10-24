/*******************************************************************************
* config.h: this file is part of the EZmock library.

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

#ifndef __CONFIG_H__
#define __CONFIG_H__

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/

#ifndef M_PI
#define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#define DOUBLE_EPSILON  1e-16   /* ~ machine epsilon for double numbers */
#define DOUBLE_TOL      1e-8    /* tolerance for double number comparison */

#ifdef SINGLE_PREC
#define REAL_MIN        FLT_MIN
#define REAL_TOL        1e-5
#else
#define REAL_MIN        DBL_MIN
#define REAL_TOL        1e-10
#endif

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/

#define EZMOCK_MAX_GRID_SIZE            65536
#define EZMOCK_MAX_LINEAR_PK_VALUE      1073741824
#define EZMOCK_MIN_NUM_TRACER           128
#define EZMOCK_MAX_NUM_TRACER_PER_CELL  1024
#define EZMOCK_MAX_K_FOR_BAO_ENHANCE    0.5

/* Number of histogram bins of log2(rho_t) for speeding up PDF mapping. */
#define EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN      64

/*============================================================================*\
                            Default model parameters
\*============================================================================*/

/* density saturation (rho_sat) */
#define EZMOCK_MODEL_RHO_SATURATION     100
/* width of Gaussian random numbers for the stochaostic bias */
#define EZMOCK_MODEL_GAUSS_WIDTH        10

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/

/* k should vary the fastest to reduce cache miss. */
#define IDX(Ng,i,j,k)      (((size_t) (i) * (Ng) + (j)) * (Ng) + (k))

#endif

