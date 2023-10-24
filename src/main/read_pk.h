/*******************************************************************************
* read_pk.h: this file is part of the EZmock program.

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

#ifndef __READ_PK_H__
#define __READ_PK_H__

#include "load_conf.h"

/*============================================================================*\
              Data structure for the linear matter power spectrum
\*============================================================================*/

typedef struct {
  int n;                /* number of wavenumber bins                  */
  double *k;            /* array for wavenumbers                      */
  double *Plin;         /* array for the linear matter power spectrum */
  double *Pnw;          /* array for the non-wiggle power spectrum    */
} PK;


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
PK *read_pk(const CONF *conf);

/******************************************************************************
Function `pk_destroy`:
  Deconstruct the structure for the linear power spectrum.
Arguments:
  * `pk`:       structure for the linear power spectrum.
******************************************************************************/
void pk_destroy(PK *pk);

#endif
