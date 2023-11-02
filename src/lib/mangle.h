/*******************************************************************************
* mangle.h: this file is part of the EZmock program.

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

#ifndef __MANGLE_H__
#define __MANGLE_H__

/*******************************************************************************
  A simplified implementation of the `polyid` command of mangle:
  https://space.mit.edu/~molly/mangle/
  which reports a polygon that a point lies inside, given a polygon format
  mask file, optionally with the `simple` pixelization scheme.
  Ref: https://doi.org/10.1111/j.1365-2966.2008.13296.x
*******************************************************************************/

/*============================================================================*\
                          Data structure for polygons
\*============================================================================*/

/* Cap of a polygon, see
   https://space.mit.edu/~molly/mangle/manual/polygon.html */
typedef double POLYCAP[4];

/* Data structure for polygons. */
typedef struct {
  int polyid;           /* polygon ID            */
  int pixel;            /* pixel of the polygon  */
  POLYCAP *cap;         /* caps of the polygon   */
  int ncap;             /* number of caps        */
  double weight;        /* weight of the polygon */
  double area;          /* area of the polygon   */
} POLYGON;

/* Data strucutre for a mask (a collection of polygons). */
typedef struct {
  POLYGON *poly;        /* array of all valid polygons                  */
  int npoly;            /* number of valid polygons                     */
  int res;              /* resolution for the pixelization              */
  int *pix;             /* indices of polygons with different pixel IDs */
} MANGLE;


/*============================================================================*\
                     Interfaces for applying polygon masks
\*============================================================================*/

/******************************************************************************
Function `mangle_init`:
  Initialise a mangle polygon mask.
Arguments:
  * `fname`:    name of the mangle polygon-format mask file;
  * `wmin`:     minimum weight of polygons to be kept;
  * `err`:      an integer for indicating error messages.
Return:
  Address of the structure for the mask.
******************************************************************************/
MANGLE *mangle_init(const char *fname, const double wmin, int *err);

/******************************************************************************
Function `mangle_destroy`:
  Deconstruct the structure for mangle polygon mask.
Arguments:
  * `mask`:     address of the structure for the mask.
******************************************************************************/
void mangle_destroy(MANGLE *mask);

/******************************************************************************
Function `mangle_query`:
  Find a polygon that contains a given point.
Arguments:
  * `mask`:     address of the structure for the mask;
  * `ra`:       right ascension (in degrees) of the point;
  * `dec`:      declination (in degrees) of the point.
Return:
  Pointer to the polygon that contains the point; NULL if no polygon is found.
******************************************************************************/
POLYGON *mangle_query(const MANGLE *mask, const double ra, const double dec);

/******************************************************************************
Function `mangle_errmsg`:
  Produce error message for a given error code.
Arguments:
  * `err`:      the error code
Return:
  A string with error message.
******************************************************************************/
char *mangle_errmsg(const int err);

#endif
