/*******************************************************************************
* load_conf.h: this file is part of the EZmock program.

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

#ifndef __LOAD_CONF_H__
#define __LOAD_CONF_H__

#include <stdbool.h>

/*============================================================================*\
                   Data structure for storing configurations
\*============================================================================*/

typedef struct {
  char *fconf;          /* Name of the configuration file. */
  double Lbox;          /* BOX_SIZE        */
  int Ngrid;            /* NUM_GRID        */
  long Ngal;            /* NUM_TRACER      */
  char *plin;           /* LINEAR_PK       */
  char *pnw;            /* LINEAR_PK_NW    */
  double zpk;           /* REDSHIFT_PK     */
  bool logint;          /* PK_INTERP_LOG   */
  int rng;              /* RAND_GENERATOR  */
  long seed;            /* RAND_SEED       */
  bool fixamp;          /* FIX_AMPLITUDE   */
  bool iphase;          /* INVERT_PHASE    */
  double growth2;       /* GROWTH_PK       */
  double vfac;          /* VELOCITY_FAC    */
  bool eval_growth;     /* Indicate if computing structure growth parameters. */
  double omega_m;       /* OMEGA_M         */
  double omega_nu;      /* OMEGA_NU        */
  double eos_w;         /* DE_EOS_W        */
  double redshift;      /* REDSHIFT        */
  double bao_mod;       /* BAO_ENHANCE     */
  double rho_c;         /* RHO_CRITICAL    */
  double rho_exp;       /* RHO_EXP         */
  double pdf_base;      /* PDF_BASE        */
  double sigv;          /* SIGMA_VELOCITY  */
  bool particle;        /* ATTACH_PARTICLE */
  bool cutsky;          /* CUTSKY          */
  char *y5foot;         /* DESI_Y5_FOOT    */
  char *foot;           /* DESI_NOW_FOOT   */
  char gcap;            /* GALACTIC_CAP    */
  double zmin;          /* CUTSKY_ZMIN     */
  double zmax;          /* CUTSKY_ZMAX     */
  char *output;         /* OUTPUT          */
  char *cutout;         /* OUTPUT_CUTSKY   */
  int ofmt;             /* OUTPUT_FORMAT   */
  bool header;          /* OUTPUT_HEADER   */
  int ovwrite;          /* OVERWRITE       */
  bool verbose;         /* VERBOSE         */
  int nthread;          /* OMP_NUM_THREADS */
} CONF;


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv);

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf);

#endif
