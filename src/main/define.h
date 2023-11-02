/*******************************************************************************
* define.h: this file is part of the EZmock program.

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

#ifndef __DEFINE_H__
#define __DEFINE_H__

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#define SPEED_OF_LIGHT  299792.458
#ifndef M_PI
#define M_PI            0x1.921fb54442d18p+1    /* PI */
#endif
#ifndef M_E
#define M_E             0x1.5bf0a8b145769p+1    /* e */
#endif

#define DEGREE_2_RAD    0x1.1df46a2529d39p-6    /* M_PI / 180 */
#define RAD_2_DEGREE    0x1.ca5dc1a63c1f8p+5    /* 180 / M_PI */

/*============================================================================*\
                           Definitions of data types
\*============================================================================*/
#define DOUBLE_EPSILON  1e-16   /* ~ machine epsilon for double numbers */
#define DOUBLE_TOL      1e-8    /* tolerance for double number comparison */
#ifdef SINGLE_PREC
#define REAL_OFMT       "%.6g"
#else
#define REAL_OFMT       "%.10g"
#endif

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Default value for unset parameters. */
#define DEFAULT_CODE_NAME               "EZmock"
#define DEFAULT_CONF_FILE               "ezmock.conf"
#define DEFAULT_INTERP_LOG              false
#define DEFAULT_RNG                     PRAND_RNG_MT19937
#define DEFAULT_FIX_AMPLITUDE           false
#define DEFAULT_INVERT_PHASE            false
#define DEFAULT_OMEGA_NU                0
#define DEFAULT_EOS_W                   (-1)
#define DEFAULT_BAO_ENHANCE             0
#define DEFAULT_ATTACH_PARTICLE         false
#define DEFAULT_CUTSKY                  false
#define DEFAULT_OVERWRITE               0
#define DEFAULT_OUTPUT_FORMAT           EZMOCK_OFMT_ASCII
#define DEFAULT_HEADER                  true
#define DEFAULT_VERBOSE                 true

#define EZMOCK_MAX_NGRID                65536

/* Priority of parameters from different sources. */
#define EZMOCK_PRIOR_CMD                5
#define EZMOCK_PRIOR_FILE               1

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
/* Number of subvolumes per side for the output chunks. */
#define EZMOCK_FITS_CHUNK_PER_SIDE      4
/* Parameters for redshift conversion with cubic spline interpolation. */
#define EZMOCK_ZCNVT_DZ         1e-3    /* interval of redshift (z) samples */
#define EZMOCK_ZCNVT_MIN_NSP    100     /* minimum number of z samples      */
#define EZMOCK_ZCNVT_MAX_NSP    100000  /* maximum number of z samples      */
#define EZMOCK_ZCNVT_ORDER      10      /* order for Gauss integration of z */
#define EZMOCK_ZCNVT_EXT        10      /* number of bins extended on edges */
#define EZMOCK_ZCNVT_MAX_V      3000    /* maximum peculiar velocity */

/* Right ascension range that distinguishes NGC and SGC. */
#define DESI_NGC_RA_MIN         90
#define DESI_NGC_RA_MAX         300

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define EZMOCK_PATH_SEP        '/'     /* separator for file paths     */
#define EZMOCK_FILE_CHUNK      1048576 /* chunk size for ASCII file IO */
#define EZMOCK_MAX_CHUNK       INT_MAX /* maximum allowed chunk size   */
/* Initial number of objects allocated for the catalogs.        */
#define EZMOCK_DATA_INIT_NUM   128
/* Maximum number of objects stored for each thread.            */
#define EZMOCK_DATA_THREAD_NUM 1024
/* Comment symbol for the coordinate conversion file.           */
#define EZMOCK_READ_COMMENT    '#'
/* Comment symbol for the output files.                         */
#define EZMOCK_SAVE_COMMENT    '#'

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"             /* Output format for double parameters */

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define EZMOCK_ERR_MEMORY       (-1)
#define EZMOCK_ERR_ARG          (-2)
#define EZMOCK_ERR_FILE         (-3)
#define EZMOCK_ERR_CFG          (-4)
#define EZMOCK_ERR_INIT         (-5)
#define EZMOCK_ERR_COSMO        (-6)
#define EZMOCK_ERR_PK           (-7)
#define EZMOCK_ERR_RHO          (-8)
#define EZMOCK_ERR_ASCII        (-9)
#define EZMOCK_ERR_SAVE         (-10)
#define EZMOCK_ERR_ZCVT         (-11)
#define EZMOCK_ERR_CUTSKY       (-12)
#define EZMOCK_ERR_RNG          (-13)
#define EZMOCK_ERR_UNKNOWN      (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

/* k should vary the fastest to reduce cache miss. */
#define IDX(Ng,i,j,k)      (((size_t) (i) * (Ng) + (j)) * (Ng) + (k))

#endif

