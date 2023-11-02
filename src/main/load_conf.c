/*******************************************************************************
* load_conf.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "libcfg.h"
#include "prand.h"
#include "save_res.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Check existence of configuration parameters. */
#define CHECK_EXIST_PARAM(name, cfg, var)                       \
  if (!cfg_is_set((cfg), (var))) {                              \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return EZMOCK_ERR_CFG;                                      \
  }
#define CHECK_EXIST_ARRAY(name, cfg, var, num)                  \
  if (!(num = cfg_get_size((cfg), (var)))) {                    \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return EZMOCK_ERR_CFG;                                      \
  }

/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  char *pname = (char *) args;
  if (!pname || !(*pname)) pname = DEFAULT_CODE_NAME;
  else {
    /* Get the basename of the executable. */
    char *end = strrchr(pname, EZMOCK_PATH_SEP);
    if (end) pname = end + 1;
  }

  printf("Usage: %s [OPTION]\n\
Generate the Effective Zel'dovich approximation mock (EZmock).\n\
  -h, --help\n\
        Display this message and exit\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -b, --box             " FMT_KEY(BOX_SIZE) "        Double\n\
        Set the side length of the cubic simulation box\n\
  -g, --grid            " FMT_KEY(NUM_GRID) "        Integer\n\
        Set the number of grid cells per box side for the density field\n\
  -n, --num             " FMT_KEY(NUM_TRACER) "      Long integer\n\
        Set the expected number of tracers to be generated\n\
  -p, --pk              " FMT_KEY(LINEAR_PK) "       String\n\
        Specify the filename of the input linear matter power spectrum\n\
  -P, --pk-nw           " FMT_KEY(LINEAR_PK_NW) "    String\n\
        Specify the filename of the input linear non-wiggle power spectrum\n\
  -Z, --redshift-pk     " FMT_KEY(REDSHIFT_PK) "     Double\n\
        Specify the redshift at which the power spectra are normalized\n\
      --interp-log      " FMT_KEY(PK_INTERP_LOG) "   Boolean\n\
        Indicate whether to interpolate the power spectra in log scale\n\
  -r, --rng             " FMT_KEY(RAND_GENERATOR) "  Integer\n\
        Specify the random number generation algorithm\n\
  -s, --seed            " FMT_KEY(RAND_SEED) "       Long integer\n\
        Set the seed for random number generation\n\
      --fix-amp         " FMT_KEY(FIX_AMPLITUDE) "   Boolean\n\
        Indicate whether to fix the amplitude of the Gaussian random field\n\
      --inv-phase       " FMT_KEY(INVERT_PHASE) "    Boolean\n\
        Indicate whether to invert the phase of the Gaussian random field\n\
      --growth-pk       " FMT_KEY(GROWTH_PK) "       Double\n\
        Set the normalization factor of the input linear power spectra\n\
      --vel-fac         " FMT_KEY(VELOCITY_FAC) "    Double\n\
        Set the factor for computing peculiar velocities from displacements\n\
  -m, --omega-m         " FMT_KEY(OMEGA_M) "         Double\n\
        Set the matter (without neutrino) density parameter at z = 0\n\
  -N, --omega-nu        " FMT_KEY(OMEGA_NU) "        Double\n\
        Set the neutrino density parameter at z = 0\n\
  -W, --de-w            " FMT_KEY(DE_EOS_W) "        Double\n\
        Set the dark energy equation of state\n\
  -z, --redshift        " FMT_KEY(REDSHIFT) "        Double\n\
        Set the redshift of the output simulation catalog\n\
  -E, --bao-enhance     " FMT_KEY(BAO_ENHANCE) "     Double\n\
        Set the BAO enhancement parameter\n\
  -C, --rho-c           " FMT_KEY(RHO_CRITICAL) "    Double\n\
        Set the critical density for structure formation\n\
  -X, --rho-exp         " FMT_KEY(RHO_EXP) "         Double\n\
        Set the exponential cut-off of the effective bias model\n\
  -B, --pdf-base        " FMT_KEY(PDF_BASE) "        Double\n\
        Set the base of the power-law tracer probability distribution function\n\
  -S, --sigma-v         " FMT_KEY(SIGMA_VELOCITY) "  Double\n\
        Set the standard deviation of random local peculiar motions\n\
  -A, --attach-part     " FMT_KEY(ATTACH_PARTICLE) " Boolean\n\
        Indicate whether to attach tracers to particles whenever possible\n\
  -o, --output          " FMT_KEY(OUTPUT) "          String\n\
        Specify the filename of the output tracer catalog\n"
"  -F, --format          " FMT_KEY(OUTPUT_FORMAT) "    Integer\n\
        Specify the format of the output catalog\n"
"      --output-header   " FMT_KEY(OUTPUT_HEADER) "   Boolean\n\
        Indicate whether to save configurations as header of the output file\n\
  -w, --overwrite       " FMT_KEY(OVERWRITE) "       Integer\n\
        Indicate whether to overwrite existing output files\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters.\n\
Github repository: https://github.com/cheng-zhao/EZmock.\n\
Licence: GPLv3.\n",
    pname, DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  (void) args;
  printf("# Configuration file for EZmock (default: `%s').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# For supported random number generation algorithms, see\n\
#         https://github.com/cheng-zhao/prand\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
######################\n\
#  General settings  #\n\
######################\n\n\
BOX_SIZE        = \n\
    # Double-precision number, side length of the periodic box.\n\
NUM_GRID        = \n\
    # Integer, number of grid cells per box side for the density field.\n\
NUM_TRACER      = \n\
    # Long integer, expected number of tracers to be generated.\n\
\n\
#############################################\n\
#  Specifications of the initial condition  #\n\
#############################################\n\
\n\
LINEAR_PK       = \n\
    # String, filename for the input linear matter power spectrum.\n\
    # It must be a text file with the leading two columns being k and P(k).\n\
    # Lines starting with '%c' are omitted.\n\
LINEAR_PK_NW    = \n\
    # String, filename for the input linear non-wiggle matter power spectrum.\n\
    # It is only used of `BAO_ENHANCE` is non-zero.\n\
    # It must be a text file with the leading two columns being k and P_nw(k).\n\
    # Lines starting with '%c' are omitted.\n\
REDSHIFT_PK     = \n\
    # Double-precision number, redshift at which `LINEAR_PK` is normalized.\n\
    # It is only used if `GROWTH_PK` and `VELOCITY_FAC` are not both set.\n\
PK_INTERP_LOG   = \n\
    # Boolean option, indicate whether to interpolate\n\
    # the input power spectrum in log scale (unset: %c).\n\
RAND_GENERATOR  = \n\
    # Integer, specify the random number generator (unset: %d).\n\
    # Allowed values are:\n\
    # * 0: MRG32k3a\n\
    # * 1: Mersenne Twister 19937\n\
    # See https://github.com/cheng-zhao/prand for details.\n\
RAND_SEED       = \n\
    # Integer, specify the seed for the random number generator.\n\
FIX_AMPLITUDE   = \n\
    # Boolean option, true for fixing the amplitude of the initial\n\
    # Gaussian random field (unset: %c).\n\
INVERT_PHASE    = \n\
    # Boolean option, true for inverting the phase of the initial\n\
    # Gaussian random field (unset: %c).\n\
\n\
##################################################\n\
#  Cosmological parameters (assuming flat-wCDM)  #\n\
##################################################\n\
\n\
GROWTH_PK       = \n\
    # Double-precision number, (D(z) / D(z_pk))^2, for the normalization of\n\
    # the input power spectrum.\n\
    # It is only used if `VELOCITY_FAC` is also set.\n\
VELOCITY_FAC    = \n\
    # Factor for computing peculiar velocities from Lagrangian displacements,\n\
    # i.e., f * H(a) * a / h.\n\
    # It is only used if `GROWTH_PK` is also set.\n\
OMEGA_M         = \n\
    # Double-precision number, matter (without neutrino) density parameter\n\
    # at z = 0.\n\
    # It is only used of `GROWTH_PK` and `VELOCITY_FAC` are not both set.\n\
OMEGA_NU        = \n\
    # Double-precision number, neutrino density parameter at z = 0 (unset: "
    OFMT_DBL ").\n\
    # It is only used of `GROWTH_PK` and `VELOCITY_FAC` are not both set.\n\
DE_EOS_W        = \n\
    # Double-precision number, dark energy equation of state: w (unset: "
    OFMT_DBL ").\n\
    # It is only used of `GROWTH_PK` and `VELOCITY_FAC` are not both set.\n\
REDSHIFT        = \n\
    # Double-precision number, redshift of the periodic box.\n\
    # It is only used of `GROWTH_PK` and `VELOCITY_FAC` are not both set.\n\
\n\
##########################################################\n\
#  Parameters for mock generation                        #\n\
#  Ref: section 2.2 of https://arxiv.org/abs/2007.08997  #\n\
##########################################################\n\
\n\
BAO_ENHANCE     = \n\
    # BAO enhancement parameter (eq. 6 of https://arxiv.org/abs/1409.1124).\n\
    # Positive for enhancing BAO; negative for damping BAO (unset: "
    OFMT_DBL ").\n\
RHO_CRITICAL    = \n\
    # Critical density for structure formation (eq. 16).\n\
RHO_EXP         = \n\
    # Expotential cut-off of densities for the bias model (eq. 16).\n\
PDF_BASE        = \n\
    # Base of the power law for PDF mapping (eq. 18).\n\
SIGMA_VELOCITY  = \n\
    # Standard deviation for random local peculiar motions (eq. 24).\n\
ATTACH_PARTICLE = \n\
    # Boolean option, true for attaching tracers to DM particles whenever\n\
    # possible (unset: %c).\n\
\n\
##############################\n\
#  Settings for the outputs  #\n\
##############################\n\
\n\
OUTPUT          = \n\
    # String, name of the output mock catalog.\n"
"OUTPUT_FORMAT   = \n\
    # Integer, format of the output catalog (unset: %d). Allowed values are:\n\
    # %d: ASCII text file"
#ifdef WITH_CFITSIO
";\n    # %d: FITS binary table"
  #if EZMOCK_FITS_CHUNK_PER_SIDE <= 1
    ".\n"
  #else
    ", with %d subboxes.\n"
  #endif
#else
    ".\n"
#endif
"OUTPUT_HEADER   = \n\
    # Boolean option, true for saving configurations in the output (unset: %c).\n\
OVERWRITE       = \n\
    # Integer, indicate whether to overwrite existing files (unset: %d).\n\
    # Allowed values are:\n\
    # * 0: quit the program when an output file exist;\n\
    # * positive: force overwriting output files whenever possible;\n\
    # * negative: notify at most this number of times for existing files.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether to show detailed outputs (unset: %c).\n",
  DEFAULT_CONF_FILE, EZMOCK_READ_COMMENT, EZMOCK_READ_COMMENT,
  DEFAULT_INTERP_LOG ? 'T' : 'F', DEFAULT_RNG,
  DEFAULT_FIX_AMPLITUDE ? 'T' : 'F', DEFAULT_INVERT_PHASE ? 'T' : 'F',
  (double) DEFAULT_OMEGA_NU, (double) DEFAULT_EOS_W,
  (double) DEFAULT_BAO_ENHANCE, DEFAULT_ATTACH_PARTICLE ? 'T' : 'F',
  DEFAULT_OUTPUT_FORMAT, EZMOCK_OFMT_ASCII,
#ifdef WITH_CFITSIO
  EZMOCK_OFMT_FITS,
#if EZMOCK_FITS_CHUNK_PER_SIDE > 1
  EZMOCK_FITS_CHUNK_PER_SIDE * EZMOCK_FITS_CHUNK_PER_SIDE *
  EZMOCK_FITS_CHUNK_PER_SIDE,
#endif
#endif
  DEFAULT_HEADER ? 'T' : 'F',
  DEFAULT_OVERWRITE, DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fconf = conf->plin = conf->pnw = conf->output = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const cfg_func_t funcs[] = {
    {   'h',        "help",             usage,          NULL},
    {   't',    "template",     conf_template,          NULL}
  };

  /* Configuration parameters. */
  const cfg_param_t params[] = {
    {'c', "conf"         , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf   },
    {'b', "box"          , "BOX_SIZE"       , CFG_DTYPE_DBL , &conf->Lbox    },
    {'g', "grid"         , "NUM_GRID"       , CFG_DTYPE_INT , &conf->Ngrid   },
    {'n', "num"          , "NUM_TRACER"     , CFG_DTYPE_LONG, &conf->Ngal    },
    {'p', "pk"           , "LINEAR_PK"      , CFG_DTYPE_STR , &conf->plin    },
    {'P', "pk-nw"        , "LINEAR_PK_NW"   , CFG_DTYPE_STR , &conf->pnw     },
    {'Z', "redshift-pk"  , "REDSHIFT_PK"    , CFG_DTYPE_DBL , &conf->zpk     },
    { 0 , "interp-log"   , "PK_INTERP_LOG"  , CFG_DTYPE_BOOL, &conf->logint  },
    {'r', "rng"          , "RAND_GENERATOR" , CFG_DTYPE_INT , &conf->rng     },
    {'s', "seed"         , "RAND_SEED"      , CFG_DTYPE_LONG, &conf->seed    },
    { 0 , "fix-amp"      , "FIX_AMPLITUDE"  , CFG_DTYPE_BOOL, &conf->fixamp  },
    { 0 , "inv-phase"    , "INVERT_PHASE"   , CFG_DTYPE_BOOL, &conf->iphase  },
    { 0 , "growth-pk"    , "GROWTH_PK"      , CFG_DTYPE_DBL , &conf->growth2 },
    { 0 , "vel-fac"      , "VELOCITY_FAC"   , CFG_DTYPE_DBL , &conf->vfac    },
    {'m', "omega-m"      , "OMEGA_M"        , CFG_DTYPE_DBL , &conf->omega_m },
    {'N', "omega-nu"     , "OMEGA_NU"       , CFG_DTYPE_DBL , &conf->omega_nu},
    {'W', "de-w"         , "DE_EOS_W"       , CFG_DTYPE_DBL , &conf->eos_w   },
    {'z', "redshift"     , "REDSHIFT"       , CFG_DTYPE_DBL , &conf->redshift},
    {'E', "bao-enhance"  , "BAO_ENHANCE"    , CFG_DTYPE_DBL , &conf->bao_mod },
    {'C', "rho-c"        , "RHO_CRITICAL"   , CFG_DTYPE_DBL , &conf->rho_c   },
    {'X', "rho-exp"      , "RHO_EXP"        , CFG_DTYPE_DBL , &conf->rho_exp },
    {'B', "pdf-base"     , "PDF_BASE"       , CFG_DTYPE_DBL , &conf->pdf_base},
    {'S', "sigma-v"      , "SIGMA_VELOCITY" , CFG_DTYPE_DBL , &conf->sigv    },
    {'A', "attach-part"  , "ATTACH_PARTICLE", CFG_DTYPE_BOOL, &conf->particle},
    {'o', "output"       , "OUTPUT"         , CFG_DTYPE_STR , &conf->output  },
    {'F', "format"       , "OUTPUT_FORMAT"  , CFG_DTYPE_INT , &conf->ofmt    },
    { 0 , "output-header", "OUTPUT_HEADER"  , CFG_DTYPE_BOOL, &conf->header  },
    {'w', "overwrite"    , "OVERWRITE"      , CFG_DTYPE_INT , &conf->ovwrite },
    {'v', "verbose"      , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose }
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, sizeof(funcs) / sizeof(funcs[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, sizeof(params) / sizeof(params[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, EZMOCK_PRIOR_CMD, &optidx)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (!cfg_is_set(cfg, &conf->fconf)) conf->fconf = DEFAULT_CONF_FILE;
  if (cfg_read_file(cfg, conf->fconf, EZMOCK_PRIOR_FILE)) P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR("the input " FMT_KEY(%s) " is not set\n", key);
    return EZMOCK_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
    return EZMOCK_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_output`:
  Check whether an output file can be written.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_output(char *fname, const char *key, const int ovwrite) {
  if (!fname || *fname == '\0') {
    P_ERR("the output " FMT_KEY(%s) " is not set\n", key);
    return EZMOCK_ERR_CFG;
  }

  /* Check if the file exists. */
  if (!access(fname, F_OK)) {
    /* not overwriting */
    if (ovwrite == 0) {
      P_ERR("the output " FMT_KEY(%s) " exists: `%s'\n", key, fname);
      return EZMOCK_ERR_FILE;
    }
    /* force overwriting */
    else if (ovwrite > 0) {
      P_WRN("the output " FMT_KEY(%s) " will be overwritten: `%s'\n",
          key, fname);
    }
    /* ask for decision */
    else {
      P_WRN("the output " FMT_KEY(%s) " exists: `%s'\n", key, fname);
      char confirm = 0;
      for (int i = 0; i != ovwrite; i--) {
        fprintf(stderr, "Are you going to overwrite it? (y/n): ");
        if (scanf("%c", &confirm) != 1) continue;
        int c;
        while((c = getchar()) != '\n' && c != EOF) continue;
        if (confirm == 'n') {
          P_ERR("cannot write to the file\n");
          return EZMOCK_ERR_FILE;
        }
        else if (confirm == 'y') break;
      }
      if (confirm != 'y') {
        P_ERR("too many failed inputs\n");
        return EZMOCK_ERR_FILE;
      }
    }

    /* Check file permission for overwriting. */
    if (access(fname, W_OK)) {
      P_ERR("cannot write to file `%s'\n", fname);
      return EZMOCK_ERR_FILE;
    }
  }
  /* Check the path permission. */
  else {
    char *end;
    if ((end = strrchr(fname, EZMOCK_PATH_SEP)) != NULL) {
      *end = '\0';
      if (access(fname, X_OK)) {
        P_ERR("cannot access the directory `%s'\n", fname);
        return EZMOCK_ERR_FILE;
      }
      *end = EZMOCK_PATH_SEP;
    }
  }
  return 0;
}

#ifdef WITH_CFITSIO
/******************************************************************************
Function `check_output_fits`:
  Check the extension of the output filename, and if the files can be written.
Arguments:
  * `fname`:    the filename to be tested;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files.
Return:
  Zero if the filename ends with the desired extension name.
******************************************************************************/
static int check_output_fits(char *fname, const char *key, const int ovwrite) {
  if (!fname || *fname == '\0') {
    P_ERR("the output " FMT_KEY(%s) " is not set\n", key);
    return EZMOCK_ERR_CFG;
  }

  /* Check if the filename ends with ".fits" or ".fits.gz". */
  size_t len = strlen(fname);
  size_t elen = strlen(".fits");
  if (len < elen || strncmp(fname + len - elen, ".fits", elen)) {
    elen = strlen(".fits.gz");
    if (len < elen || strncmp(fname + len - elen, ".fits.gz", elen)) {
      P_ERR(FMT_KEY(%s) " must end with \".fits\" or \".fits.gz\"\n", key);
      return EZMOCK_ERR_CFG;
    }
  }

  if (EZMOCK_FITS_CHUNK_PER_SIDE <= 1) {
    return check_output(fname, key, ovwrite);
  }
  else if (EZMOCK_FITS_CHUNK_PER_SIDE > 10) {
    P_ERR("too many file chunks.\n"
        "Please reset `EZMOCK_FITS_CHUNK_PER_SIDE` in `src/main/define.h`\n");
    return EZMOCK_ERR_CFG;
  }
  else {        /* save catalogue in chunks */
    int num = EZMOCK_FITS_CHUNK_PER_SIDE * EZMOCK_FITS_CHUNK_PER_SIDE *
        EZMOCK_FITS_CHUNK_PER_SIDE;
    char *fchunk = malloc(len + 5);     /* add digits and null termination */
    if (!fchunk) {
      P_ERR("failed to allocate memory for checking " FMT_KEY(%s) "\n", key);
      return EZMOCK_ERR_MEMORY;
    }
    memcpy(fchunk, fname, len - elen);

    for (int i = 0; i < num; i++) {
      int n = snprintf(fchunk + len - elen, 5, ".%d", i);
      if (n <= 0) {
        P_ERR("failed to generate filename for chunks of the output catalog\n");
        free(fchunk);
        return EZMOCK_ERR_UNKNOWN;
      }
      strncpy(fchunk + len - elen + n, fname + len - elen, elen + 1);

      if ((n = check_output(fchunk, key, ovwrite))) return n;
    }

    free(fchunk);
  }

  return 0;
}
#endif

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int e;

  /* Check BOX_SIZE. */
  CHECK_EXIST_PARAM(BOX_SIZE, cfg, &conf->Lbox);
  if (conf->Lbox <= 0) {
    P_ERR(FMT_KEY(BOX_SIZE) " must be > 0\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check NUM_GRID. */
  CHECK_EXIST_PARAM(NUM_GRID, cfg, &conf->Ngrid);
  if (conf->Ngrid <= 1) {
    P_ERR(FMT_KEY(NUM_GRID) " must be > 1\n");
    return EZMOCK_ERR_CFG;
  }
  else if (conf->Ngrid > EZMOCK_MAX_NGRID) {
    P_ERR(FMT_KEY(NUM_GRID) " cannot exceed the pre-defined limit: %d\n",
        EZMOCK_MAX_NGRID);
    return EZMOCK_ERR_CFG;
  }
  else if (conf->Ngrid & 1) {
    P_ERR(FMT_KEY(NUM_GRID) " must be an even number\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check NUM_TRACER. */
  CHECK_EXIST_PARAM(NUM_TRACER, cfg, &conf->Ngal);
  if (conf->Ngal <= 0) {
    P_ERR(FMT_KEY(NUM_TRACER) " must be > 0\n");
    return EZMOCK_ERR_CFG;
  }
  else if (conf->Ngal / (long) conf->Ngrid > (long) conf->Ngrid * conf->Ngrid) {
    P_ERR(FMT_KEY(NUM_TRACER) " must be <= " FMT_KEY(NUM_GRID) "^3\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check LINEAR_PK. */
  CHECK_EXIST_PARAM(LINEAR_PK, cfg, &conf->plin);
  if ((e = check_input(conf->plin, "LINEAR_PK"))) return e;

  /* Check PK_INTERP_LOG. */
  if (!cfg_is_set(cfg, &conf->logint)) conf->logint = DEFAULT_INTERP_LOG;

  /* Check RAND_GENERATOR. */
  if (!cfg_is_set(cfg, &conf->rng)) conf->rng = DEFAULT_RNG;
  switch (conf->rng) {
    case PRAND_RNG_MT19937:
    case PRAND_RNG_MRG32K3A:
      break;
    default:
      P_ERR("invalid " FMT_KEY(RAND_GENERATOR) ": %d\n", conf->rng);
      return EZMOCK_ERR_CFG;
  }

  /* Check RAND_SEED. */
  CHECK_EXIST_PARAM(RAND_SEED, cfg, &conf->seed);
  if (conf->seed <= 0) {
    P_ERR(FMT_KEY(RAND_SEED) " must be > 0\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check FIX_AMPLITUDE. */
  if (!cfg_is_set(cfg, &conf->fixamp)) conf->fixamp = DEFAULT_FIX_AMPLITUDE;

  /* Check INVERT_PHASE. */
  if (!cfg_is_set(cfg, &conf->iphase)) conf->iphase = DEFAULT_INVERT_PHASE;

  /* Check GROWTH_PK and VELOCITY_FAC. */
  if (cfg_is_set(cfg, &conf->growth2) && cfg_is_set(cfg, &conf->vfac)) {
    conf->eval_growth = false;

    if (conf->growth2 <= 0) {
      P_ERR(FMT_KEY(GROWTH_PK) " must be > 0\n");
      return EZMOCK_ERR_CFG;
    }
    if (conf->vfac <= 0) {
      P_ERR(FMT_KEY(VELOCITY_FAC) " must be > 0\n");
      return EZMOCK_ERR_CFG;
    }
  }
  else {
    conf->eval_growth = true;

    /* Check REDSHIFT_PK. */
    CHECK_EXIST_PARAM(REDSHIFT_PK, cfg, &conf->zpk);
    if (conf->zpk < 0) {
      P_ERR(FMT_KEY(REDSHIFT_PK) " must be >= 0\n");
      return EZMOCK_ERR_CFG;
    }

    /* Check OMEGA_M. */
    CHECK_EXIST_PARAM(OMEGA_M, cfg, &conf->omega_m);
    if (conf->omega_m <= 0 || conf->omega_m > 1) {
      P_ERR(FMT_KEY(OMEGA_M) " must be > 0 and <= 1\n");
      return EZMOCK_ERR_CFG;
    }

    /* Check OMEGA_NU. */
    if (!cfg_is_set(cfg, &conf->omega_nu)) conf->omega_nu = DEFAULT_OMEGA_NU;
    if (conf->omega_nu < 0 || conf->omega_nu >= 1) {
      P_ERR(FMT_KEY(OMEGA_NU) " must be >= 0 and < 1\n");
      return EZMOCK_ERR_CFG;
    }

    /* Check DE_EOS_W. */
    if (!cfg_is_set(cfg, &conf->eos_w)) conf->eos_w = DEFAULT_EOS_W;
    if (conf->eos_w > -0x1.5555555555555p-2) {
      P_ERR(FMT_KEY(DE_EOS_W) " must be <= -1/3\n");
      return EZMOCK_ERR_CFG;
    }

    /* Check REDSHIFT. */
    CHECK_EXIST_PARAM(REDSHIFT, cfg, &conf->redshift);
    if (conf->redshift < 0) {
      P_ERR(FMT_KEY(REDSHIFT) " must be >= 0\n");
      return EZMOCK_ERR_CFG;
    }
  }

  /* Check BAO_ENHANCE. */
  if (!cfg_is_set(cfg, &conf->bao_mod)) conf->bao_mod = DEFAULT_BAO_ENHANCE;
  if (conf->bao_mod != 0) {
    /* Check LINEAR_PK_NW. */
    CHECK_EXIST_PARAM(LINEAR_PK_NW, cfg, &conf->pnw);
    if ((e = check_input(conf->pnw, "LINEAR_PK_NW"))) return e;
  }

  /* Check RHO_CRITICAL. */
  CHECK_EXIST_PARAM(RHO_CRITICAL, cfg, &conf->rho_c);
  if (conf->rho_c < 0) {
    P_ERR(FMT_KEY(RHO_CRITICAL) " must be >= 0\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check RHO_EXP. */
  CHECK_EXIST_PARAM(RHO_EXP, cfg, &conf->rho_exp);
  if (conf->rho_exp <= 0) {
    P_ERR(FMT_KEY(RHO_EXP) " must be > 0\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check PDF_BASE. */
  CHECK_EXIST_PARAM(PDF_BASE, cfg, &conf->pdf_base);
  if (conf->pdf_base <= 0 || conf->pdf_base >= 1) {
    P_ERR(FMT_KEY(PDF_BASE) " must be > 0 and < 1\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check SIGMA_VEL. */
  CHECK_EXIST_PARAM(SIGMA_VEL, cfg, &conf->sigv);
  if (conf->sigv < 0) {
    P_ERR(FMT_KEY(SIGMA_VEL) " must be >= 0\n");
    return EZMOCK_ERR_CFG;
  }

  /* Check ATTACH_PARTICLE. */
  if (!cfg_is_set(cfg, &conf->particle))
    conf->particle = DEFAULT_ATTACH_PARTICLE;

  /* Check OVERWRITE. */
  if (!cfg_is_set(cfg, &conf->ovwrite)) conf->ovwrite = DEFAULT_OVERWRITE;

  /* Check OUTPUT. */
  CHECK_EXIST_PARAM(OUTPUT, cfg, &conf->output);

  /* Check OUTPUT_FORMAT. */
  if (!cfg_is_set(cfg, &conf->ofmt)) conf->ofmt = DEFAULT_OUTPUT_FORMAT;
  switch (conf->ofmt) {
    case EZMOCK_OFMT_ASCII:
      /* Check OUTPUT. */
      if ((e = check_output(conf->output, "OUTPUT", conf->ovwrite))) return e;
      break;
    case EZMOCK_OFMT_FITS:
#ifdef WITH_CFITSIO
      /* Check OUTPUT. */
      if ((e = check_output_fits(conf->output, "OUTPUT", conf->ovwrite)))
        return e;
      break;
#else
      P_ERR("FITS-format support not enabled.\n"
          "Please recompile the code with `WITH_FITS = T` in `options.mk`\n");
      return EZMOCK_ERR_CFG;
#endif
    default:
      P_ERR("invalid " FMT_KEY(OUTPUT_FORMAT) ": %d\n", conf->ofmt);
      return EZMOCK_ERR_CFG;
  }

  /* Check OUTPUT_HEADER. */
  if (!cfg_is_set(cfg, &conf->header)) conf->header = DEFAULT_HEADER;

  /* Check VERBOSE. */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;

#ifdef OMP
  conf->nthread = omp_get_max_threads();
#else
  conf->nthread = 1;
#endif
  if (conf->nthread <= 0) {
    P_ERR("invalid number of OpenMP threads\n");
    return EZMOCK_ERR_CFG;
  }

  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations.
******************************************************************************/
static void conf_print(const CONF *conf) {
  /* Configuration file */
  printf("\n  CONFIG_FILE     = %s", conf->fconf);

  /* General settings. */
  printf("\n  BOX_SIZE        = " OFMT_DBL, conf->Lbox);
  printf("\n  NUM_GRID        = %d", conf->Ngrid);
  printf("\n  NUM_TRACER      = %ld", conf->Ngal);

  /* Initial condition. */
  printf("\n  LINEAR_PK       = %s", conf->plin);
  if (conf->bao_mod != 0)
    printf("\n  LINEAR_PK_NW    = %s", conf->pnw);
  if (conf->eval_growth)
    printf("\n  REDSHIFT_PK     = " OFMT_DBL, conf->zpk);
  printf("\n  PK_INTERP_LOG   = %c", conf->logint ? 'T' : 'F');
  const char *rng_name[2] = {"MRG32K3A", "MT19937"};
  printf("\n  RAND_GENERATOR  = %d (%s)", conf->rng, rng_name[conf->rng]);
  printf("\n  RAND_SEED       = %ld", conf->seed);
  printf("\n  FIX_AMPLITUDE   = %c", conf->fixamp ? 'T' : 'F');
  printf("\n  INVERT_PHASE    = %c", conf->iphase ? 'T' : 'F');

  /* Cosmological parameters. */
  if (!conf->eval_growth) {
    printf("\n  GROWTH_PK       = " OFMT_DBL, conf->growth2);
    printf("\n  VELOCITY_FAC    = " OFMT_DBL, conf->vfac);
  }
  else {
    printf("\n  OMEGA_M         = " OFMT_DBL, conf->omega_m);
    printf("\n  OMEGA_NU        = " OFMT_DBL, conf->omega_nu);
    printf("\n  DE_EOS_W        = " OFMT_DBL, conf->eos_w);
    printf("\n  REDSHIFT        = " OFMT_DBL, conf->redshift);
  }

  /* Mock generation parameters. */
  printf("\n  BAO_ENHANCE     = " OFMT_DBL, conf->bao_mod);
  printf("\n  RHO_CRITICAL    = " OFMT_DBL, conf->rho_c);
  printf("\n  RHO_EXP         = " OFMT_DBL, conf->rho_exp);
  printf("\n  PDF_BASE        = " OFMT_DBL, conf->pdf_base);
  printf("\n  SIGMA_VELOCITY  = " OFMT_DBL, conf->sigv);
  printf("\n  ATTACH_PARTICLE = %c", conf->particle ? 'T' : 'F');

  /* Output. */
  printf("\n  OUTPUT          = %s", conf->output);
  const char *ofmt_name[2] = {"ASCII", "FITS"};
  printf("\n  OUTPUT_FORMAT   = %d (%s)", conf->ofmt, ofmt_name[conf->ofmt]);
  printf("\n  OUTPUT_HEADER   = %c", conf->header ? 'T' : 'F');
  printf("\n  OVERWRITE       = %d\n", conf->ovwrite);
#ifdef OMP
  printf("  OMP_NUM_THREADS = %d\n", conf->nthread);
#endif
}


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
CONF *load_conf(const int argc, char *const *argv) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose) conf_print(conf);

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  printf(FMT_DONE);
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  if (conf->plin) free(conf->plin);
  if (conf->pnw) free(conf->pnw);
  if (conf->output) free(conf->output);
  free(conf);
}

