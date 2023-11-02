/*******************************************************************************
* save_res.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "save_res.h"
#include "write_file.h"
#include "structs.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

#ifdef WITH_CFITSIO
#include <fitsio.h>
#endif

/* Shortcut for writing a line to the file. */
#define WRITE_LINE(...)                                         \
  if (output_writeline(__VA_ARGS__)) {                          \
    output_destroy(ofile); return EZMOCK_ERR_FILE;              \
  }

/* Shortcut for swapping two numbers. */
#define SWAP(x,y,t)                                             \
  (t) = (x); (x) = (y); (y) = (t);

/* Shortcut for cfitsio error handling. */
#define FITS_ABORT {                                            \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (fp) fits_close_file(fp, &status);                         \
  return EZMOCK_ERR_FILE;                                       \
}

/*============================================================================*\
                    Functions for saving the tracer catalogue
\*============================================================================*/

/******************************************************************************
Function `save_box_ascii`:
  Write the periodic tracer catalogue to an ASCII file.
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
static int save_box_ascii(CONF *conf, const real *x, const real *y,
    const real *z, const real *vx, const real *vy, const real *vz,
    const size_t ndata) {
  /* Initialise the interface for writing the file. */
  OFILE *ofile = output_init();
  if (!ofile) return EZMOCK_ERR_FILE;
  if (output_newfile(ofile, conf->output)) {
    output_destroy(ofile);
    return EZMOCK_ERR_FILE;
  }

  /* Write header. */
  if (conf->header) {
    WRITE_LINE(ofile, "%c BOX_SIZE=" OFMT_DBL " , NUM_GRID=%d , "
        "FIX_AMPLITUDE=%c , INVERT_PHASE=%c\n",
        EZMOCK_SAVE_COMMENT, conf->Lbox, conf->Ngrid,
        conf->fixamp ? 'T' : 'F', conf->iphase ? 'T' : 'F');
    WRITE_LINE(ofile, "%c GROWTH_PK=" OFMT_DBL " , VELOCITY_FAC=" OFMT_DBL "\n",
        EZMOCK_SAVE_COMMENT, conf->growth2, conf->vfac);
    const char *rng_name[2] = {"MRG32K3A", "MT19937"};
    WRITE_LINE(ofile, "%c BAO_ENHANCE=" OFMT_DBL " , RHO_CRITICAL=" OFMT_DBL
        " , RHO_EXP=" OFMT_DBL " , PDF_BASE=" OFMT_DBL " , SIGMA_VELOCITY="
        OFMT_DBL "\n%c RAND_GENERATOR=%s , RAND_SEED=%ld\n",
        EZMOCK_SAVE_COMMENT, conf->bao_mod, conf->rho_c, conf->rho_exp,
        conf->pdf_base, conf->sigv, EZMOCK_SAVE_COMMENT,
        rng_name[conf->rng], conf->seed);
    WRITE_LINE(ofile, "%c x(1) y(2) z(3) vx(4) vy(5) vz(6)\n",
        EZMOCK_SAVE_COMMENT);
  }

  /* Write the catalogue. */
  for (size_t i = 0; i < ndata; i++) {
    WRITE_LINE(ofile, REAL_OFMT " " REAL_OFMT " " REAL_OFMT " " REAL_OFMT " "
        REAL_OFMT " " REAL_OFMT "\n", x[i], y[i], z[i], vx[i], vy[i], vz[i]);
  }

  output_destroy(ofile);
  return 0;
}

/******************************************************************************
Function `save_cutsky_ascii`:
  Write the cut-sky tracer catalogue to an ASCII file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     instance of the cut-sky catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int save_cutsky_ascii(const CONF *conf, const CDATA *data) {
  /* Initialise the interface for writing the file. */
  OFILE *ofile = output_init();
  if (!ofile) return EZMOCK_ERR_FILE;
  if (output_newfile(ofile, conf->cutout)) {
    output_destroy(ofile);
    return EZMOCK_ERR_FILE;
  }

  /* Write header. */
  if (conf->header) {
    WRITE_LINE(ofile, "%c BOX_SIZE=" OFMT_DBL " , NUM_GRID=%d , "
        "FIX_AMPLITUDE=%c , INVERT_PHASE=%c\n",
        EZMOCK_SAVE_COMMENT, conf->Lbox, conf->Ngrid,
        conf->fixamp ? 'T' : 'F', conf->iphase ? 'T' : 'F');
    WRITE_LINE(ofile, "%c GROWTH_PK=" OFMT_DBL " , VELOCITY_FAC=" OFMT_DBL "\n"
        "%c OMEGA_M=" OFMT_DBL " , OMEGA_NU=" OFMT_DBL " , DE_EOS_W=" OFMT_DBL
        "\n",
        EZMOCK_SAVE_COMMENT, conf->growth2, conf->vfac,
        EZMOCK_SAVE_COMMENT, conf->omega_m, conf->omega_nu, conf->eos_w);
    const char *rng_name[2] = {"MRG32K3A", "MT19937"};
    WRITE_LINE(ofile, "%c BAO_ENHANCE=" OFMT_DBL " , RHO_CRITICAL=" OFMT_DBL
        " , RHO_EXP=" OFMT_DBL " , PDF_BASE=" OFMT_DBL " , SIGMA_VELOCITY="
        OFMT_DBL "\n%c RAND_GENERATOR=%s , RAND_SEED=%ld\n",
        EZMOCK_SAVE_COMMENT, conf->bao_mod, conf->rho_c, conf->rho_exp,
        conf->pdf_base, conf->sigv, EZMOCK_SAVE_COMMENT,
        rng_name[conf->rng], conf->seed);
    WRITE_LINE(ofile, "%c REDSHIFT=" OFMT_DBL " , REDSHIFT_MIN=" OFMT_DBL
        " , REDSHIFT_MAX=" OFMT_DBL " , GALACTIC_CAP=%cGC\n",
        EZMOCK_SAVE_COMMENT, conf->redshift, conf->zmin, conf->zmax,
        conf->gcap);
    if (conf->foot) {
      WRITE_LINE(ofile, "%c RA(1) DEC(2) Z(3) Z_COSMO(4) RAN_NUM_0_1(5) "
          "STATUS(6)\n", EZMOCK_SAVE_COMMENT);
    }
    else {
      WRITE_LINE(ofile, "%c RA(1) DEC(2) Z(3) Z_COSMO(4) RAN_NUM_0_1(5)\n",
          EZMOCK_SAVE_COMMENT);
    }
  }

  /* Write the catalogue. */
  if (conf->foot) {
    for (size_t i = 0; i < data->n; i++) {
      WRITE_LINE(ofile, REAL_OFMT " " REAL_OFMT " " REAL_OFMT " " REAL_OFMT " "
          REAL_OFMT " %d\n", data->x[0][i], data->x[1][i], data->x[2][i],
          data->x[3][i], data->rand[i], data->status[i])
    }
  }
  else {
    for (size_t i = 0; i < data->n; i++) {
      WRITE_LINE(ofile, REAL_OFMT " " REAL_OFMT " " REAL_OFMT " " REAL_OFMT " "
          REAL_OFMT "\n", data->x[0][i], data->x[1][i], data->x[2][i],
          data->x[3][i], data->rand[i]);
    }
  }

  output_destroy(ofile);
  return 0;
}

#ifdef WITH_CFITSIO
/******************************************************************************
Function `save_chunk_fits`:
  Write (a chunk of) the data catalogue to a FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `fname`:    name of the file storing the data chunk;
  * `x`:        array for the x coordinates;
  * `y`:        array for the y coordinates;
  * `z`:        array for the z coordinates;
  * `vx`:       array for the peculiar velocities along the x direction;
  * `vy`:       array for the peculiar velocities along the y direction;
  * `vz`:       array for the peculiar velocities along the z direction;
  * `ndata`:    number of tracers in the chunk;
  * `ntot`:     number of tracers in total;
  * `chk_id`:   index of the chunk;
  * `nchk`:     number of chunks in total.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int save_chunk_fits(CONF *conf, const char *fname, real *x, real *y,
    real *z, real *vx, real *vy, real *vz, const size_t ndata,
    const size_t ntot, int chk_id, int nchk) {
  fitsfile *fp = NULL;
  int status = 0;

  /* Create a FITS file. */
  if (fits_create_file(&fp, fname, &status)) FITS_ABORT;

  /* Create an empty table. */
  char *names[] = {"x", "y", "z", "vx", "vy", "vz"};
  char *units[] = {"Mpc/h", "Mpc/h", "Mpc/h", "km/s", "km/s", "km/s"};
  char *cfmt[6];
  int dtypes[6];
  switch (sizeof(real)) {
    case 4:
      for (int i = 0; i < 6; i++) {
        cfmt[i] = "1E";
        dtypes[i] = TFLOAT;
      }
      break;
    case 8:
      for (int i = 0; i < 6; i++) {
        cfmt[i] = "1D";
        dtypes[i] = TDOUBLE;
      }
      break;
    default:
      P_ERR("invalid data type for floating point numbers\n");
      return EZMOCK_ERR_UNKNOWN;
  }

  if (fits_create_tbl(fp, BINARY_TBL, 0, 6, names, cfmt, units, NULL, &status))
    FITS_ABORT;

  /* Write the header unit. */
  if (conf->header) {
    size_t ntracer = ntot;
    char *rng_name[2] = {"MRG32K3A", "MT19937"};
    if (fits_write_key(fp, TLONGLONG, "NUM_TOT" , &ntracer,
            "total number of tracers"                           , &status) ||
        fits_write_key(fp, TINT     , "CHUNK_ID", &chk_id,
            "index of the chunk of the catalog (from 1)"        , &status) ||
        fits_write_key(fp, TINT     , "NCHUNK"  , &nchk,
            "total number of chunks"                            , &status) ||
        fits_write_key(fp, TDOUBLE  , "BOX_SIZE", &conf->Lbox,
            "side length of the periodic box"                   , &status) ||
        fits_write_key(fp, TINT     , "NUM_GRID", &conf->Ngrid,
            "number of grid cells per box size"                 , &status) ||
        fits_write_key(fp, TLOGICAL , "FIX_AMP" , &conf->fixamp,
            "fix amplitude of initial white noise"              , &status) ||
        fits_write_key(fp, TLOGICAL , "INV_PH"  , &conf->iphase,
            "invert phase of initial white noise"               , &status) ||
        fits_write_key(fp, TDOUBLE  , "GROW_PK" , &conf->growth2,
            "normalization factor of the power spectrum"        , &status) ||
        fits_write_key(fp, TDOUBLE  , "VEL_FAC" , &conf->vfac, "factor for "
            "computing peculiar velocity from displacement"     , &status) ||
        fits_write_key(fp, TDOUBLE  , "BAO_MOD" , &conf->bao_mod,
            "BAO enhancement parameter"                         , &status) ||
        fits_write_key(fp, TDOUBLE  , "RHO_CUT" , &conf->rho_c,
            "critical density parameter"                        , &status) ||
        fits_write_key(fp, TDOUBLE  , "RHO_EXP" , &conf->rho_exp,
            "expotential cut-off of the bias model"             , &status) ||
        fits_write_key(fp, TDOUBLE  , "PDF_BASE", &conf->pdf_base,
            "base of the power-law tracer PDF"                  , &status) ||
        fits_write_key(fp, TDOUBLE  , "SIGMA_V" , &conf->sigv,
            "standard devition of random local motion"          , &status) ||
        fits_write_key(fp, TSTRING  , "RAN_GEN" , rng_name[conf->rng],
            "random number generator"                           , &status) ||
        fits_write_key(fp, TLONGLONG, "RAN_SEED", &conf->seed,
            "random seed"                                       , &status))
      FITS_ABORT;
  }

  /* Write the tracer catalogue. */
  int cnums[6] = {1, 2, 3, 4, 5, 6};
  real *nulval[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
  real *data[6];
  data[0] = x;
  data[1] = y;
  data[2] = z;
  data[3] = vx;
  data[4] = vy;
  data[5] = vz;

  if (fits_write_cols(fp, 6, dtypes, cnums, 1, ndata, (void **) data,
      (void **) nulval, &status)) FITS_ABORT;

  if (fits_close_file(fp, &status)) {
    P_ERR("cfitsio error: ");
    fits_report_error(stderr, status);
    return EZMOCK_ERR_FILE;
  }

  return 0;
}

/******************************************************************************
Function `save_box_fits`:
  Split the periodic tracer catalogue into chunks, and save them as FITS files.
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
static int save_box_fits(CONF *conf, real *x, real *y, real *z,
    real *vx, real *vy, real *vz, const size_t ndata) {
  const int nside = EZMOCK_FITS_CHUNK_PER_SIDE;

  /* Save the catalogue as a whole if no chunk is needed. */
  if (nside <= 1) {
    return save_chunk_fits(conf, conf->output, x, y, z, vx, vy, vz, ndata,
        ndata, 1, 1);
  }

  /* Allocate memory for counting tracers per chunk. */
  const int nchunk = nside * nside * nside;
  size_t *cnt = calloc(nchunk * 2 + 1, sizeof(size_t));
  if (!cnt) {
    P_ERR("failed to allocate memory for splitting the catalog into chunks\n");
    return EZMOCK_ERR_MEMORY;
  }
  size_t *idx = cnt + nchunk;

  /* First round: count the number of tracers per chunk. */
  const real fac = nside / conf->Lbox;
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
  for (size_t i = 0; i < ndata; i++) {
    int ix = x[i] * fac;
    int iy = y[i] * fac;
    int iz = z[i] * fac;
    if (ix >= nside) ix = nside - 1;
    if (iy >= nside) iy = nside - 1;
    if (iz >= nside) iz = nside - 1;

    int ichunk = (ix * nside + iy) * nside + iz;
#ifdef OMP
#pragma omp atomic
#endif
    cnt[ichunk]++;
  }

  /* Compute the starting and ending indices of each chunk. */
  for (int i = 1; i <= nchunk; i++) idx[i] = idx[i - 1] + cnt[i - 1];
  if (idx[nchunk] != ndata) {
    P_ERR("failed to generate chunks of the tracer catalog\n");
    free(cnt);
    return EZMOCK_ERR_UNKNOWN;
  }
  memset(cnt, 0, nchunk * sizeof(size_t));

  /* Second round: reorder the catalogue by chunk id. */
  int idst = -1;
  for (int isrc = 0; isrc < nchunk; isrc++) {
    for (size_t i = idx[isrc] + cnt[isrc]; i < idx[isrc + 1]; i++) {
      if (idst < 0) {
        int ix = x[i] * fac;
        int iy = y[i] * fac;
        int iz = z[i] * fac;
        if (ix >= nside) ix = nside - 1;
        if (iy >= nside) iy = nside - 1;
        if (iz >= nside) iz = nside - 1;
        idst = (ix * nside + iy) * nside + iz;
      }
      if (idst == isrc) {       /* the tracer is already in the right chunk */
        cnt[idst]++;
        idst = -1;
        continue;
      }

      /* Find a tracer in the destination chunk for swapping. */
      size_t j;
      int jdst = 0;
      for (j = idx[idst] + cnt[idst]; j < idx[idst + 1]; j++) {
        int jx = x[j] * fac;
        int jy = y[j] * fac;
        int jz = z[j] * fac;
        if (jx >= nside) jx = nside - 1;
        if (jy >= nside) jy = nside - 1;
        if (jz >= nside) jz = nside - 1;
        jdst = (jx * nside + jy) * nside + jz;

        /* Skip destination tracers that are in the right chunk. */
        if (jdst == idst) cnt[idst]++;
        else break;
      }
      if (j == idx[idst + 1]) {
        P_ERR("failed to split the catalog into chunks\n");
        free(cnt);
        return EZMOCK_ERR_UNKNOWN;
      }

      /* Swap the tracers. */
      real tmp;
      SWAP(x[i], x[j], tmp);
      SWAP(y[i], y[j], tmp);
      SWAP(z[i], z[j], tmp);
      SWAP(vx[i], vx[j], tmp);
      SWAP(vy[i], vy[j], tmp);
      SWAP(vz[i], vz[j], tmp);
      cnt[idst]++;
      idst = jdst;
      i--;
    }
  }

#ifdef EZMOCK_DEBUG
  for (int isrc = 0; isrc < nchunk; isrc++) {
    for (size_t i = idx[isrc] + cnt[isrc]; i < idx[isrc + 1]; i++) {
      int ix = x[i] * fac;
      int iy = y[i] * fac;
      int iz = z[i] * fac;
      if (ix >= nside) ix = nside - 1;
      if (iy >= nside) iy = nside - 1;
      if (iz >= nside) iz = nside - 1;
      int idst = (ix * nside + iy) * nside + iz;
      if (idst != isrc) {
        P_ERR("not all tracers are matched\n");
        free(cnt); return EZMOCK_ERR_UNKNOWN;
      }
    }
  }
#endif

  /* Generate filenames of chunks. */
  size_t len = strlen(conf->output);
  size_t elen = strlen(".fits");
  if (len < elen || strncmp(conf->output + len - elen, ".fits", elen)) {
    elen = strlen(".fits.gz");
    if (len < elen || strncmp(conf->output + len - elen, ".fits.gz", elen)) {
      P_ERR("invalid output file name: %s\n", conf->output);
      return EZMOCK_ERR_CFG;
    }
  }

  char *fname = malloc(len + 6);        /* add digits and null termination */
  if (!fname) {
    P_ERR("failed to allocate memory for chunks of the catalog\n");
    free(cnt);
    return EZMOCK_ERR_MEMORY;
  }

  /* Set the basename with symbol for force writing. */
  fname[0] = '!';
  memcpy(fname + 1, conf->output, len - elen);

  for (int i = 0; i < nchunk; i++) {
    /* Generate the filename of this chunk. */
    int n = snprintf(fname + 1 + len - elen, 5, ".%d", i);
    if (n <= 0) {
      P_ERR("failed to generate filenames for chunks of the output catalog\n");
      free(fname); free(cnt);
      return EZMOCK_ERR_UNKNOWN;
    }
    strncpy(fname + 1 + len - elen + n, conf->output + len - elen, elen + 1);

    /* Save chunk of the data catalogue. */
    int e = save_chunk_fits(conf, fname, x + idx[i], y + idx[i], z + idx[i],
        vx + idx[i], vy + idx[i], vz + idx[i], cnt[i], ndata, i + 1, nchunk);
    if (e) return e;
  }

  free(fname);
  free(cnt);
  return 0;
}

/******************************************************************************
Function `save_cutsky_fits`:
  Write the cut-sky tracer catalogue to a FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     instance of the cut-sky catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int save_cutsky_fits(CONF *conf, CDATA *data) {
  /* Allocate memory for the filename. */
  const size_t len = strlen(conf->cutout) + 1;
  char *fname = malloc(len + 1);
  if (!fname) {
    P_ERR("failed to allocate memory for saving the catalog\n");
    return EZMOCK_ERR_MEMORY;
  }
  fname[0] = '!';
  strncpy(fname + 1, conf->cutout, len);

  /* Create the FITS file. */
  fitsfile *fp = NULL;
  int status = 0;
  if (fits_create_file(&fp, fname, &status)) {
    free(fname);
    FITS_ABORT;
  }
  free(fname);

  /* Create an empty table. */
  char *names[] = {"RA", "DEC", "Z", "Z_COSMO", "RAN_NUM_0_1", "STATUS"};
  char *units[] = {"deg", "deg", NULL, NULL, NULL, NULL};
  char *cfmt[6];
  int dtypes[6];
  switch (sizeof(real)) {
    case 4:
      for (int i = 0; i < 5; i++) {
        cfmt[i] = "1E";
        dtypes[i] = TFLOAT;
      }
      break;
    case 8:
      for (int i = 0; i < 5; i++) {
        cfmt[i] = "1D";
        dtypes[i] = TDOUBLE;
      }
      break;
    default:
      P_ERR("invalid data type for floating point numbers\n");
      return EZMOCK_ERR_UNKNOWN;
  }
  switch (sizeof *(data->status)) {
    case 1:
      cfmt[5] = "1B";
      dtypes[5] = TBYTE;
      break;
    case 2:
      cfmt[5] = "1I";
      dtypes[5] = TUSHORT;
      break;
    case 4:
      cfmt[5] = "1J";
      dtypes[5] = TINT32BIT;
      break;
    case 8:
      cfmt[5] = "1K";
      dtypes[5] = TULONGLONG;
      break;
    default:
      P_ERR("invalid data type for the bit-code of cut-sky catalog\n");
      return EZMOCK_ERR_UNKNOWN;
  }

  const int ncol = (conf->foot) ? 6 : 5;
  if (fits_create_tbl(fp, BINARY_TBL, 0, ncol, names, cfmt, units, NULL,
      &status)) FITS_ABORT;

  /* Write the header unit. */
  if (conf->header) {
    char *rng_name[2] = {"MRG32K3A", "MT19937"};
    char gcap[] = "?GC";
    if (conf->gcap == 'N' || conf->gcap == 'n') gcap[0] = 'N';
    else gcap[0] = 'S';

    if (fits_write_key(fp, TDOUBLE  , "BOX_SIZE", &conf->Lbox,
            "side length of the periodic box"                   , &status) ||
        fits_write_key(fp, TINT     , "NUM_GRID", &conf->Ngrid,
            "number of grid cells per box size"                 , &status) ||
        fits_write_key(fp, TLOGICAL , "FIX_AMP" , &conf->fixamp,
            "fix amplitude of initial white noise"              , &status) ||
        fits_write_key(fp, TLOGICAL , "INV_PH"  , &conf->iphase,
            "invert phase of initial white noise"               , &status) ||
        fits_write_key(fp, TDOUBLE  , "GROW_PK" , &conf->growth2,
            "normalization factor of the power spectrum"        , &status) ||
        fits_write_key(fp, TDOUBLE  , "VEL_FAC" , &conf->vfac,
            "ratio of peculiar velocity to displacement"        , &status) ||
        fits_write_key(fp, TDOUBLE  , "OMEGA_M" , &conf->omega_m,
            "matter (w/o neutrino) density parameter"           , &status) ||
        fits_write_key(fp, TDOUBLE  , "OMEGA_NU", &conf->omega_nu,
            "neutrino density parameter"                        , &status) ||
        fits_write_key(fp, TDOUBLE  , "DE_W"    , &conf->eos_w,
            "dark energy equation of state"                     , &status) ||
        fits_write_key(fp, TDOUBLE  , "BAO_MOD" , &conf->bao_mod,
            "BAO enhancement parameter"                         , &status) ||
        fits_write_key(fp, TDOUBLE  , "RHO_CUT" , &conf->rho_c,
            "critical density parameter"                        , &status) ||
        fits_write_key(fp, TDOUBLE  , "RHO_EXP" , &conf->rho_exp,
            "expotential cut-off of the bias model"             , &status) ||
        fits_write_key(fp, TDOUBLE  , "PDF_BASE", &conf->pdf_base,
            "base of the power-law tracer PDF"                  , &status) ||
        fits_write_key(fp, TDOUBLE  , "SIGMA_V" , &conf->sigv,
            "standard devition of random local motion"          , &status) ||
        fits_write_key(fp, TSTRING  , "RAN_GEN" , rng_name[conf->rng],
            "random number generator"                           , &status) ||
        fits_write_key(fp, TLONGLONG, "RAN_SEED", &conf->seed,
            "random seed"                                       , &status) ||
        fits_write_key(fp, TDOUBLE  , "REDSHIFT", &conf->redshift,
            "redshift snapshot of the catalog"                  , &status) ||
        fits_write_key(fp, TDOUBLE  , "Z_MIN"   , &conf->zmin,
            "minimum redshift of cut-sky catalog"               , &status) ||
        fits_write_key(fp, TDOUBLE  , "Z_MAX"   , &conf->zmax,
            "maximum redshift of cut-sky catalog"               , &status) ||
        fits_write_key(fp, TDOUBLE  , "NUM_DENS", &data->num_dens,
            "number density of the tracers"                     , &status) ||
        fits_write_key(fp, TSTRING  , "GAL_CAP" , gcap,
            "galactic cap"                                      , &status))
      FITS_ABORT;
  }

  /* Write the cut-sky catalog. */
  int cnums[6] = {1, 2, 3, 4, 5, 6};
  void *nulval[6] = {NULL, NULL, NULL, NULL, NULL, NULL};
  void *dat[6];
  for (int i = 0; i < 4; i++) dat[i] = data->x[i];
  dat[4] = data->rand;
  dat[5] = data->status;

  if (fits_write_cols(fp, ncol, dtypes, cnums, 1, data->n, dat, nulval,
      &status)) FITS_ABORT;

  if (fits_close_file(fp, &status)) {
    P_ERR("cfitsio error: ");
    fits_report_error(stderr, status);
    return EZMOCK_ERR_FILE;
  }

  return 0;
}
#endif


/*============================================================================*\
                          Interface for catalog saving
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
    real *vx, real *vy, real *vz, const size_t ndata) {
  printf("Saving the cubic mock catalog ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return EZMOCK_ERR_INIT;
  }
  if (!x || !y || !z || !vx || !vy || !vz || ndata == 0) {
    P_ERR("the tracer catalogue is not generated\n");
    return EZMOCK_ERR_INIT;
  }
  if (conf->verbose) printf("\n  Filename: `%s'\n", conf->output);
  fflush(stdout);

#ifdef WITH_CFITSIO
  switch (conf->ofmt) {
    case EZMOCK_OFMT_ASCII:
      if (save_box_ascii(conf, x, y, z, vx, vy, vz, ndata))
        return EZMOCK_ERR_SAVE;
      break;
    case EZMOCK_OFMT_FITS:
      if (save_box_fits(conf, x, y, z, vx, vy, vz, ndata))
        return EZMOCK_ERR_SAVE;
      break;
    default:
      P_ERR("invalid output file format: %d\n", conf->ofmt);
      return EZMOCK_ERR_CFG;
  }
#else
  if (save_box_ascii(conf, x, y, z, vx, vy, vz, ndata)) return EZMOCK_ERR_SAVE;
#endif

  printf(FMT_DONE);
  return 0;
}

/******************************************************************************
Function `save_cutsky`:
  Write the cut-sky tracer catalogue to a file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     instance of the cut-sky catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_cutsky(CONF *conf, CDATA *data) {
  printf("Saving the cutsky catalog ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return EZMOCK_ERR_INIT;
  }
  if (!data) {
    P_ERR("the cut-sky catalog is not generated\n");
    return EZMOCK_ERR_INIT;
  }
  if (conf->verbose) printf("\n  Filename: `%s'\n", conf->cutout);
  fflush(stdout);

#ifdef WITH_CFITSIO
  switch (conf->ofmt) {
    case EZMOCK_OFMT_ASCII:
      if (save_cutsky_ascii(conf, data)) return EZMOCK_ERR_SAVE;
      break;
    case EZMOCK_OFMT_FITS:
      if (save_cutsky_fits(conf, data)) return EZMOCK_ERR_SAVE;
      break;
    default:
      P_ERR("invalid output file format: %d\n", conf->ofmt);
      return EZMOCK_ERR_CFG;
  }
#else
  if (save_cutsky_ascii(conf, data)) return EZMOCK_ERR_SAVE;
#endif

  printf(FMT_DONE);
  return 0;
}
