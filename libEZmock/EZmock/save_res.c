/*******************************************************************************
* save_res.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023-2025 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]
 
*******************************************************************************/

#include "EZmock.h"
#include "structs.h"
#include "config.h"
#include "errmsg.h"
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

/* Shortcut for writing a line to the file. */
#define WRITE_LINE(...)                                         \
  if (write_ascii_line(__VA_ARGS__)) {                          \
    free(chunk); fclose(fp); return EZMOCK_ERR_FILE;            \
  }

/*============================================================================*\
                    Function for writing to an ASCII file
\*============================================================================*/

/******************************************************************************
Function `write_ascii_line`:
  Write a line to the buffer and save it to the file if necessary.
Arguments:
  * `fp`:       pointer to the file to be written to;
  * `chunk`:    pointer to the buffer for writing;
  * `size`:     occupied size of the buffer;
  * `format`:   a string specifying how the data is interpreted;
  * `...`:      arguments specifying data to be saved.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int write_ascii_line(FILE *fp, char *chunk, size_t *size,
    const char *restrict format, ...) {
  if (!fp) return EZMOCK_ERR_FILE;
  if (!format || !(*format)) return EZMOCK_ERR_ARG_NULL;

  /* Get the size of the line to be written. */
  va_list args;
  va_start(args, format);
  int n = vsnprintf(chunk + *size, EZMOCK_FILE_CHUNK - *size, format, args);
  va_end(args);
  if (n < 0 || n >= EZMOCK_FILE_CHUNK) return EZMOCK_ERR_FILE;

  if (*size + n >= EZMOCK_FILE_CHUNK) {
    if (fwrite(chunk, (*size) * sizeof(char), 1, fp) != 1) return EZMOCK_ERR_FILE;
    *size = 0;

    va_start(args, format);
    n = vsnprintf(chunk, EZMOCK_FILE_CHUNK, format, args);
    va_end(args);
    if (n >= EZMOCK_FILE_CHUNK) return EZMOCK_ERR_FILE;
  }

  *size += n;
  return 0;
}

/*============================================================================*\
                  Interface for saving the tracer catalogue
\*============================================================================*/

/******************************************************************************
Function `EZmock_write_ascii`:
  Write the tracer catalogue to an ASCII file.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `x`:        array for the x coordinates;
  * `y`:        array for the y coordinates;
  * `z`:        array for the z coordinates;
  * `vx`:       array for the peculiar velocities along the x direction;
  * `vy`:       array for the peculiar velocities along the y direction;
  * `vz`:       array for the peculiar velocities along the z direction;
  * `ntracer`:  number of tracers to be written;
  * `rsd_fac`:  factor of the z-velocity to be added to the z coordinate for
                redshift-space distortion; write the redshift-space coordinates
                if `rsd_fac` > 0, otherwise write the real-space coordinates
                and velocities;
  * `header`:   indicate if the header is to be written;
  * `fname`:    name of the output file;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int EZmock_write_ascii(EZMOCK *ez, const real *x, const real *y, const real *z,
    const real *vx, const real *vy, const real *vz, const size_t ntracer,
    const double rsd_fac, const bool header, const char *fname, int *err) {
  /* Validate arguments. */
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return (*err = EZMOCK_ERR_ARG_EZ);
  if (!x || !(*x) || !y || !(*y) || !z || !(*z) ||
      !vx || !(*vx) || !vy || !(*vy) || !vz || !(*vz)) {
    return (*err = EZMOCK_ERR_ARG_NULL);
  }
  if (!ntracer) return (*err = EZMOCK_ERR_ARG_NDATA);
  if (rsd_fac > ((EZMOCK_CONF *) ez->conf)->Lbox)
    return (*err = EZMOCK_ERR_ARG_RSD);
  if (!fname || !(*fname)) return (*err = EZMOCK_ERR_ARG_NULL);

  /* Allocate buffer memory for file writing. */
  char *chunk = calloc(EZMOCK_FILE_CHUNK, sizeof(char));
  if (!chunk) return (*err = EZMOCK_ERR_MEMORY);
  size_t size = 0;

  /* Open the file for writing. */
  FILE *fp = fopen(fname, "w");
  if (!fp) {
    free(chunk);
    return (*err = EZMOCK_ERR_FILE);
  }

  /* Write the header. */
  if (header) {
    EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
    EZMOCK_RNG *rng = (EZMOCK_RNG *) ez->rng;
    EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
    EZMOCK_PAR *par = (EZMOCK_PAR *) ez->par;
    EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

    WRITE_LINE(fp, chunk, &size, "%c BOX_SIZE=" OFMT_DBL " , NUM_GRID=%d , "
        "FIX_AMPLITUDE=%c , INVERT_PHASE=%c\n",
        EZMOCK_SAVE_COMMENT, conf->Lbox, conf->Ng,
        mesh->fixamp ? 'T' : 'F', mesh->iphase ? 'T' : 'F');
    WRITE_LINE(fp, chunk, &size, "%c GROWTH_PK=" OFMT_DBL " , VELOCITY_FAC="
        OFMT_DBL "\n",
        EZMOCK_SAVE_COMMENT, cosmo->growth2, cosmo->vfac);
    const char *rng_name[2] = {"MRG32K3A", "MT19937"};
    WRITE_LINE(fp, chunk, &size, "%c BAO_ENHANCE=" OFMT_DBL " , RHO_CRITICAL="
        OFMT_DBL " , RHO_EXP=" OFMT_DBL " , PDF_BASE=" OFMT_DBL " , "
        "SIGMA_VELOCITY=" OFMT_DBL "\n%c RAND_GENERATOR=%s , RAND_SEED=%ld\n",
        EZMOCK_SAVE_COMMENT, par->bao_enhance, par->rho_c, par->rho_exp,
        par->pdf_base, par->sigma_v, EZMOCK_SAVE_COMMENT,
        rng_name[rng->ran], rng->seed);
  }

  if (rsd_fac > 0) {   /* save the redshift-space coordinates */
    WRITE_LINE(fp, chunk, &size, "%c x(1) y(2) z_rsd(3)\n",
        EZMOCK_SAVE_COMMENT);

    EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
    for (size_t i = 0; i < ntracer; i++) {
      double z_rsd = z[i] + rsd_fac * vz[i];
      if (z_rsd >= conf->Lbox) z_rsd -= conf->Lbox;
      if (z_rsd < 0) z_rsd += conf->Lbox;
      WRITE_LINE(fp, chunk, &size, REAL_OFMT " " REAL_OFMT " "
          REAL_OFMT "\n", x[i], y[i], z_rsd);
    }
  }
  else {               /* save the real-space coordinates and velocities */
    WRITE_LINE(fp, chunk, &size, "%c x(1) y(2) z(3) vx(4) vy(5) vz(6)\n",
        EZMOCK_SAVE_COMMENT);
    for (size_t i = 0; i < ntracer; i++) {
      WRITE_LINE(fp, chunk, &size, REAL_OFMT " " REAL_OFMT " "
          REAL_OFMT " " REAL_OFMT " " REAL_OFMT " " REAL_OFMT "\n",
          x[i], y[i], z[i], vx[i], vy[i], vz[i]);
    }
  }

  /* Write the buffer to the file. */
  if (fwrite(chunk, size * sizeof(char), 1, fp) != 1) {
    free(chunk); fclose(fp);
    return (*err = EZMOCK_ERR_FILE);
  }

  free(chunk); fclose(fp);
  return (*err = EZMOCK_SUCCESS);
}