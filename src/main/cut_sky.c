/*******************************************************************************
* cut_sky.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "define.h"
#include "mangle.h"
#include "prand.h"
#include "cut_sky.h"

#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                 Functions for processing the cut-sky catalogue
\*============================================================================*/

/******************************************************************************
Function `cutsky_init`:
  Initialise the cut-sky catalogue;
Arguments:
  * `ndata`:    number of tracers in the cubic mock.
Return:
  Instance of the cut-sky catalogue on success; NULL on error.
******************************************************************************/
static CDATA *cutsky_init(const size_t ndata) {
  CDATA *data = malloc(sizeof(CDATA));
  if (!data) {
    P_ERR("failed to allocate memory for the cut-sky catalog\n");
    return NULL;
  }

  for (int i = 0; i < 4; i++) data->x[i] = NULL;
  data->rand = NULL;
  data->status = NULL;
  data->n = 0;
  data->cap = ndata;    /* initialise with the number in the cubic catalogue */

  if (ndata) {
    for (int i = 0; i < 4; i++) {
      if (!(data->x[i] = malloc(ndata * sizeof(real)))) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        for (int j = 0; j < i; j++) free(data->x[j]);
        free(data);
        return NULL;
      }
    }
  }

  return data;
}

/******************************************************************************
Function `cutsky_append`:
  Append an object to the cut-sky catalogue.
Arguments:
  * `data`:     the cut-sky catalogue;
  * `ra`:       the right-acension of the tracer;
  * `dec`:      the declination of the tracer;
  * `z_cosmo`:  the real-space redshift of the tracer;
  * `z`:        the redshift-space redshift of the tracer.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cutsky_append(CDATA *data, const real ra, const real dec,
    const real z_cosmo, const real z) {
  /* Enlarge the catalogue if necessary. */
  if (data->n == data->cap) {
    if (data->cap >= (SIZE_MAX >> 1)) {
      P_ERR("unexpected size of the cut-sky catalog\n");
      return EZMOCK_ERR_CUTSKY;
    }
    data->cap <<= 1;
    for (int i = 0; i < 4; i++) {
      real *tmp = realloc(data->x[i], data->cap * sizeof(real));
      if (!tmp) {
        P_ERR("failed to allocate memory for the cut-sky catalog\n");
        return EZMOCK_ERR_MEMORY;
      }
      data->x[i] = tmp;
    }
  }

  data->x[0][data->n] = ra;
  data->x[1][data->n] = dec;
  data->x[2][data->n] = z;
  data->x[3][data->n] = z_cosmo;
  data->n++;
  return 0;
}

/******************************************************************************
Function `cutsky_run`:
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
  Instance of the cut-sky catalogue on success; NULL on error.
******************************************************************************/
static CDATA *cutsky_run(const CONF *conf, const ZCVT *cvt, real *x, real *y,
    real *z, real *vx, real *vy, real *vz, const size_t ndata) {
  /* Setup the DESI Y5 footprint. */
  int err = 0;
  MANGLE *foot = mangle_init(conf->y5foot, 0, &err);
  if (err) {
    P_ERR("%s\n", mangle_errmsg(err));
    return NULL;
  }

  /* Compute the number of box duplications needed on each side. */
  const int n_dup = ceil(sqrt(cvt->d2max) / conf->Lbox);

  /* Shift RA by 60 degree. */
  const double ra_shift = 60;
  const bool is_ngc = (conf->gcap == 'N' || conf->gcap == 'n');

#ifdef CUTSKY_DEBUG
  FILE *fp = fopen("test_cutsky.dat", "w");
  if (!fp) {
    P_ERR("error saving cut-sky catalog\n");
    exit(1);
  }
  const int Ng = 2 * n_dup;
#endif

#ifdef OMP
  /* Create thread-private cut-sky catalogue. */
  CDATA **pdata = malloc(conf->nthread * sizeof(CDATA *));
  if (!pdata) {
    P_ERR("failed to allocate memory for the cut-sky catalog\n");
    mangle_destroy(foot); return NULL;
  }
  for (int i = 0; i < conf->nthread; i++) {
    if (!(pdata[i] = cutsky_init(ndata / conf->nthread))) {
      for (int j = 0; j < i; j++) free(pdata[j]);
      free(pdata); mangle_destroy(foot);
      return NULL;
    }
  }

  const size_t pnum = ndata / conf->nthread;
  const int rem = ndata % conf->nthread;
#else
  /* Create cut-sky catalogue. */
  CDATA *data = cutsky_init(ndata);
  if (!data) {
    mangle_destroy(foot); return NULL;
  }
#endif

  /* Start processing the catalogue with box duplication. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;
    for (size_t n = istart; n < iend; n++) {
#else
    for (size_t n = 0; n < ndata; n++) {
#endif
      for (int i = -n_dup; i < n_dup; i++) {
        double xx = x[n] + i * conf->Lbox;
        for (int j = -n_dup; j < n_dup; j++) {
          double yy = y[n] + j * conf->Lbox;
          for (int k = -n_dup; k < n_dup; k++) {
            double zz = z[n] + k * conf->Lbox;

            /* Compute and trim the radial distance with tolerance for RSD. */
            double d2 = xx * xx + yy * yy + zz * zz;
            if (d2 > cvt->d2max || d2 < cvt->d2min) continue;

            /* Compute the line-of-sight velocity. */
            double dist = sqrt(d2);
            double vel = (vx[n] * xx + vy[n] * yy + vz[n] * zz) / dist;

            /* Convert squared distance to redshift. */
            double z_real = convert_z(cvt, d2);
            double z_red = z_real + vel * (1 + z_real) / SPEED_OF_LIGHT;
            if (z_red < conf->zmin || z_red > conf->zmax) continue;

            /* Compute sky coordinates. */
            double ra, dec;
            if (dist < DOUBLE_TOL) ra = dec = 0;
            else {
              dec = asin(zz / dist) * RAD_2_DEGREE;
              ra = atan2(yy, xx) * RAD_2_DEGREE + ra_shift;
              if (ra < 0) ra += 360;
            }

            /* Pre-select NGC/SGC. */
            if ((ra > DESI_NGC_RA_MIN && ra < DESI_NGC_RA_MAX) != is_ngc)
              continue;

            /* Trim survey footprint. */
            POLYGON *poly = mangle_query(foot, ra, dec);
            if (!poly) continue;

#ifdef CUTSKY_DEBUG
#ifdef OMP
#pragma omp critical
#endif
            fprintf(fp, OFMT_DBL " " OFMT_DBL " " OFMT_DBL " " OFMT_DBL " "
                OFMT_DBL " " OFMT_DBL " %d\n",
                ra, dec, z_real, x[n], y[n], z[n], (i * Ng + j) * Ng + k);
#endif

#ifdef OMP
            if (cutsky_append(pdata[tid], ra, dec, z_real, z_red)) {
              err = EZMOCK_ERR_MEMORY;
              break;
            }
#else
            if (cutsky_append(data, ra, dec, z_real, z_red)) {
              mangle_destroy(foot);
              cutsky_destroy(data);
              return NULL;
            }
#endif
          }     /* for k */
        }       /* for j */
      }         /* for i */
    }           /* for n */
#ifdef OMP
  }             /* omp parallel */
#endif

  mangle_destroy(foot);

#ifdef CUTSKY_DEBUG
  fclose(fp);
  exit(0);
#endif

#ifdef OMP
  if (err) {
    for (int i = 0; i < conf->nthread; i++) free(pdata[i]);
    free(pdata);
    return NULL;
  }

  /* Gather thread-private data. */
  CDATA *data = cutsky_init(0);
  if (!data) {
    P_ERR("failed to allocate memory for the cut-sky catalog\n");
    for (int i = 0; i < conf->nthread; i++) free(pdata[i]);
    free(pdata); return NULL;
  }
  size_t num = 0;
  for (int i = 0; i < conf->nthread; i++) num += pdata[i]->n;
  if (!num) {
    P_ERR("empty cut-sky catalog\n");
    for (int i = 0; i < conf->nthread; i++) free(pdata[i]);
    free(pdata); cutsky_destroy(data);
    return NULL;
  }
  data->n = data->cap = num;

  for (int i = 0; i < 4; i++) {
    if (!(data->x[i] = malloc(data->cap * sizeof(real)))) {
      P_ERR("failed to allocate memory for the cut-sky catalog\n");
      for (int j = 0; j < conf->nthread; j++) free(pdata[j]);
      free(pdata); cutsky_destroy(data);
      return NULL;
    }
    num = 0;
    for (int j = 0; j < conf->nthread; j++) {
      if (pdata[j]->n)
        memcpy(data->x[i] + num, pdata[j]->x[i], pdata[j]->n * sizeof(real));
      num += pdata[j]->n;
      free(pdata[j]->x[i]);
      pdata[j]->x[i] = NULL;
    }
  }

  for (int i = 0; i < conf->nthread; i++) free(pdata[i]);
  free(pdata);
#else
  if (!data->n) {
    P_ERR("empty cut-sky catalog\n");
    cutsky_destroy(data);
    return NULL;
  }
  /* Reduce memory cost if applicable. */
  if (data->n < data->cap) {
    for (int i = 0; i < 4; i++) {
      real *tmp = realloc(data->x[i], data->n * sizeof(real));
      if (tmp) data->x[i] = tmp;
    }
  }
#endif

  data->num_dens = ndata / pow(conf->Lbox, 3);
  return data;
}

/******************************************************************************
Function `cutsky_post_process`:
  Post-process the cut-sky catalogue.
Arguments:
  * `data`:     instance of the cut-sky catalogue;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int cutsky_post_process(CDATA *data, const CONF *conf) {
  /* Allocate memory. */
  if ((conf->foot &&
      !(data->status = calloc(data->n, sizeof *(data->status)))) ||
      !(data->rand = malloc(data->n * sizeof *(data->rand)))) {
    P_ERR("failed to allocate memory for the cut-sky catalog\n");
    return EZMOCK_ERR_MEMORY;
  }

  /* Mark the footprint of the current data release. */
  if (conf->foot) {
    int err = 0;
    MANGLE *foot = mangle_init(conf->foot, 0, &err);
    if (err) {
      P_ERR("%s\n", mangle_errmsg(err));
      return EZMOCK_ERR_CUTSKY;
    }

    /* Check footprint. */
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
    for (size_t i = 0; i < data->n; i++) {
      POLYGON *poly = mangle_query(foot, data->x[0][i], data->x[1][i]);
      if (poly) data->status[i] = 1;
    }

    mangle_destroy(foot);
  }
  if (conf->verbose)
    printf("  Footprint of the current data release applied\n");

  /* Generate randoms for downsampling. */
  int err = 0;
  prand_t *rng = prand_init(conf->rng, conf->seed, conf->nthread, 0, &err);
  if (PRAND_IS_ERROR(err) || PRAND_IS_WARN(err)) {
    P_ERR("failed to initialise the random number generator\n");
    return EZMOCK_ERR_UNKNOWN;
  }

#ifdef OMP
  const size_t pnum = data->n / conf->nthread;
  const int rem = data->n % conf->nthread;

  /* Distribute tracers to threads. */
#pragma omp parallel num_threads(conf->nthread)
  {
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

  /* Jump ahead random states with a base step of 2^32. */
    rng->reset(rng->state_stream[tid], conf->seed, (1ULL << 32) + istart, &err);
#pragma omp barrier             /* synchronize `err` */
    if (!(PRAND_IS_ERROR(err))) {
      for (size_t i = istart; i < iend; i++)
        data->rand[i] = rng->get_double(rng->state_stream[tid]);
    }
  }
  if (PRAND_IS_ERROR(err)) {
    P_ERR("failed to jump ahead the random state\n");
    return EZMOCK_ERR_RNG;
  }
#else
  /* Jump ahead the random state with a step of 2^32. */
  rng->reset(rng->state, conf->seed, 1ULL << 32, &err);
  if (PRAND_IS_ERROR(err)) {
    P_ERR("failed to jump ahead the random state\n");
    return EZMOCK_ERR_RNG;
  }
  for (size_t i = 0; i < data->n; i++) {
    data->rand[i] = rng->get_double(rng->state);
  }
#endif

  if (conf->verbose) printf("  Column with uniform random numbers added\n");

  prand_destroy(rng);
  return 0;
}


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
    real *vx, real *vy, real *vz, const size_t ndata) {
  printf("Constructing the cut-sky catalog ...");
  if (!conf) {
    P_ERR("configurations are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Initialise the cut-sky catalogue. */
  CDATA *data;
  if (!(data = cutsky_run(conf, cvt, x, y, z, vx, vy, vz, ndata))) {
    cutsky_destroy(data); return NULL;
  }
  if (conf->verbose)
    printf("  Cut-sky catalog constructed with %zu tracers\n", data->n);

  if (data->n == 0) {
    P_ERR("no tracer in the cut-sky catalog\n");
    cutsky_destroy(data); return NULL;
  }

  if (cutsky_post_process(data, conf)) {
    cutsky_destroy(data); return NULL;
  }

  printf(FMT_DONE);
  return data;
}

/******************************************************************************
Function `cutsky_init`:
  Deconstruct the cut-sky catalogue;
Arguments:
  * `data`:     instance of the cut-sky catalogue.
******************************************************************************/
void cutsky_destroy(CDATA *data) {
  if (!data) return;
  for (int i = 0; i < 4; i++) {
    if (data->x[i]) free(data->x[i]);
  }
  if (data->status) free(data->status);
  if (data->rand) free(data->rand);
  free(data);
}
