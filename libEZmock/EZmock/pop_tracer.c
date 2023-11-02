/*******************************************************************************
* pop_tracer.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com> [MIT license]

*******************************************************************************/

#include "EZmock.h"
#include "structs.h"
#include "config.h"
#include "errmsg.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

#ifdef OMP
#include <omp.h>
#define OMP_ATOMIC      _Pragma("omp atomic")
#else
#define OMP_ATOMIC
#endif

/*============================================================================*\
                       Definitions for template functions
\*============================================================================*/

/******************************************************************************
Function `qselect`:
  Move the nth largest/smallest elements of an array A, to the left of A[n].
Arguments:
  * `A`:        pointer to the array;
  * `n`:        the number of smallest elements to be moved;
  * `len`:      length of the array `A`;
  * `buf`:      pointer to the temporary space for swapping elements;
  * `arg`:      additional arguments for the comparison function.
******************************************************************************/

/* Data type of the array to be selected. */
#define QSELECT_DTYPE                   size_t
/* Macro for element comparison. */
#define QSELECT_COMPARE(a,b,m)                                          \
  ((real *) (m))[*(b)] - ((real *) (m))[*(a)]
/* Macro for swapping two elements. */
#define QSELECT_SWAP(a,b) {                                             \
  *((size_t *) buf) = *(a); *(a) = *(b); *(b) = *((size_t *) buf);      \
}

#include "qselect.c"

/*============================================================================*\
              Functions for sampling tracers in the density field
\*============================================================================*/

/******************************************************************************
Function `gauss_rand`:
  Generate a pair of Gaussian random numbers using the Box-Muller transform.
Arguments:
  * `rng`:      interface of the random number generator;
  * `state`:    the random state;
  * `sigma`:    standard deviation of the Gaussian distribution;
  * `rand1`:    the first random number produced;
  * `rand2`:    the second random number produced.
******************************************************************************/
static inline void gauss_rand(prand_t *rng, void *state, const double sigma,
    double *rand1, double *rand2) {
  /* Generate the modulus and angle. */
  double r = -2 * log(rng->get_double_pos(state));
  double theta = rng->get_double(state) * 2 * M_PI;

  /* Box-Muller transform. */
  *rand1 = sigma * sqrt(r) * cos(theta);
  *rand2 = sigma * sqrt(r) * sin(theta);
}

/******************************************************************************
Function `eval_pdf`:
  Evaluate the power-law probability distribution function (PDF) of tracers.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `nexp`:     expected number of tracers to be generated;
  * `pbase`:    base of the expected probability distribution function;
  * `pdf`:      instance of the EZmock tracer PDF;
  * `nmax`:     the evaluated maximum number of tracers in a cell;
  * `ncell`:    the evaluated number of cells that host tracers;
  * `ntracer`:  actual number of tracers to be generated;
  * `err`:      integer storing the error code.
******************************************************************************/
static void eval_pdf(EZMOCK *ez, const size_t nexp, const real pbase,
    EZMOCK_PDF *pdf, int *err) {
  /* Find nmax and pdf_A iteratively. */
  size_t n;
  double A, sum;
  double power = pbase;
  A = sum = 0;
  for (n = 1; n <= nexp; n++) {
    sum += power * n;
    A = nexp / sum;
    power *= pbase;
    if (A * power < 0.5) break;
  }
  if (n > EZMOCK_MAX_NUM_TRACER_PER_CELL) {
    *err = EZMOCK_ERR_RHO_PBASE; return;
  }

  size_t *num = malloc(n * sizeof(size_t));
  if (!num) {
    *err = EZMOCK_ERR_MEMORY; return;
  }

  size_t nc, nobj;      /* total number of cells and tracers, respectively */
  nc = nobj = 0;
  for (size_t i = 0; i < n; i++) {
    num[i] = round(A * pow(pbase, n - i));
    nc += num[i];
    nobj += num[i] * (n - i);
  }

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  if (nc > (size_t) conf->Ng * conf->Ng * conf->Ng || nobj == 0) {
    free(num); *err = EZMOCK_ERR_PAR_PBASE; return;
  }

  pdf->pdf = num;
  pdf->nmax = n;
  pdf->ncell = nc;
  pdf->ntot = nobj;
}

/******************************************************************************
Function `bias_model`:
  Create the tracer number density field with the effective bias model.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `pdf`:      instance of the EZmock tracer PDF;
  * `rho_c`:    the density threshold parameter rho_critical;
  * `rho_sat`:  the density saturation parameter rho_saturation;
  * `rho_exp`:  the density when the exponential bias starts to saturate;
  * `lambda`:   width of the Gaussian random distribution for stochaostic bias;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static void bias_model(EZMOCK *ez, EZMOCK_PDF *pdf, const real rho_c,
    const real rho_sat, const real rho_exp, const real lambda, int *err) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_RNG *erng = (EZMOCK_RNG *) ez->rng;

  prand_t *rng = erng->rng;
  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;

  /* Allocate memory for the tracer number density field if necessary. */
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  if (!mesh->rhot) {
    if (mesh->rhok2) mesh->rhot = (FFT_REAL *) mesh->rhok2;
    else if (!(mesh->rhot = FFT_MALLOC(num * sizeof(FFT_REAL)))) {
      *err = EZMOCK_ERR_MEMORY; return;
    }
  }

  const real thres = (rho_c == 0) ? REAL_MIN : rho_c;   /* skip 0 density */
  const double iexp = -1 / rho_exp;

  /* Estimate the maximum tracer number density for histogramming. */
  const real rho_max = 1 + lambda * 5;  /* 5-sigma limit of Gaussian randoms */
  const int imax = ilogb(rho_max);      /* expected maximum floor(log2(rhot)) */
  pdf->ilogb_max = imax;

#ifdef OMP
  size_t *ifirst = NULL;        /* first valid cell indices per thread */
  size_t *cnt = NULL;           /* number of valid cells for each thread */
  if (!(ifirst = calloc(conf->nthread, sizeof(size_t))) ||
      !(cnt = calloc(conf->nthread, sizeof(size_t)))) {
    if (ifirst) free(ifirst);
    *err = EZMOCK_ERR_MEMORY; return;
  }

  const size_t pnum = num / conf->nthread;
  const int rem = num % conf->nthread;

  /* Apply density thresholds and count cells to be processed per thread. */
#pragma omp parallel num_threads(conf->nthread)
  {
    /* First round: distribute grid points to threads uniformly. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the density field with OpenMP threads. */
    for (size_t i = istart; i < iend; i++) {
      if (mesh->rho[i] < thres) mesh->rhot[i] = 0;
      else {
        mesh->rhot[i] = (mesh->rho[i] > rho_sat) ? rho_sat : mesh->rho[i];
        ifirst[tid] = i;
        cnt[tid] = 1;
        break;
      }
    }
    for (size_t i = ifirst[tid] + 1; i < iend; i++) {
      if (mesh->rho[i] < thres) mesh->rhot[i] = 0;
      else {
        mesh->rhot[i] = (mesh->rho[i] > rho_sat) ? rho_sat : mesh->rho[i];
        cnt[tid]++;
      }
    }
  }

  /* Reassign even number of valid cells to threads, as Gaussian random numbers
     are generated in pairs. */
  for (int i = 0; i < conf->nthread - 1; i++) {
    /* Borrow one cell from the next thread, if the number of cells is odd. */
    if ((cnt[i] & 1) == 1) {
      ifirst[i + 1]++;
      cnt[i]++;
      cnt[i + 1]--;
    }
  }

  /* Compute the step size for jumping ahead the random state. */
  for (int i = 1; i < conf->nthread; i++)  cnt[i] += cnt[i - 1];
  if (cnt[conf->nthread - 1] < pdf->ncell) {
    free(ifirst); free(cnt);
    *err = EZMOCK_ERR_RHO_LESS; return;
  }

  /* Jump ahead the random states with a base step of Ng^3 * 2. */
#pragma omp parallel for num_threads(conf->nthread)
  for (int i = 0; i < conf->nthread; i++) {
    const size_t offset = (i == 0) ? 0 : cnt[i - 1];
    rng->reset(rng->state_stream[i], erng->seed, (num << 1) + offset, err);
  }
  free(cnt);
  if (PRAND_IS_ERROR(*err)) {
    free(ifirst);
    *err = EZMOCK_ERR_RNG_JUMP; return;
  }

  /* Apply the effective bias model. */
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Second round: each thread is assigned even number of valid cells. */
    const int tid = omp_get_thread_num();
    const size_t iend = (tid == conf->nthread - 1) ? num : ifirst[tid + 1];

    for (size_t i = ifirst[tid]; i < iend; i++) {
      if (mesh->rhot[i] == 0) continue;

      /* Gaussian random numbers for the stochaostic bias. */
      double g1, g2;
      gauss_rand(rng, rng->state_stream[tid], lambda, &g1, &g2);
      g1 = (g1 >= 0) ? 1 + g1 : exp(g1);
      g2 = (g2 >= 0) ? 1 + g2 : exp(g2);

      /* Map dark matter density to tracer number density. */
      double tmp = (1 - exp((double) mesh->rhot[i] * iexp)) * g1;
      mesh->rhot[i] = tmp;

      /* Update the histogram of log2(rhot). */
      int idx = imax - ilogb(mesh->rhot[i]);    /* large densities first */
      if (idx < 0) idx = 0;
      else if (idx >= EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN)
        idx = EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN - 1;
#pragma omp atomic
      pdf->hist[idx]++;

      /* Use the second random number for the next tracer density. */
      while (++i < iend) {
        if (mesh->rhot[i] == 0) continue;
        tmp = (1 - exp((double) mesh->rhot[i] * iexp)) * g2;
        mesh->rhot[i] = tmp;

        /* Update the histogram of log2(rhot). */
        idx = imax - ilogb(mesh->rhot[i]);      /* large densities first */
        if (idx < 0) idx = 0;
        else if (idx >= EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN)
          idx = EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN - 1;
#pragma omp atomic
        pdf->hist[idx]++;

        break;
      }
    }
  }
  free(ifirst);
#else
  /* Jump ahead the random state with a step of Ng^3 * 2. */
  rng->reset(rng->state, erng->seed, num << 1, err);
  if (PRAND_IS_ERROR(*err)) {
    *err = EZMOCK_ERR_RNG_JUMP; return;
  }

  /* Apply the effective bias model directly. */
  size_t cnt = 0;
  for (size_t i = 0; i < num; i++) {
    if (mesh->rho[i] < thres) {
      mesh->rhot[i] = 0;
      continue;
    }
    mesh->rhot[i] = (mesh->rho[i] > rho_sat) ? rho_sat : mesh->rho[i];
    cnt++;

    /* Gaussian random numbers for the stochaostic bias. */
    double g1, g2;
    gauss_rand(rng, rng->state, lambda, &g1, &g2);
    g1 = (g1 >= 0) ? 1 + g1 : exp(g1);
    g2 = (g2 >= 0) ? 1 + g2 : exp(g2);

    /* Map dark matter density to tracer number density. */
    mesh->rhot[i] = (1 - exp(mesh->rhot[i] * iexp)) * g1;

    /* Update the histogram of log2(rhot). */
    int idx = imax - ilogb(mesh->rhot[i]);    /* large densities first */
    if (idx < 0) idx = 0;
    else if (idx >= EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN)
      idx = EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN - 1;
    pdf->hist[idx]++;

    /* Use the second random number for the next tracer density. */
    while (++i < num) {
      if (mesh->rho[i] < thres) {
        mesh->rhot[i] = 0;
        continue;
      }
      mesh->rhot[i] = (mesh->rho[i] > rho_sat) ? rho_sat : mesh->rho[i];
      cnt++;
      mesh->rhot[i] = (1 - exp(mesh->rhot[i] * iexp)) * g2;

      /* Update the histogram of log2(rhot). */
      idx = imax - ilogb(mesh->rhot[i]);      /* large densities first */
      if (idx < 0) idx = 0;
      else if (idx >= EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN)
        idx = EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN - 1;
      pdf->hist[idx]++;

      break;
    }
  }
  if (cnt < pdf->ncell) {
    *err = EZMOCK_ERR_RHO_LESS; return;
  }
#endif
}

/******************************************************************************
Function `n_tracer_in_cell`:
  Evaluate the number of tracers to be assigned to each cell.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `pdf`:      instance of the EZmock tracer PDF;
  * `att_part`: indicate whether to attach tracers to particles;
  * `err`:      integer storing the error code.
******************************************************************************/
static void n_tracer_in_cell(EZMOCK *ez, EZMOCK_PDF *pdf, const bool att_part,
    int *err) {
  /* Get the density threshold given the histogram and desired cell number. */
  int imax;     /* maximum index of the density histogram to be considered */
  size_t cnt = 0;
  for (imax = 0; imax < EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN; imax++) {
    cnt += pdf->hist[imax];
    if (cnt >= pdf->ncell) break;
  }
  const real rho_thres = (imax == EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN - 1) ?
       REAL_MIN : pow(2, pdf->ilogb_max - imax);

  /* Allocate memory for indices of cells to be used. */
  size_t *cidx = malloc(cnt * sizeof(size_t));
  if (!cidx) {
    *err = EZMOCK_ERR_MEMORY; return;
  }

  /* Compute offsets of different histogram bins in the cell index array. */
  size_t *offset = calloc((imax + 1), sizeof(size_t));
  if (!offset) {
    free(cidx); *err = EZMOCK_ERR_MEMORY; return;
  }
  for (int i = 1; i <= imax; i++)
    offset[i] = offset[i - 1] + pdf->hist[i - 1];

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;

  /* Record indices of cells following the order of histogram bins. */
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
  for (size_t i = 0; i < num; i++) {
    if (mesh->rhot[i] < rho_thres) continue;
    /* Get the histogram bin of the tracer density. */
    int idx = pdf->ilogb_max - ilogb(mesh->rhot[i]);
    if (idx < 0) idx = 0;
    else if (idx > imax) idx = imax;
#ifdef OMP
#pragma omp critical
#endif
    cidx[offset[idx]++] = i;
  }
  free(offset);

  /* Duplicate the PDF. */
  size_t *cpdf = malloc(pdf->nmax * sizeof(size_t));
  if (!cpdf) {
    free(cidx); *err = EZMOCK_ERR_MEMORY; return;
  }
  memcpy(cpdf, pdf->pdf, pdf->nmax * sizeof(size_t));

  /* Adjust the order of cell indices according to the PDF bins. */
  cnt = 0;              /* number of cells in the correct order so far */
  size_t used = 0;
  int ipdf = 0;
  for (int i = 0; i <= imax; i++) {
    size_t nrest = pdf->hist[i] - used;
    /* Cells are already in the right order if they are in the same PDF bin. */
    if (nrest <= cpdf[ipdf]) {
      cpdf[ipdf] -= nrest;
      if (cpdf[ipdf] == 0) ipdf++;
      cnt += nrest;
      used = 0;
    }
    /* Select the largest `cpdf[ipdf]` densities from the histogram bin. */
    else {
      size_t buf;
      qselect(cidx + cnt, cpdf[ipdf], nrest, &buf, mesh->rhot);

      /* Check the minimum density of processed cells. */
      size_t ecnt = cnt + cpdf[ipdf];
      real min = mesh->rhot[cidx[ecnt - 1]];
      for (size_t k = cnt; k < ecnt - 1; k++)
        if (min > mesh->rhot[cidx[k]]) min = mesh->rhot[cidx[k]];
      /* Include extra cells with the same density. */
      for (size_t k = ecnt; k < cnt + nrest; k++) {
        if (mesh->rhot[cidx[k]] == min) {
          buf = cidx[ecnt];
          cidx[ecnt++] = cidx[k];
          cidx[k] = buf;
          pdf->pdf[ipdf]++;
          pdf->ncell++;
          pdf->ntot += pdf->nmax - ipdf;
        }
      }

      used += ecnt - cnt;
      cnt = ecnt;
      ipdf++;                           /* go to the next PDF bin */
      if (used < pdf->hist[i]) i--;     /* continue with the same hist bin */
      else used -= pdf->hist[i];
    }
    if (ipdf == pdf->nmax) break;
  }
  free(cpdf);

  /* Sanity check. */
  if (ipdf != pdf->nmax || cnt != pdf->ncell) {
    free(cidx); *err = EZMOCK_ERR_UNKNOWN; return;
  }

#ifdef EZMOCK_DEBUG
  /* Validate the order of cells: `cidx` should be sorted by
     the number of tracers per cell in descending order. */
  real *rmin, *rmax;
  size_t *cdf;
  if (!(rmin = malloc(pdf->nmax * sizeof(real))) ||
      !(rmax = malloc(pdf->nmax * sizeof(real))) ||
      !(cdf = calloc(pdf->nmax + 1, sizeof(size_t)))) {
    printf("Memory error!\n");
    exit(EZMOCK_ERR_MEMORY);
  }
  for (int i = 0; i < pdf->nmax; i++)
    cdf[i + 1] = cdf[i] + pdf->pdf[i];
  for (int i = 0; i < pdf->nmax; i++) {
    rmin[i] = mesh->rhot[cidx[cdf[i]]];
    rmax[i] = mesh->rhot[cidx[cdf[i]]];
    for (size_t j = cdf[i] + 1; j < cdf[i + 1]; j++) {
      if (rmin[i] > mesh->rhot[cidx[j]]) rmin[i] = mesh->rhot[cidx[j]];
      if (rmax[i] < mesh->rhot[cidx[j]]) rmax[i] = mesh->rhot[cidx[j]];
    }
  }
  for (int i = 1; i < pdf->nmax; i++) {
    if (rmin[i - 1] < rmax[i]) {
      printf("Wrong order of cell indices!\n");
      exit(EZMOCK_ERR_UNKNOWN);
    }
  }
  free(rmin); free(rmax); free(cdf);
#endif

  /* Recycle the tracer density field for numbers of particles and tracers.*/
  memset(mesh->rhot, 0, num * sizeof(FFT_REAL));
  /* Note that sizeof(FFT_REAL) >= sizeof(uint16_t) * 2. */
  mesh->nt = (int16_t *) mesh->rhot;
  mesh->np = (uint16_t *) mesh->rhot + num;

  /* Evaluate the numbers of tracers in each cell. */
  size_t iend = 0;
  for (int i = 0; i < pdf->nmax; i++) {
    size_t istart = iend;
    iend += pdf->pdf[i];
    int16_t nt = pdf->nmax - i;
#ifdef OMP
    /* Use multiple threads only if there are enough tasks per threads. */
    if (pdf->pdf[i] > (((size_t) conf->nthread) << 2)) {
#pragma omp parallel for num_threads(conf->nthread)
      for (size_t j = istart; j < iend; j++) mesh->nt[cidx[j]] = nt;
    }
    else
#endif
    for (size_t j = istart; j < iend; j++) mesh->nt[cidx[j]] = nt;
  }
  free(cidx);

  if (att_part) {
    /* Count the numbers of particles in each cell. */
    const real igs = conf->Ng / conf->Lbox;             /* inverse grid size */
#ifdef OMP
    const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
    const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#pragma omp parallel num_threads(conf->nthread)
    {
      /* Distribute mesh grids to threads. */
      const int tid = omp_get_thread_num();
      const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
      const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
      const size_t iend = istart + pcnt;

      /* 1st round: count the number of particles per grid cell. */
      for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
        int i = idx_ij / conf->Ng;
        int j = idx_ij % conf->Ng;
#else
    /* 1st round: count the number of particles per grid cell. */
    for (int i = 0; i < conf->Ng; i++) {
      size_t idx_i = (size_t) i * conf->Ng;
      for (int j = 0; j < conf->Ng; j++) {
        size_t idx_ij = idx_i + j;
#endif
        size_t idx0 = idx_ij * conf->Ng;
        for (int k = 0; k < conf->Ng; k++) {
          size_t idx = idx0 + k;
          real px = mesh->psi[0][idx] * igs + i;
          real py = mesh->psi[1][idx] * igs + j;
          real pz = mesh->psi[2][idx] * igs + k;

          /* Periodc boundary condition. */
          if (px < 0) px += conf->Ng;
          else if (px >= conf->Ng) px -= conf->Ng;
          if (py < 0) py += conf->Ng;
          else if (py >= conf->Ng) py -= conf->Ng;
          if (pz < 0) pz += conf->Ng;
          else if (pz >= conf->Ng) pz -= conf->Ng;

          /* Find the nearest grid point (NGP) of particle. */
          int ix = (int) (px + 0.5);
          int iy = (int) (py + 0.5);
          int iz = (int) (pz + 0.5);
          if (ix >= conf->Ng) ix = 0;
          if (iy >= conf->Ng) iy = 0;
          if (iz >= conf->Ng) iz = 0;

          /* Skip cells that do not host tracers. */
          size_t pidx = IDX(conf->Ng, ix, iy, iz);
          if (mesh->nt[pidx] == 0) continue;

          if (mesh->np[pidx] == UINT16_MAX) {
            *err = EZMOCK_ERR_RHO_NPART;
#ifdef OMP
            break;
#else
            return;
#endif
          }
OMP_ATOMIC
          mesh->np[pidx]++;
        }       /* for k */
      }         /* for idx_ij (OMP) or for j (no OMP) */
    }           /* omp parallel (OMP) or for i (no OMP) */
  }             /* if att_part */
}

/******************************************************************************
Function `disp_cic`:
  Interpolate the displacement field with cloud-in-cell (CIC).
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `dim`:      dimension (component) of the displacement to be interpolated;
  * `out`:      the interpolated displacement field.
******************************************************************************/
static void disp_cic(EZMOCK *ez, const int dim, FFT_REAL *out) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  const real igs = conf->Ng / conf->Lbox;       /* inverse grid size */

  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;
  memset(out, 0, num * sizeof(FFT_REAL));

#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the displacement field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int i = idx_ij / conf->Ng;
      int j = idx_ij % conf->Ng;
#else
  /* Traverse the displacement field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      size_t idx0 = idx_ij * conf->Ng;
      for (int k = 0; k < conf->Ng; k++) {
        size_t idx = idx0 + k;
        real x = mesh->psi[0][idx] * igs + i;
        real y = mesh->psi[1][idx] * igs + j;
        real z = mesh->psi[2][idx] * igs + k;

        int x0 = (int) floor(x);
        int y0 = (int) floor(y);
        int z0 = (int) floor(z);

        /* Weights for neighbours. */
        real wx1 = (x - x0);
        real wx0 = 1 - wx1;
        real wy1 = y - y0;
        real wy0 = 1 - wy1;
        real wz1 = z - z0;
        real wz0 = 1 - wz1;
        wx0 *= mesh->psi[dim][idx];
        wx1 *= mesh->psi[dim][idx];

        /* Periodic boundary condition. */
        if (x0 < 0) x0 += conf->Ng;
        else if (x0 >= conf->Ng) x0 -= conf->Ng;
        if (y0 < 0) y0 += conf->Ng;
        else if (y0 >= conf->Ng) y0 -= conf->Ng;
        if (z0 < 0) z0 += conf->Ng;
        else if (z0 >= conf->Ng) z0 -= conf->Ng;

        int x1 = (x0 == conf->Ng - 1) ? 0 : x0 + 1;
        int y1 = (y0 == conf->Ng - 1) ? 0 : y0 + 1;
        int z1 = (z0 == conf->Ng - 1) ? 0 : z0 + 1;

OMP_ATOMIC
        out[IDX(conf->Ng,x0,y0,z0)] += wx0 * wy0 * wz0;
OMP_ATOMIC
        out[IDX(conf->Ng,x0,y0,z1)] += wx0 * wy0 * wz1;
OMP_ATOMIC
        out[IDX(conf->Ng,x0,y1,z0)] += wx0 * wy1 * wz0;
OMP_ATOMIC
        out[IDX(conf->Ng,x0,y1,z1)] += wx0 * wy1 * wz1;
OMP_ATOMIC
        out[IDX(conf->Ng,x1,y0,z0)] += wx1 * wy0 * wz0;
OMP_ATOMIC
        out[IDX(conf->Ng,x1,y0,z1)] += wx1 * wy0 * wz1;
OMP_ATOMIC
        out[IDX(conf->Ng,x1,y1,z0)] += wx1 * wy1 * wz0;
OMP_ATOMIC
        out[IDX(conf->Ng,x1,y1,z1)] += wx1 * wy1 * wz1;
      }         /* for k */
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for idx (no OMP) */

  /* Normalize the interpolated displacement field. */
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
  for (size_t i = 0; i < num; i++) {
    if (mesh->rho[i] > REAL_TOL) out[i] /= mesh->rho[i];
  }
}

/******************************************************************************
Function `assign_vel`:
  Assign peculiar velocities to tracers with a displacement field on grids.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `psi`:      the interpolated displacement field along a given dimension;
  * `n`:        number of tracers to be assigned velocities;
  * `x`:        x coordinates of the tracers;
  * `y`:        y coordinates of the tracers;
  * `z`:        z coordinates of the tracers;
  * `v`:        component of tracer peculiar velocities along a given dimension.
******************************************************************************/
static void assign_vel(EZMOCK *ez, const FFT_REAL *psi, const size_t n,
    const real *x, const real *y, const real *z, real *v) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  const real igs = conf->Ng / conf->Lbox;       /* inverse grid size */

#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
  for (size_t i = 0; i < n; i++) {
    real xx = x[i] * igs;
    real yy = y[i] * igs;
    real zz = z[i] * igs;

    int x0 = (int) xx;
    int y0 = (int) yy;
    int z0 = (int) zz;

    /* Weights for neighbours. */
    real wx1 = xx - x0;
    real wx0 = 1 - wx1;
    real wy1 = yy - y0;
    real wy0 = 1 - wy1;
    real wz1 = zz - z0;
    real wz0 = 1 - wz1;

    if (x0 == conf->Ng) x0 = 0;
    if (y0 == conf->Ng) y0 = 0;
    if (z0 == conf->Ng) z0 = 0;
    int x1 = (x0 == conf->Ng - 1) ? 0 : x0 + 1;
    int y1 = (y0 == conf->Ng - 1) ? 0 : y0 + 1;
    int z1 = (z0 == conf->Ng - 1) ? 0 : z0 + 1;

    /* Weight the velocity contribution by the dark matter density. */
    const real r000 = mesh->rho[IDX(conf->Ng,x0,y0,z0)] * wx0 * wy0 * wz0;
    const real r001 = mesh->rho[IDX(conf->Ng,x0,y0,z1)] * wx0 * wy0 * wz1;
    const real r010 = mesh->rho[IDX(conf->Ng,x0,y1,z0)] * wx0 * wy1 * wz0;
    const real r011 = mesh->rho[IDX(conf->Ng,x0,y1,z1)] * wx0 * wy1 * wz1;
    const real r100 = mesh->rho[IDX(conf->Ng,x1,y0,z0)] * wx1 * wy0 * wz0;
    const real r101 = mesh->rho[IDX(conf->Ng,x1,y0,z1)] * wx1 * wy0 * wz1;
    const real r110 = mesh->rho[IDX(conf->Ng,x1,y1,z0)] * wx1 * wy1 * wz0;
    const real r111 = mesh->rho[IDX(conf->Ng,x1,y1,z1)] * wx1 * wy1 * wz1;
    const real rho = r000 + r001 + r010 + r011 + r100 + r101 + r110 + r111;

    v[i] = psi[IDX(conf->Ng,x0,y0,z0)] * r000
        + psi[IDX(conf->Ng,x0,y0,z1)] * r001
        + psi[IDX(conf->Ng,x0,y1,z0)] * r010
        + psi[IDX(conf->Ng,x0,y1,z1)] * r011
        + psi[IDX(conf->Ng,x1,y0,z0)] * r100
        + psi[IDX(conf->Ng,x1,y0,z1)] * r101
        + psi[IDX(conf->Ng,x1,y1,z0)] * r110
        + psi[IDX(conf->Ng,x1,y1,z1)] * r111;

    v[i] *= cosmo->vfac / rho;
  }
}

/******************************************************************************
Function `generate_tracers`:
  Create the tracer catalogue given the cells hosting tracers and tracer PDF
  without parallelisation.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `pdf`:      instance of the EZmock tracer PDF;
  * `att_part`: indicate whether to attach tracers to particles;
  * `x`:        array to be filled x coordinates of the tracers;
  * `y`:        array to be filled y coordinates of the tracers;
  * `z`:        array to be filled z coordinates of the tracers;
  * `vx`:       array to be filled velocities of the tracers along x;
  * `vy`:       array to be filled velocities of the tracers along y;
  * `vz`:       array to be filled velocities of the tracers along z;
  * `err`:      integer storing the error code.
Return:
  The tracer catalogue as (x,y,z,vx,vy,vz) arrays on success; NULL on error.
******************************************************************************/
static void generate_tracers(EZMOCK *ez, const EZMOCK_PDF *pdf,
    const bool att_part, real **x, real **y, real **z,
    real **vx, real **vy, real **vz, int *err) {
  /* Jump ahead the random state with a base step of Ng^3 * 4. */
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_RNG *erng = (EZMOCK_RNG *) ez->rng;
  size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;
  prand_t *rng = erng->rng;
  rng->reset(rng->state, erng->seed, num << 2, err);
  if (PRAND_IS_ERROR(*err)) {
    *err = EZMOCK_ERR_RNG_JUMP; return;
  }

  /* Allocate memory for the tracer positions and velocities. */
  real *cat[6];
  for (int i = 0; i < 6; i++) {
    if (!(cat[i] = malloc(pdf->ntot * sizeof(real)))) {
      for (int j = 0; j < i; j++) free(cat[j]);
      *err = EZMOCK_ERR_MEMORY; return;
    }
  }
  real *pos[3], *vel[3];
  for (int i = 0; i < 3; i++) {
    pos[i] = cat[i];
    vel[i] = cat[i + 3];
  }

  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  const real gs = conf->Lbox / conf->Ng;                /* grid size */
  size_t iv = 0;        /* starting index of tracers without particles */

  if (att_part) {
#ifndef EZMOCK_VELOCITY_ALL_INTERP
    EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
#endif
    const real igs = conf->Ng / conf->Lbox;             /* inverse grid size */
    size_t cnt = 0;

    /* Assign tracers to particles: this part has to be run serially. */
    for (int i = 0; i < conf->Ng; i++) {
      size_t idx_i = (size_t) i * conf->Ng;
      for (int j = 0; j < conf->Ng; j++) {
        size_t idx0 = (idx_i + j) * conf->Ng;
        for (int k = 0; k < conf->Ng; k++) {
          size_t idx = idx0 + k;
          real px = mesh->psi[0][idx] * igs + i;
          real py = mesh->psi[1][idx] * igs + j;
          real pz = mesh->psi[2][idx] * igs + k;

          /* Periodc boundary condition. */
          if (px < 0) px += conf->Ng;
          else if (px >= conf->Ng) px -= conf->Ng;
          if (py < 0) py += conf->Ng;
          else if (py >= conf->Ng) py -= conf->Ng;
          if (pz < 0) pz += conf->Ng;
          else if (pz >= conf->Ng) pz -= conf->Ng;

          /* Find the nearest grid point (NGP) of particle. */
          int ix = (int) (px + 0.5);
          int iy = (int) (py + 0.5);
          int iz = (int) (pz + 0.5);
          if (ix >= conf->Ng) ix = 0;
          if (iy >= conf->Ng) iy = 0;
          if (iz >= conf->Ng) iz = 0;

          /* Skip cells that do not host tracers. */
          size_t pidx = IDX(conf->Ng, ix, iy, iz);
          if (mesh->nt[pidx] == 0) continue;

          /* Randomly select particles if not all of them assigned tracers. */
          if (mesh->nt[pidx] < mesh->np[pidx]) {
            /* Randomly select `m` elements from array with length `n`.
               See section 12.2 of Programming Pearls 2nd Edition (Bentley). */
            double ran = rng->get_double(rng->state);
            if (ran * mesh->np[pidx] >= mesh->nt[pidx]) {
              mesh->np[pidx]--;
              continue;
            }
          }

          /* Now assign the tracer to the particle. */
          mesh->nt[pidx]--;
          mesh->np[pidx]--;
          pos[0][cnt] = px * gs;
          pos[1][cnt] = py * gs;
          pos[2][cnt] = pz * gs;
  #ifndef EZMOCK_VELOCITY_ALL_INTERP
          vel[0][cnt] = mesh->psi[0][idx] * cosmo->vfac;
          vel[1][cnt] = mesh->psi[1][idx] * cosmo->vfac;
          vel[2][cnt] = mesh->psi[2][idx] * cosmo->vfac;
  #endif
          cnt++;
        }
      }
    }

    iv = cnt;                   /* index of tracers without particles */
    if (iv == pdf->ntot) {
      *x = pos[0];
      *y = pos[1];
      *z = pos[2];
      *vx = vel[0];
      *vy = vel[1];
      *vz = vel[2];
      return;
    }
  }     /* if att_part */

#ifdef OMP
  /* Allocate memory for jumping ahead random states. */
  size_t *rcnt = calloc(conf->nthread, sizeof(size_t));
  if (!rcnt) {
    *err = EZMOCK_ERR_MEMORY;
    for (int i = 0; i < 6; i++) free(cat[i]);
    return;
  }
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;

  /* Generate tracers without host particles. */
#pragma omp parallel num_threads(conf->nthread)
  {
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Count the number of tracers to be generated at each thread. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      size_t idx0 = idx_ij * conf->Ng;
      for (int k = 0; k < conf->Ng; k++) {
        size_t idx = idx0 + k;
        if (mesh->nt[idx] > 0) rcnt[tid] += mesh->nt[idx];
      }
    }
#pragma omp barrier             /* synchronise `rcnt` */
#pragma omp single
    {
      for (int i = 1; i < conf->nthread; i++) rcnt[i] += rcnt[i - 1];
      for (int i = conf->nthread - 1; i > 0; i--) rcnt[i] = rcnt[i - 1];
      rcnt[0] = 0;
    }
#pragma omp barrier             /* synchronise `rcnt` */
    /* Jump ahead random states with a base step of Ng^3 * 5. */
    size_t offset = rcnt[tid];
    /* 3 random numbers for uniform random positions. */
    rng->reset(rng->state_stream[tid], erng->seed, num * 5 + offset * 3, err);
    if (PRAND_IS_ERROR(*err)) *err = EZMOCK_ERR_RNG_JUMP;
#pragma omp barrier             /* synchronise `err` */

    if (*err == EZMOCK_SUCCESS) {
      offset += iv;
      for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
        int i = idx_ij / conf->Ng;
        int j = idx_ij % conf->Ng;
        size_t idx0 = idx_ij * conf->Ng;
#else
  /* Jump ahead the random state with a step of Ng^3 * 5. */
  rng->reset(rng->state, erng->seed, num * 5, err);
  if (PRAND_IS_ERROR(*err)) {
    *err = EZMOCK_ERR_RNG_JUMP;
    for (int i = 0; i < 6; i++) free(cat[i]);
    return;
  }
  size_t cnt = 0;

  /* Generate tracers without host particles. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx0 = (idx_i + j) * conf->Ng;
#endif

        for (int k = 0; k < conf->Ng; k++) {
          size_t idx = idx0 + k;
          if (mesh->nt[idx] == 0) continue;
          for (int n = 0; n < mesh->nt[idx]; n++) {
#ifdef OMP
            double r = rng->get_double(rng->state_stream[tid]) - 0.5;
#else
            double r = rng->get_double(rng->state) - 0.5;
#endif
            real xx = (i + r) * gs;
#ifdef OMP
            r = rng->get_double(rng->state_stream[tid]) - 0.5;
#else
            r = rng->get_double(rng->state) - 0.5;
#endif
            real yy = (j + r) * gs;
#ifdef OMP
            r = rng->get_double(rng->state_stream[tid]) - 0.5;
#else
            r = rng->get_double(rng->state) - 0.5;
#endif
            real zz = (k + r) * gs;

            if (xx < 0) xx += conf->Lbox;
            if (yy < 0) yy += conf->Lbox;
            if (zz < 0) zz += conf->Lbox;

#ifdef OMP
            pos[0][offset] = xx;
            pos[1][offset] = yy;
            pos[2][offset] = zz;
            offset++;
#else
            pos[0][cnt] = xx;
            pos[1][cnt] = yy;
            pos[2][cnt] = zz;
            cnt++;
#endif
          }     /* for n */
        }       /* for k */
      }         /* for idx_ij (OMP) or for j (no OMP) */
    }           /* *err == EZMOCK_SUCCESS (OMP) or for i (no OMP) */
#ifdef OMP
    /* Make sure that all tracers are generated without duplicates. */
    if (tid == conf->nthread - 1) {
      if (offset != pdf->ntot) *err = EZMOCK_ERR_UNKNOWN;
    }
    else if (offset != rcnt[tid + 1] + iv) *err = EZMOCK_ERR_UNKNOWN;
  }             /* omp parallel */

  free(rcnt);
  if (*err != EZMOCK_SUCCESS) {
    for (int i = 0; i < 6; i++) free(cat[i]);
    return;
  }
#else
  /* Now the positions of all tracers should have been generated. */
  if (cnt != pdf->ntot) {
    for (int i = 0; i < 6; i++) free(cat[i]);
    *err = EZMOCK_ERR_UNKNOWN; return;
  }
#endif

#ifdef EZMOCK_VELOCITY_ALL_INTERP
  iv = 0;
#endif
  /* Assign peculiar velocities to randomly sampled tracers with CIC. */
  for (int dim = 0; dim < 3; dim++) {
    /* Paint the displacement field on mesh with CIC. */
    disp_cic(ez, dim, mesh->rhot);
    /* Assign peculiar velocities. */
    assign_vel(ez, mesh->rhot, pdf->ntot - iv, pos[0] + iv, pos[1] + iv,
        pos[2] + iv, vel[dim] + iv);
  }

  *x = pos[0];
  *y = pos[1];
  *z = pos[2];
  *vx = vel[0];
  *vy = vel[1];
  *vz = vel[2];
}

/******************************************************************************
Function `add_vel_scatter`:
  Add Gaussian scatter to the peculiar velocities of tracers.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `vx`:       velocities of the tracers along x direction;
  * `vy`:       velocities of the tracers along y direction;
  * `vz`:       velocities of the tracers along z direction;
  * `n`:        number of tracers in total;
  * `sigma_v`:  width of the Gaussian distribution of randoms;
  * `err`:      integer storing the error code.
Return:
  The tracer catalogue as (x,y,z,vx,vy,vz) arrays on success; NULL on error.
******************************************************************************/
static void add_vel_scatter(EZMOCK *ez, real *vx, real *vy, real *vz,
    const size_t n, const real sigma_v, int *err) {
  EZMOCK_RNG *erng = (EZMOCK_RNG *) ez->rng;
  prand_t *rng = erng->rng;
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;

#ifdef OMP
  const size_t pnum = n / conf->nthread;
  const int rem = n % conf->nthread;
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute even numbers of tracers to threads. */
    const int tid = omp_get_thread_num();
    size_t pcnt, istart, iend;
    if (tid < rem) {
      pcnt = pnum + (tid & 1) + ((pnum ^ tid) & 1);
      istart = (pnum + 1) * tid - (tid & (~pnum) & 1);
    }
    else {
      pcnt = pnum + (pnum & 1) * ((((tid ^ rem) & 1) << 1) - 1);
      istart = pnum * tid + rem - (rem & (~pnum) & 1)
          - ((tid ^ rem) & pnum & 1);
    }
    if (tid == conf->nthread - 1) pcnt = n - istart;
    iend = istart + pcnt;

    /* Jump ahead random states with a base step of Ng^3 * 8. */
    rng->reset(rng->state_stream[tid], erng->seed, (num << 3) + istart * 3,
        err);
    if (PRAND_IS_ERROR(*err)) *err = EZMOCK_ERR_RNG_JUMP;
#pragma omp barrier             /* synchronize `err` */

    if (*err == EZMOCK_SUCCESS) {
      for (size_t i = istart; i < iend; i++) {
        double g1, g2;
        gauss_rand(rng, rng->state_stream[tid], sigma_v, &g1, &g2);
        vx[i] += g1;
        vy[i] += g2;
        gauss_rand(rng, rng->state_stream[tid], sigma_v, &g1, &g2);
        vz[i++] += g1;
        if(i == iend) break;
        vx[i] += g2;
        gauss_rand(rng, rng->state_stream[tid], sigma_v, &g1, &g2);
        vy[i] += g1;
        vz[i] += g2;
      }
    }
  }
#else
  /* Jump ahead the random state with a base step of Ng^3 * 8. */
  rng->reset(rng->state, erng->seed, num << 3, err);
  if (PRAND_IS_ERROR(*err)) {
    *err = EZMOCK_ERR_RNG_JUMP; return;
  }

  for (size_t i = 0; i < n; i++) {
    double g1, g2;
    gauss_rand(rng, rng->state, sigma_v, &g1, &g2);
    vx[i] += g1;
    vy[i] += g2;
    gauss_rand(rng, rng->state, sigma_v, &g1, &g2);
    vz[i++] += g1;
    if(i == n) break;
    vx[i] += g2;
    gauss_rand(rng, rng->state, sigma_v, &g1, &g2);
    vy[i] += g1;
    vz[i] += g2;
  }
#endif
}


/*============================================================================*\
                Interface for constructing the tracer catalogue
\*============================================================================*/

/******************************************************************************
Function `EZmock_populate_tracer`:
  Create the tracer catalogue based on the effective bias model.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `params`:   free parameters of the effective bias model:
                density_cut, density_exp, pdf_base, random_motion;
  * `nexp`:     expected number of tracers to be generated;
  * `att_part`: indicate whether to attach tracers to particles;
  * `ntracer`:  actual number of tracers generated;
  * `x`:        array to be filled x coordinates of the tracers;
  * `y`:        array to be filled y coordinates of the tracers;
  * `z`:        array to be filled z coordinates of the tracers;
  * `vx`:       array to be filled velocities of the tracers along x;
  * `vy`:       array to be filled velocities of the tracers along y;
  * `vz`:       array to be filled velocities of the tracers along z;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int EZmock_populate_tracer(EZMOCK *ez, const real *params, const size_t nexp,
    const bool att_part, size_t *ntracer, real **x, real **y, real **z,
    real **vx, real **vy, real **vz, int *err) {
  /* Validate arguments. */
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return (*err = EZMOCK_ERR_ARG_EZ);
  if (!params) return (*err = EZMOCK_ERR_ARG_PARAMS);
  if (nexp < EZMOCK_MIN_NUM_TRACER) return (*err = EZMOCK_ERR_ARG_LOWEXP);
  if (!ntracer) return (*err = EZMOCK_ERR_ARG_NDATA);
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  if (nexp >= (size_t) conf->Ng * conf->Ng * conf->Ng)
    return (*err = EZMOCK_ERR_ARG_HIEXP);

  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  if (!mesh->psi[0] || !mesh->psi[1] || !mesh->psi[2] || !mesh->rho)
    return (*err = EZMOCK_ERR_ARG_MESH);
  if (mesh->rho_replaced) return (*err = EZMOCK_ERR_RHO_USED);

  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  if (cosmo->vfac == HUGE_VAL) return (*err = EZMOCK_ERR_ARG_COSMO);

  /* Load effective bias model parameters. */
  const real rho_c = params[0];
  const real rho_exp = params[1];
  const real pbase = params[2];
  const real sigma_v = params[3];
  if (!isfinite(rho_c) || rho_c < 0) return (*err = EZMOCK_ERR_PAR_RCUT);
  if (!isfinite(rho_exp) || rho_exp <= 0) return (*err = EZMOCK_ERR_PAR_REXP);
  if (!isfinite(pbase) || pbase <= 0 || pbase >= 1)
    return (*err = EZMOCK_ERR_PAR_PBASE);
  if (!isfinite(sigma_v) || sigma_v < 0) return (*err = EZMOCK_ERR_PAR_SIGV);

  /* Evaluate the PDF of tracers given the base number. */
  EZMOCK_PDF *pdf = EZmock_pdf_init();
  if (!pdf) return (*err = EZMOCK_ERR_MEMORY);
  eval_pdf(ez, nexp, pbase, pdf, err);
  if (*err != EZMOCK_SUCCESS) {
    EZmock_pdf_destroy(pdf); return *err;
  }

  /* Modify the density field based on the bias model. */
  bias_model(ez, pdf, rho_c, EZMOCK_MODEL_RHO_SATURATION, rho_exp,
      EZMOCK_MODEL_GAUSS_WIDTH, err);
  if (*err != EZMOCK_SUCCESS) {
    EZmock_pdf_destroy(pdf); return *err;
  }

  /* Find indices of cells hosting tracers, sorted by the number of tracers
     per cell in descending order. */
  n_tracer_in_cell(ez, pdf, att_part, err);
  if (*err != EZMOCK_SUCCESS) {
    EZmock_pdf_destroy(pdf); return *err;
  }

  /* Populate tracers. */
  generate_tracers(ez, pdf, att_part, x, y, z, vx, vy, vz, err);
  if (*err != EZMOCK_SUCCESS) {
    EZmock_pdf_destroy(pdf); return *err;
  }

  *ntracer = pdf->ntot;
  EZmock_pdf_destroy(pdf);

  /* Add Gaussian scatter to the peculiar velocities. */
  add_vel_scatter(ez, *vx, *vy, *vz, *ntracer, sigma_v, err);
  if (*err != EZMOCK_SUCCESS) {
    free(*x); free(*y); free(*z);
    free(*vx); free(*vy); free(*vz);
    return *err;
  }

  return EZMOCK_SUCCESS;
}
