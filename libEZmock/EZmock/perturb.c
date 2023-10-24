/*******************************************************************************
* perturb.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com> [GPLv3 license]
 
*******************************************************************************/

/* Macros for the template functions. */
#if !defined(EZMOCK_PT_WHITENOISE) && defined(EZMOCK_PT_LOGPK) && \
  defined(EZMOCK_PT_FIXAMP) && defined(EZMOCK_PT_IPHASE)

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c,d)         a##b##c##d
#endif

#ifndef EZMOCK_PT_FUNCNAME
  #define EZMOCK_PT_FUNCNAME(a,b,c,d)   CONCAT_FNAME(a,b,c,d)
#endif

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef EZMOCK_LOGPK_NAME
  #undef EZMOCK_LOGPK_NAME
#endif
#ifdef EZMOCK_FIXAMP_NAME
  #undef EZMOCK_FIXAMP_NAME
#endif
#ifdef EZMOCK_IPHASE_NAME
  #undef EZMOCK_IPHASE_NAME
#endif

#if     EZMOCK_PT_LOGPK == 1
  #define EZMOCK_LOGPK_NAME     _logpk
#elif   EZMOCK_PT_LOGPK == 0
  #define EZMOCK_LOGPK_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_LOGPK`"
#endif

#if     EZMOCK_PT_FIXAMP == 1
  #define EZMOCK_FIXAMP_NAME    _fixamp
#elif   EZMOCK_PT_FIXAMP == 0
  #define EZMOCK_FIXAMP_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_FIXAMP`"
#endif

#if     EZMOCK_PT_IPHASE == 1
  #define EZMOCK_IPHASE_NAME    _iphase
#elif   EZMOCK_PT_IPHASE == 0
  #define EZMOCK_IPHASE_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_IPHASE`"
#endif


/*============================================================================*\
                   Function for displacement field generation
\*============================================================================*/

/******************************************************************************
Function `EZmock_ZA_disp<EZMOCK_LOGPK_NAME><EZMOCK_FIXAMP_NAME>
    <EZMOCK_PT_IPHASE>`:
  Generate the Zel'dovich displacement field given the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/
static void EZMOCK_PT_FUNCNAME(EZmock_ZA_disp, EZMOCK_LOGPK_NAME,
    EZMOCK_FIXAMP_NAME, EZMOCK_IPHASE_NAME)
    (EZMOCK *ez, FFT_PLAN plan, int *err) {
  /* Initialise the random number generator. */
  EZMOCK_RNG *erng = (EZMOCK_RNG *) ez->rng;
  prand_t *rng = erng->rng;

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  const int Ngh = conf->Ng >> 1;
  const int Ngk = Ngh + 1;
#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#endif

  const double kfac = M_PI * 2 / conf->Lbox;
#if EZMOCK_PT_LOGPK == 1
  const double logkfac = log(kfac);
#endif
  const double Pfac = pow(conf->Lbox, -1.5);

  /* First dimension: generate initial condition and compute psi[0]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;
    /* Jump ahead the random states. Note that i=j=k=0 is skipped. */
  #if EZMOCK_PT_FIXAMP == 1
    rng->reset(rng->state_stream[tid], erng->seed,
        (tid == 0) ? 0 : istart * Ngk - 1, err);
  #else
    rng->reset(rng->state_stream[tid], erng->seed,
        (tid == 0) ? 0 : (istart * Ngk - 1) << 1, err);
  #endif
    if (PRAND_IS_ERROR(*err)) *err = EZMOCK_ERR_RNG_JUMP;
#pragma omp barrier     /* synchronize `err` */

    /* Sample the Fourier space density field with OpenMP. */
    if (*err == EZMOCK_SUCCESS) {
      for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
        int i = idx_ij / conf->Ng;
        int j = idx_ij % conf->Ng;
        double ki = (i <= Ngh) ? i : i - conf->Ng;
        int ni = (i == 0) ? 0 : conf->Ng - i;           /* index of -i */
        double kni = (ni <= Ngh) ? ni : ni - conf->Ng;
#else
    /* Sample the Fourier space density field sequentially. */
    for (int i = 0; i < conf->Ng; i++) {
      size_t idx_i = (size_t) i * conf->Ng;
      double ki = (i <= Ngh) ? i : i - conf->Ng;
      int ni = (i == 0) ? 0 : conf->Ng - i;             /* index of -i */
      double kni = (ni <= Ngh) ? ni : ni - conf->Ng;
      for (int j = 0; j < conf->Ng; j++) {
        size_t idx_ij = idx_i + j;
#endif
        double kj = (j <= Ngh) ? j : j - conf->Ng;
        int nj = (j == 0) ? 0 : conf->Ng - j;           /* index of -j */
        double kij = ki * ki + kj * kj;
        size_t idx0 = idx_ij * Ngk;

        /* Expand k loop and treat 0 and Nyquist frequencies separately. */
        /* k = 0 */
        if (idx0 == 0) {                /* i = j = k = 0 */
          mesh->rhok[0][0] = mesh->rhok[0][1] = 0;
          mesh->rhok2[0][0] = mesh->rhok2[0][1] = 0;
        }
        else {
          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
  #else
          double amp = rng->get_double_pos(rng->state);
  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
#endif

          /* Ensure conjugation on the k = 0 plane. */
          /* Fill the field with the lower half plane. */
          if ((i == ni && j <= Ngh) || (i != ni && i <= Ngh)) {
            size_t idx = idx0;
            size_t nidx = ((size_t) ni * conf->Ng + nj) * Ngk;  /* (-i,-j,0) */
            double ksq = kij;
#if EZMOCK_PT_LOGPK == 1
            double kmod = log(ksq) * 0.5 + logkfac;     /* log(sqrt(ksq)) */
#else
            double kmod = kfac * sqrt(ksq);             /* sqrt(ksq) */
#endif
            /* Obtain the power with interpolation. */
            double P = pk_interp(pk, kmod);
#if EZMOCK_PT_FIXAMP == 0
            amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
            P = exp(P * 0.5) * sqrt(amp);
  #else
            P = exp(P * 0.5);                           /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
            P = sqrt(P * amp);
  #else
            P = sqrt(P);
  #endif
#endif
            P *= Pfac;

            /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
            phase += M_PI;
#endif
            P /= (ksq * kfac);
            mesh->rhok[idx][0] = -P * sin(phase);
            mesh->rhok[idx][1] = P * cos(phase);
            mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
            mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;

            mesh->rhok[nidx][0] = -mesh->rhok[idx][0];
            mesh->rhok[nidx][1] = mesh->rhok[idx][1];
            mesh->rhok2[nidx][0] = mesh->rhok[nidx][0] * kni;
            mesh->rhok2[nidx][1] = mesh->rhok[nidx][1] * kni;
          }             /* check of the lower half plane */
        }               /* i = j = k = 0 */

        /* 0 < k < k_ny */
        for (int k = 1; k <= ((conf->Ng - 1) >> 1); k++) {
          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
  #else
          double amp = rng->get_double_pos(rng->state);
  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
#endif

          size_t idx = idx0 + k;
          double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
          double kmod = log(ksq) * 0.5 + logkfac;       /* log(sqrt(ksq)) */
#else
          double kmod = kfac * sqrt(ksq);               /* sqrt(ksq) */
#endif
          /* Obtain the power with interpolation. */
          double P = pk_interp(pk, kmod);

#if EZMOCK_PT_FIXAMP == 0
          amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
          P = exp(P * 0.5) * sqrt(amp);
  #else
          P = exp(P * 0.5);                             /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
          P = sqrt(P * amp);
  #else
          P = sqrt(P);
  #endif
#endif
          P *= Pfac;

          /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
          phase += M_PI;
#endif
          P /= (ksq * kfac);
          mesh->rhok[idx][0] = -P * sin(phase);
          mesh->rhok[idx][1] = P * cos(phase);
          mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
          mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;
        }               /* for k */

        /* k = k_ny */
        if ((conf->Ng & 1) == 0) {
          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
  #else
          double amp = rng->get_double_pos(rng->state);
  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
#endif

          /* Ensure conjugation on the k = k_ny plane. */
          /* Fill the field with the lower half plane. */
          if ((i == ni && j <= Ngh) || (i != ni && i <= Ngh)) {
            int k = Ngh;
            size_t idx = idx0 + k;
            size_t nidx = ((size_t) ni * conf->Ng + nj) * Ngk + k;
            double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
            double kmod = log(ksq) * 0.5 + logkfac;     /* log(sqrt(ksq)) */
#else
            double kmod = kfac * sqrt(ksq);             /* sqrt(ksq) */
#endif
            /* Obtain the power with interpolation. */
            double P = pk_interp(pk, kmod);
#if EZMOCK_PT_FIXAMP == 0
            amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
            P = exp(P * 0.5) * sqrt(amp);
  #else
            P = exp(P * 0.5);                           /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
            P = sqrt(P * amp);
  #else
            P = sqrt(P);
  #endif
#endif
            P *= Pfac;

            /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
            phase += M_PI;
#endif
            P /= (ksq * kfac);
            mesh->rhok[idx][0] = -P * sin(phase);
            mesh->rhok[idx][1] = P * cos(phase);
            mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
            mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;

            mesh->rhok[nidx][0] = -mesh->rhok[idx][0];
            mesh->rhok[nidx][1] = mesh->rhok[idx][1];
            mesh->rhok2[nidx][0] = mesh->rhok[nidx][0] * kni;
            mesh->rhok2[nidx][1] = mesh->rhok[nidx][1] * kni;
          }     /* check of the lower half plane */
        }       /* Ng is even */
      }         /* for idx_ij (OMP) or for j (no OMP) */
    }           /* *err == EZMOCK_SUCCESS (OMP) or for i (no OMP) */

#ifdef OMP
  }             /* omp parallel */
  if (*err != EZMOCK_SUCCESS) return;
#endif

  /* The imaginary parts of white noise field should be 0 at:
     (0,0,0), (0,0,Ngh), (0,Ngh,0), (0,Ngh,Ngh),
     (Ngh,0,0), (Ngh,0,Ngh), (Ngh,Ngh,0), (Ngh,Ngh,Ngh) */
  if ((conf->Ng & 1) == 0) {
    size_t idx = Ngh;
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = (size_t) Ngh * Ngk;
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = (size_t) Ngh * (Ngk + 1);
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = (size_t) Ngh * conf->Ng * Ngk;
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = Ngh * ((size_t) conf->Ng * Ngk + 1);
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = (size_t) Ngh * (conf->Ng + 1) * Ngk;
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
    idx = ((size_t) Ngh * conf->Ng + Ngh) * Ngk + Ngh;
    mesh->rhok[idx][0] = mesh->rhok2[idx][0] = 0;
  }

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[0]);

  /* Second and third dimension: compute psi[1] and psi[2]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int j = idx_ij % conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      double kj = (j <= Ngh) ? j : j - conf->Ng;
      size_t idx0 = idx_ij * Ngk;

      /* 0 <= k <= k_ny */
      for (int k = 0; k < Ngk; k++) {
        size_t idx = idx0 + k;
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * kj;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * kj;
        mesh->rhok[idx][0] *= k;
        mesh->rhok[idx][1] *= k;
      }
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[1]);
  FFT_EXEC_C2R(plan, mesh->rhok, mesh->psi[2]);
}


#undef EZMOCK_LOGPK_NAME
#undef EZMOCK_FIXAMP_NAME
#undef EZMOCK_IPHASE_NAME

#undef EZMOCK_PT_LOGPK
#undef EZMOCK_PT_FIXAMP
#undef EZMOCK_PT_IPHASE


/******************************************************************************/

/* Macros for the template functions. */
#elif defined(EZMOCK_PT_WHITENOISE) && defined(EZMOCK_PT_LOGPK) && \
  !defined(EZMOCK_PT_FIXAMP) && !defined(EZMOCK_PT_IPHASE)

/* Macros for generating function names. */
#ifndef CONCAT_FNAME2
  #define CONCAT_FNAME2(a,b)            a##b
#endif

#ifndef EZMOCK_PT_FUNCNAME2
  #define EZMOCK_PT_FUNCNAME2(a,b)      CONCAT_FNAME2(a,b)
#endif

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef EZMOCK_LOGPK_NAME
  #undef EZMOCK_LOGPK_NAME
#endif

#if     EZMOCK_PT_LOGPK == 1
  #define EZMOCK_LOGPK_NAME     _logpk
#elif   EZMOCK_PT_LOGPK == 0
  #define EZMOCK_LOGPK_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_LOGPK`"
#endif


/*============================================================================*\
                   Function for displacement field generation
\*============================================================================*/

/******************************************************************************
Function `EZmock_ZA_disp_wn<EZMOCK_LOGPK_NAME>`:
  Generate the Zel'dovich displacement field given the white noise field and
  the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/
static void EZMOCK_PT_FUNCNAME2(EZmock_ZA_disp_wn, EZMOCK_LOGPK_NAME)
    (EZMOCK *ez, FFT_PLAN plan, int *err) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  const int Ngk = (conf->Ng >> 1) + 1;
#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#endif

  const double kfac = M_PI * 2 / conf->Lbox;
#if EZMOCK_PT_LOGPK == 1
  const double logkfac = log(kfac);
#endif
  const double Pfac = pow(conf->Lbox, -1.5);

  /* First dimension: compute psi[0]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int i = idx_ij / conf->Ng;
      int j = idx_ij % conf->Ng;
      double ki = (i <= (conf->Ng >> 1)) ? i : i - conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    double ki = (i <= (conf->Ng >> 1)) ? i : i - conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      double kj = (j <= (conf->Ng >> 1)) ? j : j - conf->Ng;
      double kij = ki * ki + kj * kj;
      size_t idx0 = idx_ij * Ngk;

      /* Expand k loop and treat 0 frequency separately. */
      /* k = 0 */
      size_t idx = idx0;
      if (idx == 0) {
        mesh->rhok[0][0] = mesh->rhok[0][1] = 0;
        mesh->rhok2[0][0] = mesh->rhok2[0][1] = 0;
      }
      else {
        double ksq = kij;
#if EZMOCK_PT_LOGPK == 1
        double kmod = log(ksq) * 0.5 + logkfac;         /* log(sqrt(ksq)) */
#else
        double kmod = kfac * sqrt(ksq);                 /* sqrt(ksq) */
#endif
        /* Obtain the power with interpolation. */
        double P = pk_interp(pk, kmod);

#if EZMOCK_PT_LOGPK == 1
        P = exp(P * 0.5);
#else
        P = sqrt(P);
#endif
        P *= Pfac / (ksq * kfac);
        mesh->rhok[idx][0] = -P * mesh->rhok2[idx][1];
        mesh->rhok[idx][1] = P * mesh->rhok2[idx][0];
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;
      }

      /* 0 < k <= k_ny */
      for (int k = 1; k < Ngk; k++) {
        idx = idx0 + k;
        double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
        double kmod = log(ksq) * 0.5 + logkfac;         /* log(sqrt(ksq)) */
#else
        double kmod = kfac * sqrt(ksq);                 /* sqrt(ksq) */
#endif
        /* Obtain the power with interpolation. */
        double P = pk_interp(pk, kmod);

#if EZMOCK_PT_LOGPK == 1
        P = exp(P * 0.5);
#else
        P = sqrt(P);
#endif
        P *= Pfac / (ksq * kfac);
        mesh->rhok[idx][0] = -P * mesh->rhok2[idx][1];
        mesh->rhok[idx][1] = P * mesh->rhok2[idx][0];
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;
      }         /* for k */
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[0]);

  /* Second and third dimension: compute psi[1] and psi[2]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int j = idx_ij % conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      double kj = (j <= (conf->Ng >> 1)) ? j : j - conf->Ng;
      size_t idx0 = idx_ij * Ngk;

      /* 0 <= k <= k_ny */
      for (int k = 0; k < Ngk; k++) {
        size_t idx = idx0 + k;
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * kj;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * kj;
        mesh->rhok[idx][0] *= k;
        mesh->rhok[idx][1] *= k;
      }
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[1]);
  FFT_EXEC_C2R(plan, mesh->rhok, mesh->psi[2]);

  /* Renormalise the displacement fields. */
  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;
  const real norm = 1 / sqrt(num);
  for (int j = 0; j < 3; j++) {
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
    for (size_t i = 0; i < num; i++) mesh->psi[j][i] *= norm;
  }
}


#undef EZMOCK_LOGPK_NAME

#undef EZMOCK_PT_WHITENOISE
#undef EZMOCK_PT_LOGPK

#endif
