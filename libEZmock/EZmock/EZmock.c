/*******************************************************************************
* EZmock.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "EZmock.h"
#include "structs.h"
#include "config.h"
#include "errmsg.h"

#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                       Interfaces of the EZmock generator
\*============================================================================*/

/******************************************************************************
Function `EZmock_init`:
  Initialise the EZmock instance.
Arguments:
  * `Lbox`:     side length of the simulation box;
  * `Ngrid`:    number of grids per side for the density field;
  * `randgen`:  random number generation algorithm;
  * `seed`:     random seed;
  * `nthread`:  number of OpenMP threads, 0 for omp_get_max_threads();
  * `err`:      integer storing the error code.
Return:
  Interface of EZmock generator.
******************************************************************************/
EZMOCK *EZmock_init(const double Lbox, const int Ngrid, const int randgen,
    const uint64_t seed, const int nthread, int *err) {
  /* Validate arguments. */
  if (!err || *err != EZMOCK_SUCCESS) return NULL;
  if (Lbox <= 0) {
    *err = EZMOCK_ERR_ARG_LBOX; return NULL;
  }
  if (Ngrid <= 0 || Ngrid > EZMOCK_MAX_GRID_SIZE) {
    *err = EZMOCK_ERR_ARG_NGRID; return NULL;
  }
  switch (randgen) {
    case PRAND_RNG_MRG32K3A:
    case PRAND_RNG_MT19937:
      break;
    default:
      *err = EZMOCK_ERR_ARG_RNG; return NULL;
  }
  if (seed == 0) {
    *err = EZMOCK_ERR_ARG_SEED; return NULL;
  }
  if (nthread < 0) {
    *err = EZMOCK_ERR_ARG_THREAD; return NULL;
  }

  EZMOCK *ez = calloc(1, sizeof(EZMOCK));
  if (!ez) {
    *err = EZMOCK_ERR_MEMORY; return NULL;
  }

  /* Setup anonymous structures. */
  ez->conf = ez->rng = ez->cosmo = ez->pk = ez->par = ez->mesh = NULL;
  if (!(ez->conf = malloc(sizeof(EZMOCK_CONF))) ||
      !(ez->rng = malloc(sizeof(EZMOCK_RNG)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
#ifdef OMP
  conf->nthread = (nthread == 0) ? omp_get_max_threads() : nthread;
#else
  conf->nthread = 1;
#endif
  EZMOCK_RNG *rng = (EZMOCK_RNG *) ez->rng;
  rng->rng = prand_init(randgen, seed, conf->nthread, 0, err);
  if (PRAND_IS_ERROR(*err) || PRAND_IS_WARN(*err)) {
    *err = EZMOCK_ERR_RNG_SET;
    EZmock_destroy(ez); return NULL;
  }
  rng->ran = randgen;
  rng->seed = seed;

  if (!(ez->cosmo = malloc(sizeof(EZMOCK_COSMO)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  cosmo->growth2 = cosmo->vfac = HUGE_VAL;

  if (!(ez->par = malloc(sizeof(EZMOCK_PAR)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }
  EZMOCK_PAR *par = (EZMOCK_PAR *) ez->par;
  par->bao_enhance = par->rho_c = par->rho_exp = HUGE_VAL;
  par->pdf_base = par->sigma_v = HUGE_VAL;

  if (!(ez->mesh = calloc(1, sizeof(EZMOCK_MESH)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }

  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  mesh->rho_replaced = mesh->psi_ref = false;
  mesh->fixamp = mesh->iphase = false;
  mesh->psi[0] = mesh->psi[1] = mesh->psi[2] = NULL;
  mesh->rhok = mesh->rhok2 = NULL;
  mesh->rho = mesh->rhot = NULL;

  conf->Lbox = Lbox;
  conf->Ng = Ngrid;

  return ez;
}


/******************************************************************************
Function `EZmock_destroy`:
  Release memory allocated for the EZmock generator interface.
Arguments:
  * `ez`:       instance of the EZmock generator.
******************************************************************************/
void EZmock_destroy(EZMOCK *ez) {
  if (!ez) return;
  if (ez->conf) free(ez->conf);
  EZmock_rng_destroy(ez->rng);
  if (ez->cosmo) free(ez->cosmo);
  EZmock_pk_destroy(ez->pk);
  if (ez->par) free(ez->par);
  EZmock_mesh_destroy(ez->mesh);
  free(ez);
}
