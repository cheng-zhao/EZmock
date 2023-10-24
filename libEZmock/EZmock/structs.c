/*******************************************************************************
* structs.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "structs.h"
#include "config.h"
#include <stdlib.h>

/*============================================================================*\
                  Functions for deconstructing data structures
\*============================================================================*/

/******************************************************************************
Function `EZmock_rng_destroy`:
  Release memory allocated for the EZmock random number generator.
Arguments:
  * `rng`:      instance of the random number generator.
******************************************************************************/
void EZmock_rng_destroy(void *rng) {
  if (!rng) return;
  EZMOCK_RNG *p = (EZMOCK_RNG *) rng;
  prand_destroy(p->rng);
  free(rng);
}

/******************************************************************************
Function `EZmock_pk_destroy`:
  Release memory allocated for the EZmock linear power spectrum.
Arguments:
  * `pk`:       instance of the linear power spectrum.
******************************************************************************/
void EZmock_pk_destroy(void *pk) {
  if (!pk) return;
  EZMOCK_PK *p = (EZMOCK_PK *) pk;
  if (p->k) free(p->k);
  if (p->P) free(p->P);
  if (p->P2) free(p->P2);
  free(pk);
}

/******************************************************************************
Function `EZmock_mesh_destroy`:
  Release memory allocated for the EZmock density and displacement fields.
Arguments:
  * `mesh`:     instance of the EZmock fields.
******************************************************************************/
void EZmock_mesh_destroy(void *mesh) {
  if (!mesh) return;
  EZMOCK_MESH *p = (EZMOCK_MESH *) mesh;
  if (!p->psi_ref) {            /* free only if psis are not references */
    if (p->psi[0]) FFT_FREE(p->psi[0]);
    if (p->psi[1]) FFT_FREE(p->psi[1]);
    if (p->psi[2]) FFT_FREE(p->psi[2]);
  }
  /* `rho` and `rhot` may recycle `rhok` and `rhok2`, respectively. */
  if (p->rho && p->rho != (FFT_REAL *) p->rhok) FFT_FREE(p->rho);
  if (p->rhot && p->rhot != (FFT_REAL *) p->rhok2) FFT_FREE(p->rhot);
  if (p->rhok) FFT_FREE(p->rhok);
  if (p->rhok2) FFT_FREE(p->rhok2);
  free(mesh);
}

/******************************************************************************
Function `EZmock_pdf_init`:
  Initialise the data structure for the EZmock tracer PDF.
return:
  Instance of the EZmock tracer PDF.
******************************************************************************/
EZMOCK_PDF *EZmock_pdf_init(void) {
  EZMOCK_PDF *pdf = calloc(1, sizeof(EZMOCK_PDF));
  if (!pdf) return NULL;
  pdf->hist = calloc(EZMOCK_TRACER_DENSITY_NUM_LOGB_BIN, sizeof(size_t));
  if (!pdf->hist) {
    free(pdf); return NULL;
  }
  pdf->pdf = NULL;
  return pdf;
}

/******************************************************************************
Function `EZmock_pdf_destroy`:
  Release memory allocated for the EZmock tracer PDF.
Arguments:
  * `pdf`:      instance of the EZmock tracer PDF.
******************************************************************************/
void EZmock_pdf_destroy(EZMOCK_PDF *pdf) {
  if (!pdf) return;
  if (pdf->pdf) free(pdf->pdf);
  if (pdf->hist) free(pdf->hist);
  free(pdf);
}

