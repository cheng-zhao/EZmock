/*******************************************************************************
* run_mock.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "read_pk.h"
#include "EZmock.h"
#include "structs.h"
#include "save_res.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef OMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return EZMOCK_ERR_CFG;
  }

  int err = 0;
  printf("Initializing the EZmock generator ...");
  fflush(stdout);

  EZMOCK *ez;
  if (!(ez = EZmock_init(conf->Lbox, conf->Ngrid, conf->rng, conf->seed,
      conf->nthread, &err))) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to initialize the EZmock generator\n");
    conf_destroy(conf);
    return EZMOCK_ERR_INIT;
  }

  printf(FMT_DONE);
  printf("Setting cosmological parameters for EZmock generation ...");
  fflush(stdout);

  if (EZmock_set_cosmology(ez, conf->growth2, conf->vfac, conf->eval_growth,
      conf->redshift, conf->zpk, conf->omega_m, conf->omega_nu, conf->eos_w,
      &err)) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to set cosmological parameters for EZmock generation\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_COSMO;
  }

  printf(FMT_DONE);

  PK *pk = NULL;
  if (!(pk = read_pk(conf)) ||
      EZmock_setup_linear_pk(ez, pk->k, pk->n, pk->Plin, pk->Pnw, conf->bao_mod,
      conf->logint, &err)) {
    if (err) P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to setup the input power spectrum\n");
    conf_destroy(conf); EZmock_destroy(ez); pk_destroy(pk);
    return EZMOCK_ERR_PK;
  }
  pk_destroy(pk);

  printf("Generating the EZmock density field ...");
  fflush(stdout);

  if (EZmock_create_dens_field(ez, NULL, false, NULL, conf->fixamp, conf->iphase,
      &err)) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to construct the density field\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_RHO;
  }

  printf(FMT_DONE);
  printf("Populating EZmock tracers ...");
  fflush(stdout);

  const real params[] = {
    conf->rho_c,
    conf->rho_exp,
    conf->pdf_base,
    conf->sigv
  };

  size_t ntracer;
  real *x, *y, *z, *vx, *vy, *vz;

  if (EZmock_populate_tracer(ez, params, conf->Ngal, conf->particle, &ntracer,
      &x, &y, &z, &vx, &vy, &vz, &err)) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to sample EZmock tracers\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_RHO;
  }

  printf(FMT_DONE);
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  conf->growth2 = cosmo->growth2;
  conf->vfac = cosmo->vfac;
  EZmock_destroy(ez);

  if (save_box(conf, x, y, z, vx, vy, vz, ntracer)) {
    printf(FMT_FAIL);
    P_EXT("failed to save the output tracer catalog\n");
    conf_destroy(conf);
    free(x); free(y); free(z);
    free(vx); free(vy); free(vz);
    return EZMOCK_ERR_SAVE;
  }

  conf_destroy(conf);
  free(x); free(y); free(z);
  free(vx); free(vy); free(vz);
  return 0;
}
