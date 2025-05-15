/*******************************************************************************
* errmsg.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "EZmock.h"
#include "errmsg.h"

/******************************************************************************
Function `EZmock_errmsg`:
  Return the error message given the error code.
Arguments:
  * `err`:      the error code.
Return:
  Error message corresponding to the input error code.
******************************************************************************/
const char *EZmock_errmsg(const int err) {
  switch (err) {
    case EZMOCK_SUCCESS:        return "no error";
    case EZMOCK_ERR_MEMORY:     return "failed to allocate memory";
    case EZMOCK_ERR_ARG_ECODE:  return "error code is not set";
    case EZMOCK_ERR_ARG_EZ:     return "EZmock interface uninitialized";
    case EZMOCK_ERR_ARG_COSMO:  return "EZmock cosmology uninitialized";
    case EZMOCK_ERR_ARG_PK:     return "EZmock power spectrum uninitialized";
    case EZMOCK_ERR_ARG_MESH:   return "EZmock density fields uninitialized";
    case EZMOCK_ERR_ARG_LBOX:   return "invalid box size";
    case EZMOCK_ERR_ARG_NGRID:  return "invalid grid size";
    case EZMOCK_ERR_ARG_RNG:    return "invalid random number generator";
    case EZMOCK_ERR_ARG_SEED:   return "invalid random seed";
    case EZMOCK_ERR_ARG_THREAD: return "invalid number of threads";
    case EZMOCK_ERR_ARG_NULL:   return "null input array";
    case EZMOCK_ERR_ARG_NBIN:   return "invalid number of bins";
    case EZMOCK_ERR_ARG_PARAMS: return "EZmock calibration parameters not set";
    case EZMOCK_ERR_ARG_LOWEXP: return "too few expected number of tracers";
    case EZMOCK_ERR_ARG_HIEXP:  return "too many expected number of tracers";
    case EZMOCK_ERR_ARG_NDATA:  return "cannot return the number of tracers";
    case EZMOCK_ERR_ARG_FNAME:  return "invalid file name";
    case EZMOCK_ERR_ARG_RSD:    return "invalid RSD factor";
    case EZMOCK_ERR_COS_GROWTH: return "cannot obtain cosmic growth parameters";
    case EZMOCK_ERR_PK_FINITE:  return "infinite value detected";
    case EZMOCK_ERR_PK_NONPOS:  return "non-positive value detected";
    case EZMOCK_ERR_PK_NONASC:  return "non-ascending k array detected";
    case EZMOCK_ERR_PK_RANGE:   return "k range insufficient for density field";
    case EZMOCK_ERR_PK_NBIN:    return "too few k bins available for the field";
    case EZMOCK_ERR_PK_MODBAO:  return "invalid P(k) after BAO modification";
    case EZMOCK_ERR_RNG_SET:    return "failed to setup the random generator";
    case EZMOCK_ERR_RNG_RESET:  return "random generator reseting failed";
    case EZMOCK_ERR_RNG_JUMP:   return "random generator jumping ahead failed";
    case EZMOCK_ERR_PAR_COSMO:  return "invalid cosmological parameters";
    case EZMOCK_ERR_PAR_RCUT:   return "invalid density threshold parameter";
    case EZMOCK_ERR_PAR_REXP:   return "invalid exponential bias parameter";
    case EZMOCK_ERR_PAR_PBASE:  return "invalid base number for PDF mapping";
    case EZMOCK_ERR_PAR_SIGV:   return "invalid parameter for random motions";
    case EZMOCK_ERR_RHO_USED:   return "the density field has been overwritten";
    case EZMOCK_ERR_RHO_LESS:   return "not enough cells to host all tracers";
    case EZMOCK_ERR_RHO_PBASE:  return "PDF implies too many tracers per cell";
    case EZMOCK_ERR_RHO_NPART:  return "too many particles per cell";
#ifdef OMP
    case EZMOCK_ERR_FFT_INIT:   return "failed to initialise FFTW with OpenMP";
#endif
    case EZMOCK_ERR_FILE:       return "failed to write to the output file";
    case EZMOCK_ERR_UNKNOWN:    return "unknown error";
    default:                    return "unregistered error code";
  }
}
