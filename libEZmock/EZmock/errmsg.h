/*******************************************************************************
* errmsg.h: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <https://www.gnu.org/licenses/>.

*******************************************************************************/

#ifndef __ERRMSG_H__
#define __ERRMSG_H__

/*============================================================================*\
                           Definitions of error codes
\*============================================================================*/

typedef enum {
  EZMOCK_SUCCESS = 0,

  EZMOCK_ERR_MEMORY,

  EZMOCK_ERR_ARG_ECODE,
  EZMOCK_ERR_ARG_EZ,
  EZMOCK_ERR_ARG_COSMO,
  EZMOCK_ERR_ARG_PK,
  EZMOCK_ERR_ARG_MESH,
  EZMOCK_ERR_ARG_LBOX,
  EZMOCK_ERR_ARG_NGRID,
  EZMOCK_ERR_ARG_RNG,
  EZMOCK_ERR_ARG_SEED,
  EZMOCK_ERR_ARG_THREAD,
  EZMOCK_ERR_ARG_NULL,
  EZMOCK_ERR_ARG_NBIN,
  EZMOCK_ERR_ARG_PARAMS,
  EZMOCK_ERR_ARG_LOWEXP,
  EZMOCK_ERR_ARG_HIEXP,
  EZMOCK_ERR_ARG_NDATA,

  EZMOCK_ERR_COS_GROWTH,

  EZMOCK_ERR_PK_FINITE,
  EZMOCK_ERR_PK_NONPOS,
  EZMOCK_ERR_PK_NONASC,
  EZMOCK_ERR_PK_RANGE,
  EZMOCK_ERR_PK_NBIN,
  EZMOCK_ERR_PK_MODBAO,

  EZMOCK_ERR_RNG_SET,
  EZMOCK_ERR_RNG_RESET,
  EZMOCK_ERR_RNG_JUMP,

  EZMOCK_ERR_PAR_COSMO,
  EZMOCK_ERR_PAR_RCUT,
  EZMOCK_ERR_PAR_REXP,
  EZMOCK_ERR_PAR_PBASE,
  EZMOCK_ERR_PAR_SIGV,

  EZMOCK_ERR_RHO_USED,
  EZMOCK_ERR_RHO_LESS,
  EZMOCK_ERR_RHO_PBASE,
  EZMOCK_ERR_RHO_NPART,

#ifdef OMP
  EZMOCK_ERR_FFT_INIT,
#endif

  EZMOCK_ERR_UNKNOWN,
} EZMOCK_ERRNO;


#endif

