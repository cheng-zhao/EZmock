/*******************************************************************************
* mangle.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "mangle.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>

/*******************************************************************************
  A simplified implementation of the `polyid` command of mangle:
  https://space.mit.edu/~molly/mangle/
  which reports a polygon that a point lies inside, given a polygon format
  mask file, optionally with the `simple` pixelization scheme.
  Ref: https://doi.org/10.1111/j.1365-2966.2008.13296.x
*******************************************************************************/

/*============================================================================*\
                            Definitions of constants
\*============================================================================*/

#define MANGLE_CHUNK_INIT_SIZE  1048576
#define MANGLE_CHUNK_MAX_SIZE   INT_MAX
#define MANGLE_MAX_RES          15      /* for pixel id < INT_MAX */

#define TWOPI                   0x1.921fb54442d18p+2    /* 2 * PI */
#define DEG2RAD                 0x1.1df46a2529d39p-6    /* PI / 180 */

/*============================================================================*\
                           Definitions of error codes
\*============================================================================*/

#define MANGLE_ERR_ARGS         (-1)
#define MANGLE_ERR_MEMORY       (-2)
#define MANGLE_ERR_FILE         (-3)
#define MANGLE_ERR_FORMAT       (-4)
#define MANGLE_ERR_PIXEL        (-5)
#define MANGLE_ERR_RES          (-6)
#define MANGLE_ERR_POLYGON      (-7)
#define MANGLE_ERR_PID          (-8)
#define MANGLE_ERR_CAP          (-9)
#define MANGLE_ERR_NPOLY_MORE   (-10)
#define MANGLE_ERR_NPOLY_LESS   (-11)
#define MANGLE_ERR_NOPOLY       (-12)


/*============================================================================*\
                        Definitions for sorting polygons
\*============================================================================*/

/******************************************************************************
Function `tim_sort`:
  Sort two arrays with predefined data types and comparison rules,
  using the timsort algorithm.
Arguments:
  * `x`:        pointer to the first array;
  * `y`:        pointer to the second array;
  * `size`:     size of the input arrays.
******************************************************************************/

/* Data type for the two arrays. */
#define TIMSORT_DTYPE1                  POLYGON
#define TIMSORT_DTYPE2                  POLYGON
/* Data type of the variable for binding elements from both arrays. */
#define TIMSORT_BIND_DTYPE              POLYGON
/* Macro for changing the starting indices of the arrays. */
#define TIMSORT_OFFSET(x,y,s)           (x) += (s);
/* Macro for comparing array elements with indices i and j. */
#define TIMSORT_CMP_IDX(x,y,i,j)        ((x)[(i)].pixel - (x)[(j)].pixel)
/* Macro for comparing binded value with the arrays with index i. */
#define TIMSORT_CMP_BIND(x,y,i,b)       ((x)[(i)].pixel - (b).pixel)
/* Macro for assigning values with index i to those with index j. */
#define TIMSORT_ASSIGN_IDX(x,y,i,j)                                     \
  memcpy((x) + (j), (x) + (i), sizeof(POLYGON))
/* Macro for assigning array elements with index i, to the binded variable. */
#define TIMSORT_ASSIGN_BIND(x,y,i,b)    memcpy(&(b), (x) + (i), sizeof(POLYGON))
/* Macro for assigning binded value to array elements with index i. */
#define TIMSORT_GET_BIND(x,y,i,b)       memcpy((x) + (i), &(b), sizeof(POLYGON))
/* Macro for swapping array elements with indices i and j. */
#define TIMSORT_SWAP(x,y,i,j) {                                         \
  memcpy((y), (x) + (i), sizeof(POLYGON));                              \
  memcpy((x) + (i), (x) + (j), sizeof(POLYGON));                        \
  memcpy((x) + (j), (y), sizeof(POLYGON));                              \
}

/* Note that the second pointer serves as the buffer for swapping elements. */
#include "timsort.c"


/*============================================================================*\
                      Functions for initialising polygons
\*============================================================================*/

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of a chunk.
Arguments:
  * `chunk`:    address of the chunk;
  * `size`:     size of the chunk.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int chunk_resize(char **chunk, size_t *size) {
  if (MANGLE_CHUNK_MAX_SIZE / 2 < *size) return MANGLE_ERR_FILE;
  size_t num = *size << 1;
  char *tmp = realloc(*chunk, num * sizeof(char));
  if (!tmp) return MANGLE_ERR_MEMORY;
  *chunk = tmp;
  *size = num;
  return 0;
}

/******************************************************************************
Function `mangle_read`:
  Read polygons from a mask file.
Arguments:
  * `fname`:    name of the mangle polygon-format mask file;
  * `wmin`:     minimum weight of polygons to be kept;
  * `num`:      number of valid polygons read from file;
  * `reso`:     resolution of pixelization read from file;
  * `err`:      an integer for indicating error messages.
Return:
  Address of the array for valid polygons.
******************************************************************************/
static POLYGON *mangle_read(const char *fname, const double wmin, int *num,
    int *reso, int *err) {
  /* Allocate memory for reading the file by chunks. */
  size_t csize = MANGLE_CHUNK_INIT_SIZE;
  char *chunk = malloc(csize * sizeof(char));
  if (!chunk) {
    *err = MANGLE_ERR_MEMORY;
    return NULL;
  }

  FILE *fp = fopen(fname, "r");
  if (!fp) {
    *err = MANGLE_ERR_FILE;
    free(chunk);
    return NULL;
  }

  int npoly, res;
  npoly = res = 0;
  bool header_done = false;
  size_t nread, nrest;

  /* Start reading the file (header) by chunks. */
  nrest = 0;
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == '\0') {                 /* omit empty lines */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (*p >= '0' && *p <= '9') {             /* start with a number */
        if (npoly || sscanf(p, "%d polygons", &npoly) != 1)
          *err = MANGLE_ERR_FORMAT;
        else if (npoly <= 0) *err = MANGLE_ERR_NPOLY_LESS;
      }
      else if (strncmp(p, "pixelization", 12) == 0) {   /* pixelization */
        if (sscanf(p, "pixelization %ds", &res) != 1) *err = MANGLE_ERR_PIXEL;
        else if (res <= 0 || res > MANGLE_MAX_RES)  *err = MANGLE_ERR_RES;
      }
      else if (strncmp(p, "polygon", 7) == 0) { /* polygon data starts */
        *endl = '\n';
        nrest = end - p;
        memmove(chunk, p, nrest);
        header_done = true;                     /* finish reading header */
        break;
      }

      if (*err) {
        free(chunk); fclose(fp);
        return NULL;
      }
      /* Continue with the next line. */
      p = endl + 1;
    }
    if (header_done) break;

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      *err = chunk_resize(&chunk, &csize);
      if (*err) {
        free(chunk); fclose(fp);
        return NULL;
      }
      nrest += nread;
      continue;
    }
    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  /* Allocate memory for polygons. */
  POLYGON *poly = NULL;
  if (!npoly) *err = MANGLE_ERR_FORMAT;
  else if (!(poly = malloc(npoly * sizeof(POLYGON)))) *err = MANGLE_ERR_MEMORY;
  if (*err) {
    free(chunk); fclose(fp);
    return NULL;
  }
  for (int i = 0; i < npoly; i++) poly[i].cap = NULL;

  /* Minimum and maximum pixel id with the given resolution. */
  const int pmin = ((1 << (res << 1)) - 1) / 3;         /* (4^n - 1) / 3 */
  const int pmax = pmin << 2;

  int ipoly, icap, iread;
  ipoly = icap = iread = 0;
  int skip_ncap = 0;

  /* Continue reading the polygons. */
  nread = 0;
  do {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread && nread < csize - nrest) *end = '\n';    /* '\n' for last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == '\0') {                 /* omit empty lines */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (!icap) {      /* icap == 0: read polygon */
        skip_ncap = 0;
        if (iread >= npoly) *err = MANGLE_ERR_NPOLY_MORE;
        else {
          POLYGON *ply = poly + ipoly;
          if (sscanf(p, "polygon %d ( %d caps, %lf weight, %d pixel, %lf",
              &ply->polyid, &ply->ncap, &ply->weight, &ply->pixel, &ply->area)
              != 5) {
            *err = MANGLE_ERR_POLYGON;
          }

          iread++;
          if (ply->ncap == 0) {
            p = endl + 1;
            continue;
          }
          if (ply->pixel < pmin || ply->pixel > pmax) *err = MANGLE_ERR_PID;
          if (ply->weight < wmin) {             /* skip the polygon */
            if ((skip_ncap = ply->ncap)) icap++;
            p = endl + 1;
            continue;
          }
          /* Allocate memory for caps. */
          if (!(ply->cap = malloc(ply->ncap * sizeof(POLYCAP))))
            *err = MANGLE_ERR_MEMORY;
          icap++;
          ipoly++;
        }
      }
      else {            /* icap != 0: read caps of the polygon */
        POLYCAP cap;
        if (sscanf(p, "%lf %lf %lf %lf", cap, cap + 1, cap + 2, cap + 3)
            != 4) {
          *err = MANGLE_ERR_CAP;
        }
        else if (skip_ncap) {
          if (icap++ >= skip_ncap) icap = 0;
        }
        else {
          memcpy(poly[ipoly - 1].cap + icap - 1, cap, sizeof(POLYCAP));
          if (icap++ >= poly[ipoly - 1].ncap) icap = 0;
        }
      }

      /* Check errors. */
      if (*err) {
        free(chunk); fclose(fp);
        for (int i = 0; i < ipoly; i++) if (poly[i].cap) free(poly[i].cap);
        free(poly);
        return NULL;
      }
      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      *err = chunk_resize(&chunk, &csize);
      if (*err) {
        free(chunk); fclose(fp);
        for (int i = 0; i < ipoly; i++) if (poly[i].cap) free(poly[i].cap);
        free(poly);
        return NULL;
      }
      nrest += nread;
      continue;
    }
    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp)));

  free(chunk);
  fclose(fp);

  if (iread != npoly) *err = MANGLE_ERR_NPOLY_LESS;
  else if (!ipoly) *err = MANGLE_ERR_NOPOLY;
  if (*err) {
    for (int i = 0; i < ipoly; i++) if (poly[i].cap) free(poly[i].cap);
    free(poly);
    return NULL;
  }
  /* Reduce memory cost if applicable. */
  if (ipoly < iread) {
    POLYGON *tmp = realloc(poly, ipoly * sizeof(POLYGON));
    if (tmp) poly = tmp;
  }

  *num = ipoly;
  *reso = res;
  return poly;
}

/******************************************************************************
Function `mangle_pix_idx`:
  Create indices for the polygon pixels.
Arguments:
  * `poly`:     pointer to the array for polygons;
  * `npoly`:    number of polygons in the array;
  * `res`:      resolution of the pixelization;
  * `err`:      an integer for indicating error messages.
Return:
  Address of the array for pixel ID indices.
******************************************************************************/
static int *mangle_pix_idx(POLYGON *poly, const int npoly, const int res,
    int *err) {
  /* Compute the number of pixels given the resolution, and allocate memory. */
  const int npix = 1 << (res << 1);
  const int pmin = (npix - 1) / 3;      /* minimum pixel ID */
  int *pix = calloc(npix + 1, sizeof(int));
  if (!pix) {
    *err = MANGLE_ERR_MEMORY;
    return NULL;
  }

  /* Sort polygons based on the pixel IDs. */
  POLYGON tmp;
  tim_sort(poly, &tmp, npoly);

  /* Record indices. */
  int prev = -1;
  for (int i = 0; i < npoly; i++) {
    int j = poly[i].pixel - pmin;
    if (j != prev) {
      for (int k = prev + 1; k <= j; k++) pix[k] = i;
      prev = j;
    }
  }
  pix[npix] = npoly;

  return pix;
}


/*============================================================================*\
                          Functions for polygon query
\*============================================================================*/

/******************************************************************************
Function `mangle_inside_cap`:
  Check if the endpoint of a unit vector is inside a cap.
Arguments:
  * `cap`:      pointer to the cap to be checked;
  * `v`:        the unit vector to be checked.
Return:
  True if the endpoint is inside the cap.
******************************************************************************/
static inline bool mangle_inside_cap(const POLYCAP *cap, const double *v) {
  /* Cosine of the angle between the unit vector and the polar axis. */
  double angle = (*cap)[0] * v[0] + (*cap)[1] * v[1] + (*cap)[2] * v[2];
  if ((*cap)[3] >= 0) return angle > 1 - (*cap)[3];
  else return angle < 1 + (*cap)[3];
}

/******************************************************************************
Function `mangle_inside_poly`:
  Check if the endpoint of a unit vector is inside a polygon.
Arguments:
  * `poly`:     pointer to the polygon to be checked;
  * `v`:        the unit vector to be checked.
Return:
  True if the endpoint is inside the polygon.
******************************************************************************/
static inline bool mangle_inside_poly(const POLYGON *poly, const double *v) {
  /* Visit all caps of the polygon. */
  for (int i = 0; i < poly->ncap; i++) {
    if (!mangle_inside_cap(poly->cap + i, v)) return false;
  }
  return true;
}

/******************************************************************************
Function `mangle_query_pix`:
  Find the pixel index of a point, following the `simple` pixelization scheme.
Arguments:
  * `res`:      resolution for the pixelization;
  * `az`:       azimuth angle (in radians) of the point;
  * `el`:       elevation angle (in radians) of the point.
Return:
  Index of the pixel (starting from 0).
******************************************************************************/
static inline int mangle_query_pix(const int res, const double az,
    const double el) {
  const int nside = 1 << res;           /* number of pixels per side: 2^res */
  const int m = (az >= TWOPI) ? nside - 1 : (az / TWOPI) * nside;
  double sine = sin(el);
  const int n = (sine == 1) ? 0 : ceil((1 - sine) * 0.5 * nside) - 1;
  return (n << res) + m;
}

/******************************************************************************
Function `mangle_query_poly`:
  Find a polygon that contains a given point.
Arguments:
  * `poly`:     pointer to the array of polygons to be visited;
  * `npoly`:    number of polygons to be visited;
  * `az`:       azimuth angle (in radians) of the point;
  * `el`:       elevation angle (in radians) of the point.
Return:
  Pointer to the polygon that contains the point.
******************************************************************************/
static inline POLYGON *mangle_query_poly(POLYGON *poly, const int npoly,
    const double az, const double el) {
  double v[3];          /* unit vector given the angular direction */
  v[0] = cos(el) * cos(az);
  v[1] = cos(el) * sin(az);
  v[2] = sin(el);

  /* Visit the polygons and report the first match. */
  for (int i = 0; i < npoly; i++) {
    if (mangle_inside_poly(poly + i, v)) return poly + i;
  }

  return NULL;
}


/*============================================================================*\
                     Interfaces for applying polygon masks
\*============================================================================*/

/******************************************************************************
Function `mangle_init`:
  Initialise a mangle polygon mask.
Arguments:
  * `fname`:    name of the mangle polygon-format mask file;
  * `wmin`:     minimum weight of polygons to be kept;
  * `err`:      an integer for indicating error messages.
Return:
  Address of the structure for the mask.
******************************************************************************/
MANGLE *mangle_init(const char *fname, const double wmin, int *err) {
  if (!fname || !(*fname)) {
    if (err) *err = MANGLE_ERR_ARGS;
    return NULL;
  }

  /* Allocate memory. */
  MANGLE *mask = calloc(1, sizeof *mask);
  if (!mask) {
    *err = MANGLE_ERR_MEMORY;
    return NULL;
  }

  *err = 0;
  mask->poly = NULL;
  mask->pix = NULL;

  /* Read polygons from file. */
  if (!(mask->poly = mangle_read(fname, wmin, &mask->npoly, &mask->res, err))) {
    mangle_destroy(mask);
    return NULL;
  }

  if (!(mask->pix = mangle_pix_idx(mask->poly, mask->npoly, mask->res, err))) {
    mangle_destroy(mask);
    return NULL;
  }

  return mask;
}

/******************************************************************************
Function `mangle_destroy`:
  Deconstruct the structure for mangle polygon mask.
Arguments:
  * `mask`:     address of the structure for the mask.
******************************************************************************/
void mangle_destroy(MANGLE *mask) {
  if (!mask) return;
  if (mask->poly) {
    for (int i = 0; i < mask->npoly; i++) {
      if (mask->poly[i].cap) free(mask->poly[i].cap);
    }
    free(mask->poly);
  }
  if (mask->pix) free(mask->pix);
  free(mask);
}

/******************************************************************************
Function `mangle_query`:
  Find a polygon that contains a given point.
Arguments:
  * `mask`:     address of the structure for the mask;
  * `ra`:       right ascension (in degrees) of the point;
  * `dec`:      declination (in degrees) of the point.
Return:
  Pointer to the polygon that contains the point; NULL if no polygon is found.
******************************************************************************/
POLYGON *mangle_query(const MANGLE *mask, const double ra, const double dec) {
  /* Compute the azimuth and elevation angles in radians. */
  const double az = ra * DEG2RAD;
  const double el = dec * DEG2RAD;
  POLYGON *poly = mask->poly;
  int npoly = mask->npoly;

  /* Find the pixel index if applicable. */
  if (mask->res) {
    int idx = mangle_query_pix(mask->res, az, el);
    poly += mask->pix[idx];
    npoly = mask->pix[idx + 1] - mask->pix[idx];
  }

  /* Find the polygon that contains the point. */
  return mangle_query_poly(poly, npoly, az, el);
}

/******************************************************************************
Function `mangle_errmsg`:
  Produce error message for a given error code.
Arguments:
  * `err`:      the error code
Return:
  A string with error message.
******************************************************************************/
char *mangle_errmsg(const int err) {
  switch (err) {
    case 0:
      return "no error is detected";
    case MANGLE_ERR_ARGS:
      return "invalid argument for mask initialization";
    case MANGLE_ERR_MEMORY:
      return "failed to allocate memory for the mask";
    case MANGLE_ERR_FILE:
      return "cannot read the polygon mask file";
    case MANGLE_ERR_FORMAT:
      return "invalid polygon mask file: unexpected format";
    case MANGLE_ERR_PIXEL:
      return "invalid polygon mask file: only the simple pixelization scheme "
          "is supported";
    case MANGLE_ERR_RES:
      return "invalid polygon mask file: unexpected resolution";
    case MANGLE_ERR_POLYGON:
      return "invalid polygon mask file: unexpected polygon format";
    case MANGLE_ERR_PID:
      return "invalid polygon mask file: unexpected pixel id";
    case MANGLE_ERR_CAP:
      return "invalid polygon mask file: unexpected definition of caps";
    case MANGLE_ERR_NPOLY_MORE:
      return "invalid polygon mask file: too many polygons";
    case MANGLE_ERR_NPOLY_LESS:
      return "invalid polygon mask file: too few polygons";
    case MANGLE_ERR_NOPOLY:
      return "no valid polygon read from the file";
    default:
      return "unknown problem";
  }
}
