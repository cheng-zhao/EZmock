/*******************************************************************************
* read_ascii.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <ctype.h>

/*============================================================================*\
                      Functions for reading file by chunks
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
  /* Assume the arguments are not NULL. */
  size_t num;
  if (!(*chunk)) num = EZMOCK_FILE_CHUNK;
  else {
    if (EZMOCK_MAX_CHUNK / 2 < *size) return EZMOCK_ERR_FILE;
    num = *size << 1;
  }

  char *tmp = realloc(*chunk, num * sizeof(char));
  if (!tmp) return EZMOCK_ERR_MEMORY;

  *chunk = tmp;
  *size = num;
  return 0;
}

/******************************************************************************
Function `read_ascii_table`:
  Read the first two columns of an ASCII file as double arrays.
Arguments:
  * `fname`:    filename of the input catalog;
  * `x`:        array for the first column;
  * `y`:        array for the second column;
  * `num`:      number of lines read successfully;
  * `verb`:     indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii_table(const char *fname, double **x, double **y, size_t *num,
    const int verb) {
  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    return EZMOCK_ERR_FILE;
  }

  if (verb) printf("  Reading table from file: %s\n", fname);

  /* Prepare for the chunk. */
  char *chunk = NULL;
  size_t csize = 0;
  if (chunk_resize(&chunk, &csize)) {
    P_ERR("failed to allocate memory for reading the file by chunks.\n");
    fclose(fp);
    return EZMOCK_ERR_MEMORY;
  }

  /* Allocate memory for the data. */
  size_t max = EZMOCK_DATA_INIT_NUM;
  double *nx = malloc(max * sizeof(double));
  double *ny = malloc(max * sizeof(double));
  if (!nx || !ny) {
    P_ERR("failed to allocate memory for the table.\n");
    fclose(fp); free(chunk);
    if (nx) free(nx);
    if (ny) free(ny);
    return EZMOCK_ERR_MEMORY;
  }

  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), csize - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < csize - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == EZMOCK_READ_COMMENT || *p == '\0') {   /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      if (sscanf(p, "%lf %lf", nx + n, ny + n) != 2) {
        P_ERR("failed to read line: %s\n", p);
        fclose(fp); free(chunk); free(nx); free(ny);
        return EZMOCK_ERR_FILE;
      }

      /* Enlarge the memory for the data if necessary. */
      if (++n >= max) {
        if (SIZE_MAX / 2 < max) {
          P_ERR("too many lines in the file: `%s'.\n", fname);
          fclose(fp); free(chunk); free(nx); free(ny);
          return EZMOCK_ERR_FILE;
        }
        max <<= 1;
        double *tmp = realloc(nx, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the samples.\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return EZMOCK_ERR_MEMORY;
        }
        nx = tmp;
        tmp = realloc(ny, sizeof(double) * max);
        if (!tmp) {
          P_ERR("failed to allocate memory for the table.\n");
          fclose(fp); free(chunk); free(nx); free(ny);
          return EZMOCK_ERR_MEMORY;
        }
        ny = tmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      if (chunk_resize(&chunk, &csize)) {
        P_ERR("failed to allocate memory for reading the file by chunk.\n");
        fclose(fp); free(chunk); free(nx); free(ny);
        return EZMOCK_ERR_MEMORY;
      }
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'.\n", fname);
    fclose(fp); free(chunk); free(nx); free(ny);
    return EZMOCK_ERR_FILE;
  }

  free(chunk);
  if (fclose(fp)) P_WRN("failed to close file: `%s'.\n", fname);
  if (verb) printf("  Number of samples: %zu\n", n);

  *x = nx;
  *y = ny;
  *num = n;
  return 0;
}

