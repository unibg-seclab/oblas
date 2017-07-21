#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "oblas.h"
#include "octmat.h"

void om_resize(octmat *v, uint16_t r, uint16_t c) {
  v->rows = r;
  v->cols = c;

  void *aligned = NULL;
  v->cols_al = (c / OCTMAT_ALIGN + ((c % OCTMAT_ALIGN) ? 1 : 0)) * OCTMAT_ALIGN;

  if (posix_memalign(&aligned, OCTMAT_ALIGN, r * v->cols_al) != 0) {
    exit(ENOMEM);
  }
  ozero(aligned, 0, v->cols_al * r);
  v->data = (uint8_t *)aligned;
}

void om_copy(octmat *v1, octmat *v0) {
  v1->rows = v0->rows;
  v1->cols = v0->cols;
  v1->cols_al = v0->cols_al;

  if (!v1->data) {
    void *aligned = NULL;
    if (posix_memalign(&aligned, OCTMAT_ALIGN, v0->rows * v0->cols_al) != 0)
      exit(ENOMEM);
    v1->data = (uint8_t *)aligned;
  }
  memcpy(v1->data, v0->data, v0->rows * v0->cols_al);
}

void om_destroy(octmat *v) {
  v->rows = 0;
  v->cols = 0;
  v->cols_al = 0;
  free(v->data);
  v->data = NULL;
}

void om_print(FILE *stream, octmat m) {
  fprintf(stream, "[%ux%u]\n", m.rows, m.cols);
  for (int i = 0; i < m.rows; i++) {
    fprintf(stream, "|%3d", om_A(m, i, 0));
    for (int j = 1; j < m.cols; j++) {
      fprintf(stream, ", %3d", om_A(m, i, j));
    }
    fprintf(stream, "  |\n");
  }
}

int crsrow_find(uint16_t *a, int L, int R, uint16_t x) {
  size_t idx;
  while (L <= R) {
    idx = (L + R) / 2;
    if (a[idx] == x)
      return idx;
    else if (a[idx] < x)
      L = idx + 1;
    else if (a[idx] > x)
      R = idx - 1;
  }
  return -1;
}

uint16_t crsrow_resize(crsrow *a, uint16_t nz) {
  if (a->nz_al < nz) {
    a->nz_al = (nz / OCTMAT_ALIGN + ((nz % OCTMAT_ALIGN) ? 1 : 0)) * OCTMAT_ALIGN;

    if (1) {
      void *aligned = NULL;
      if (posix_memalign(&aligned, OCTMAT_ALIGN, sizeof(uint16_t) * a->nz_al) != 0)
        exit(ENOMEM);
      ozero(aligned, 0, sizeof(uint16_t) * a->nz_al);
      a->idxs = (uint16_t *)aligned;
    }

    if (1) {
      void *aligned = NULL;
      if (posix_memalign(&aligned, OCTMAT_ALIGN, a->nz_al) != 0)
        exit(ENOMEM);
      ozero(aligned, 0, a->nz_al);
      a->vals = (uint8_t *)aligned;
    }
  }
  a->nz = nz;

  return nz;
}

void crsrow_copy(crsrow *v1, crsrow *v0) {
  v1->nz = v0->nz;
  v1->nz_al = v0->nz_al;

  if (!v1->idxs) {
    void *aligned = NULL;
    if (posix_memalign(&aligned, OCTMAT_ALIGN, sizeof(uint16_t) * v0->nz_al) != 0)
      exit(ENOMEM);
    v1->idxs = (uint16_t *)aligned;
  }

  if (!v1->vals) {
    void *aligned = NULL;
    if (posix_memalign(&aligned, OCTMAT_ALIGN, v0->nz_al) != 0)
      exit(ENOMEM);
    v1->vals = (uint8_t *)aligned;
  }

  memcpy(v1->idxs, v0->idxs, v0->nz_al);
  memcpy(v1->vals, v0->vals, v0->nz_al);
}

void crs_freerow(crsrow *a) {
  free(a->idxs);
  free(a->vals);
  a->nz = 0;
  a->nz_al = 0;
}

uint16_t crsrow_unpack(crsrow *a, uint8_t *dense) {
  for (int idx = 0; idx < a->nz; idx++) {
    dense[a->idxs[idx]] = a->vals[idx];
  }
  return a->nz;
}

uint16_t crsrow_pack(crsrow *a, uint8_t *dense, uint16_t cols) {
  int nz = 0;
  for (int idx = 0; idx < cols; idx++) {
    nz += (dense[idx] != 0);
  }

  crsrow_resize(a, nz);

  int at = 0;
  for (int idx = 0; idx < cols; idx++) {
    if (dense[idx] != 0) {
      a->idxs[at] = idx;
      a->vals[at] = dense[idx];
      at++;
    }
  }

  return nz;
}

void crs_resize(crsmat *v, uint16_t r, uint16_t c) {
  v->row = calloc(r, sizeof(crsrow));
  v->rows = r;
}

void crs_copy(crsmat *v1, crsmat *v0) {
}

void crs_destroy(crsmat *v) {
  for(int r = 0; r < v->rows; r++) {
    crs_freerow(&v->row[r]);
  }
  free(v->row);
}

void crs_print(FILE *stream, crsmat *m) {
}

