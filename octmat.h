#ifndef OCTMAT_H
#define OCTMAT_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef OCTMAT_ALIGN
#define OCTMAT_ALIGN 16
#endif

typedef struct {
  uint8_t *data;
  uint16_t rows;
  uint16_t cols;
  uint16_t cols_al;
} octmat;

typedef struct {
  uint16_t row;
  uint16_t col;
  uint8_t val;
} coo_tuple;

typedef struct {
  uint8_t *vals;
  uint16_t *idxs;
  uint16_t nz;
  uint16_t nz_al;
} crsrow;

typedef struct {
  crsrow *row;
  uint16_t rows;
  uint16_t cols;
} crsmat;

#define OM_INITIAL                                                             \
  { .rows = 0, .cols = 0, .cols_al = 0, .data = 0 }
#define CRS_INITIAL                                                            \
  { .row = 0, .rows = 0 }

#define om_R(v, x) ((v).data + ((x) * (v).cols_al))
#define om_P(v) om_R(v, 0)
#define om_A(v, x, y) (om_R(v, x)[(y)])

void om_resize(octmat *v, uint16_t r, uint16_t c);
void om_copy(octmat *v1, octmat *v0);
void om_destroy(octmat *v);
void om_print(FILE *stream, octmat m);

int crsrow_find(uint16_t *a, int L, int R, uint16_t x);
uint16_t crsrow_resize(crsrow *a, uint16_t nz);
uint16_t crsrow_unpack(crsrow *a, uint8_t *dense);
uint16_t crsrow_pack(crsrow *a, uint8_t *dense, uint16_t cols);

void crs_resize(crsmat *v, uint16_t r, uint16_t c);
void crs_copy(crsmat *v1, crsmat *v0);
void crs_destroy(crsmat *v);
void crs_print(FILE *stream, crsmat *m);

#endif
