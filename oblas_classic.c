#include "oblas.h"
#include "octmul_hilo.h"

void ocopy(uint8_t *restrict a, uint8_t *restrict b, uint16_t i, uint16_t j,
           uint16_t k) {
  octet *ap = a + (i * ALIGNED_COLS(k));
  octet *bp = b + (j * ALIGNED_COLS(k));

  for (int idx = 0; idx < k; idx++) {
    ap[idx] = bp[idx];
  }
}

void oswaprow(uint8_t *restrict a, uint16_t i, uint16_t j, uint16_t k) {
  if (i == j)
    return;
  octet *ap = a + (i * ALIGNED_COLS(k));
  octet *bp = a + (j * ALIGNED_COLS(k));

  for (int idx = 0; idx < k; idx++) {
    OCTET_SWAP(ap[idx], bp[idx]);
  }
}

void oswapcol(octet *restrict a, uint16_t i, uint16_t j, uint16_t k,
              uint16_t l) {
  if (i == j)
    return;
  octet *ap = a;

  for (int idx = 0; idx < k; idx++, ap += ALIGNED_COLS(l)) {
    OCTET_SWAP(ap[i], ap[j]);
  }
}

void oaxpy(uint8_t *restrict a, uint8_t *restrict b, uint16_t i, uint16_t j,
           uint16_t k, uint8_t u) {
  octet *ap = a + (i * ALIGNED_COLS(k));
  octet *bp = b + (j * ALIGNED_COLS(k));

  if (u == 0)
    return;

  if (u == 1)
    return oaddrow(a, b, i, j, k);

  const octet *urow_hi = OCT_MUL_HI[u];
  const octet *urow_lo = OCT_MUL_LO[u];
  for (int idx = 0; idx < k; idx++) {
    octet b_lo = bp[idx] & 0x0f;
    octet b_hi = (bp[idx] & 0xf0) >> 4;
    ap[idx] ^= urow_hi[b_hi] ^ urow_lo[b_lo];
  }
}

void oaddrow(uint8_t *restrict a, uint8_t *restrict b, uint16_t i, uint16_t j,
             uint16_t k) {
  octet *ap = a + (i * ALIGNED_COLS(k));
  octet *bp = b + (j * ALIGNED_COLS(k));

  for (int idx = 0; idx < k; idx++) {
    ap[idx] ^= bp[idx];
  }
}

void oscal(uint8_t *restrict a, uint16_t i, uint16_t k, uint8_t u) {
  octet *ap = a + (i * ALIGNED_COLS(k));

  if (u == 0)
    return;

  const octet *urow_hi = OCT_MUL_HI[u];
  const octet *urow_lo = OCT_MUL_LO[u];
  for (int idx = 0; idx < k; idx++) {
    octet a_lo = ap[idx] & 0x0f;
    octet a_hi = (ap[idx] & 0xf0) >> 4;
    ap[idx] = urow_hi[a_hi] ^ urow_lo[a_lo];
  }
}

void ozero(uint8_t *restrict a, uint16_t i, size_t k) {
  octet *ap = a + (i * ALIGNED_COLS(k));
  for (int idx = 0; idx < k; idx++)
    ap[idx] = 0;
}

void ogemm(uint8_t *restrict a, uint8_t *restrict b, uint8_t *restrict c,
           uint16_t n, uint16_t k, uint16_t m) {
  octet *ap, *cp = c;

  for (int row = 0; row < n; row++, cp += ALIGNED_COLS(m)) {
    ap = a + (row * ALIGNED_COLS(k));

    ozero(cp, 0, m);
    for (int idx = 0; idx < k; idx++) {
      oaxpy(cp, b, 0, idx, m, ap[idx]);
    }
  }
}

void crsrow_swapcol(crsrow *a, uint16_t i, uint16_t j) {
  uint8_t tmp;
  int ix = crsrow_find(a->idxs, 0, a->nz - 1, i);
  int jx = crsrow_find(a->idxs, 0, a->nz - 1, j);

  if (ix < 0 && jx < 0) {
    return;
  } else if (ix < 0) {
    ix = a->nz;
    crsrow_resize(a, a->nz + 1);
    tmp = a->vals[jx];
    a->vals[jx] = 0;
    a->vals[ix] = tmp;
  } else if (jx < 0) {
    jx = a->nz;
    crsrow_resize(a, a->nz + 1);
    tmp = a->vals[ix];
    a->vals[ix] = 0;
    a->vals[jx] = tmp;
  } else {
    /* have both */
    tmp = a->vals[ix];
    a->vals[ix] = a->vals[jx];
    a->vals[jx] = tmp;
  }
}

void crsrow_add(crsrow *a, crsrow *b, uint16_t k) {

  octet a_dense[k];
  memset(a_dense, 0, k);

  crsrow_unpack(a, a_dense);

  for (int idx = 0; idx < b->nz; idx++) {
    a_dense[b->idxs[idx]] ^= b->vals[idx];
  }

  crsrow_pack(a, a_dense, a->nz);
}

void crsrow_axpy(crsrow *a, crsrow *b, uint8_t u, uint16_t k) {
  if (u == 0)
    return;

  octet a_dense[k];
  memset(a_dense, 0, k);

  crsrow_unpack(a, a_dense);

  const octet *urow_hi = OCT_MUL_HI[u];
  const octet *urow_lo = OCT_MUL_LO[u];
  for (int idx = 0; idx < b->nz; idx++) {
    octet b_lo = b->vals[idx] & 0x0f;
    octet b_hi = (b->vals[idx] & 0xf0) >> 4;
    a_dense[b->idxs[idx]] ^= urow_hi[b_hi] ^ urow_lo[b_lo];
  }

  crsrow_pack(a, a_dense, a->nz);
}
