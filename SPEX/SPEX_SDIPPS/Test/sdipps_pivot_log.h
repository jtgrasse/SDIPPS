//------------------------------------------------------------------------------
// SPEX_SDIPPS/Test/sdipps_pivot_log.h: simplex pivot logs for the LP tests
//------------------------------------------------------------------------------

// Shared by sdipps_lp_test.c (updates) and sdipps_lp_refactor.c (the
// refactorization baseline).
//
// A pivot log (written by QSopt_ex) holds, with % comment lines anywhere and
// 0-based indices:
//
//      A           m N nz, then nz lines "row col value"
//      b           m nz, then nz lines "row value"
//      c           N nz, then nz lines "col value"
//      l, u        bounds (read and ignored)
//      B           the starting basis: m column indices of A
//      U           one per basis change: "leaving entering"
//      x, y, d     the final exact solution (only y is kept); z (ignored)
//
// Values are exact rationals "p/q" or integers.
//
// A is made integral column by column: column j is multiplied by the LCM
// s_j of its denominators.  Every basis matrix built here is therefore B D,
// D = diag (s_j) over the basic columns, and
//
//      B x = b        becomes  (B D) x' = b,          x = D x'
//      B^T y = c_B    becomes  (B D)^T y = D c_B      (the same y)
//
// The basis is kept as basis[slot] (the column of A in column `slot` of the
// basis matrix) and slot_of[col] (-1 if col is not basic).  An entering
// column takes the slot of the leaving one.
//
// Read the log BEFORE SPEX_initialize: the reader calls GMP directly, and
// after SPEX_initialize GMP lists each block it allocates until the next
// SPEX GMP wrapper call, so many direct calls in a row become very slow.

#ifndef SDIPPS_PIVOT_LOG_H
#define SDIPPS_PIVOT_LOG_H

#include "spex_demos.h"

typedef struct
{
    int64_t m, N ;          // A is m-by-N
    int64_t nz ;            // number of entries of A
    int64_t *Ap ;           // A by column, made integral column by column:
    int64_t *Ai ;           // column j is Ai, Ax [Ap [j] ... Ap [j+1]-1]
    mpz_t *Ax ;
    mpz_t *colscale ;       // column j of A was multiplied by colscale [j]
    int64_t *colbits ;      // total bits of the (integral) entries of column j
    mpq_t *b ;              // right-hand side, dense, size m
    mpq_t *c ;              // objective, dense, size N
    mpq_t *y ;              // QSopt_ex's final dual solution, size m (or NULL)
    int64_t *basis0 ;       // starting basis, m column indices
    int64_t nupdates ;
    int64_t *leaving ;      // update k: column leaving[k] leaves the basis,
    int64_t *entering ;     // column entering[k] enters it
}
pivot_log ;

// read a pivot log; false (with a message) if it cannot be read
bool lp_read (const char *filename, pivot_log *L) ;

void lp_free (pivot_log *L) ;

//------------------------------------------------------------------------------
// the basis
//------------------------------------------------------------------------------

// basis = the starting basis; allocates basis (size m) and slot_of (size N)
bool lp_basis_start (const pivot_log *L, int64_t **basis, int64_t **slot_of) ;

// the slot of update k's leaving column, or -1 (with a message) if update k
// does not fit the current basis
int64_t lp_update_slot (const pivot_log *L, const int64_t *slot_of, int64_t k);

// apply update k to the basis (it must fit: see lp_update_slot)
void lp_basis_apply (const pivot_log *L, int64_t *basis, int64_t *slot_of,
    int64_t k) ;

//------------------------------------------------------------------------------
// SPEX matrices
//------------------------------------------------------------------------------

// B D = A (:, basis), as SPEX_CSC SPEX_MPZ
SPEX_info lp_basis_matrix (SPEX_matrix *B_handle, const pivot_log *L,
    const int64_t *basis, const SPEX_options option) ;

// column j of A (integral), as an m-by-1 SPEX_DYNAMIC_CSC SPEX_MPZ matrix
SPEX_info lp_column (SPEX_matrix *vk_handle, const pivot_log *L, int64_t j,
    const SPEX_options option) ;

// b, as an m-by-1 SPEX_DENSE SPEX_MPZ matrix (made integral, see b->scale)
SPEX_info lp_rhs_b (SPEX_matrix *b_handle, const pivot_log *L,
    const SPEX_options option) ;

// D c_B, as an m-by-1 SPEX_DENSE SPEX_MPZ matrix (made integral)
SPEX_info lp_rhs_cB (SPEX_matrix *c_handle, const pivot_log *L,
    const int64_t *basis, const SPEX_options option) ;

// SPEX LU of B D (analyze + factorize), in updatable form; the time of
// analyze + factorize is returned in *t_factor
SPEX_info lp_factorize (SPEX_factorization *F_handle, double *t_factor,
    const pivot_log *L, const int64_t *basis, const SPEX_options option) ;

//------------------------------------------------------------------------------
// solves, checks and metrics
//------------------------------------------------------------------------------

// x = solution of (B D) x = b with SPEX_update_solve, and, if do_tsolve,
// y = solution of (B D)^T y = D c_B with SPEX_update_tsolve (otherwise y is
// left NULL).  The solves are timed; D c_B is built for the current basis
// (not timed).
SPEX_info lp_solves (SPEX_matrix *x, SPEX_matrix *y, double *t_solve,
    double *t_tsolve, bool do_tsolve, SPEX_factorization F,
    const SPEX_matrix b, const pivot_log *L, const int64_t *basis,
    const SPEX_options option) ;

// Check x (and y, unless it is NULL) exactly against B D built from
// scratch.  If final_basis and the log has QSopt_ex's y, also check that y
// equals it.
SPEX_info lp_check (bool *ok, const SPEX_matrix x, const SPEX_matrix y,
    const SPEX_matrix b, const pivot_log *L, const int64_t *basis,
    bool final_basis, const SPEX_options option) ;

// Size of a factorization.  Entries are counted and measured as stored: with
// a pending scale, the stored x[p], not x[p] * scale.
typedef struct
{
    int64_t nnz_L, nnz_U ;      // nonzero values, pivots included (SPEX LU
                                // also stores entries that cancelled to
                                // zero; those are not counted)
    int64_t bits ;              // total bits of every integer F stores: the
                                // entries of L and U, the pending scales
                                // (numerator and denominator), the pivots
    int64_t max_bits_L ;        // largest off-diagonal entry of L, in bits
    int64_t max_bits_U ;        // largest off-diagonal entry of U, in bits
    int64_t max_bits_rho ;      // largest pivot (rhos), in bits
    int64_t max_bits_scale ;    // largest pending scale p/q of a vector:
                                // bits (p) + bits (q); 0 if there is none
}
lp_metrics ;

void lp_factor_metrics (SPEX_factorization F, lp_metrics *M) ;

// number of nonzero values in v (explicit zeros are not counted)
int64_t lp_vector_nnz (const SPEX_vector v) ;

// total bits of the entries of B D
int64_t lp_basis_bits (const pivot_log *L, const int64_t *basis) ;

#endif
