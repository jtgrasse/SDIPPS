//------------------------------------------------------------------------------
// SPEX_SDIPPS/Test/sdipps_pivot_log.c: simplex pivot logs for the LP tests
//------------------------------------------------------------------------------

// See sdipps_pivot_log.h.  Plain GMP is used while reading (before
// SPEX_initialize); SPEX's GMP wrappers are used afterwards.

#include <ctype.h>
#include "sdipps_pivot_log.h"

//------------------------------------------------------------------------------
// lp_free
//------------------------------------------------------------------------------

static void free_mpq_array (mpq_t *x, int64_t n)
{
    if (x == NULL) return ;
    for (int64_t i = 0 ; i < n ; i++) mpq_clear (x [i]) ;
    free (x) ;
}

void lp_free (pivot_log *L)
{
    if (L->Ax != NULL)
    {
        for (int64_t p = 0 ; p < L->nz ; p++) mpz_clear (L->Ax [p]) ;
    }
    if (L->colscale != NULL)
    {
        for (int64_t j = 0 ; j < L->N ; j++) mpz_clear (L->colscale [j]) ;
    }
    free_mpq_array (L->b, L->m) ;
    free_mpq_array (L->c, L->N) ;
    free_mpq_array (L->y, L->m) ;
    free (L->Ap) ; free (L->Ai) ; free (L->Ax) ; free (L->colscale) ;
    free (L->colbits) ; free (L->basis0) ; free (L->leaving) ;
    free (L->entering) ;
    memset (L, 0, sizeof (pivot_log)) ;
}

//------------------------------------------------------------------------------
// reading tokens
//------------------------------------------------------------------------------

// A token buffer that grows as needed (rationals can have any length).
typedef struct
{
    char *s ;
    size_t size ;
}
token_buf ;

// Read the next white-space-separated token, skipping comments (from % to
// the end of the line).  Returns false at the end of the file.
static bool next_token (FILE *f, token_buf *t)
{
    int ch ;
    for (;;)
    {
        ch = fgetc (f) ;
        if (ch == EOF) return (false) ;
        if (ch == '%')
        {
            while (ch != EOF && ch != '\n') ch = fgetc (f) ;
            if (ch == EOF) return (false) ;
            continue ;
        }
        if (!isspace (ch)) break ;
    }
    size_t k = 0 ;
    while (ch != EOF && !isspace (ch))
    {
        if (k + 1 >= t->size)
        {
            size_t size = 2 * t->size ;
            char *s = realloc (t->s, size) ;
            if (s == NULL) return (false) ;
            t->s = s ;
            t->size = size ;
        }
        t->s [k++] = (char) ch ;
        ch = fgetc (f) ;
    }
    t->s [k] = '\0' ;
    return (true) ;
}

static bool read_int (FILE *f, token_buf *t, int64_t *x)
{
    if (!next_token (f, t)) return (false) ;
    char *end ;
    (*x) = strtoll (t->s, &end, 10) ;
    return (*end == '\0') ;
}

// read "p/q" or "p" into x, in canonical form
static bool read_mpq (FILE *f, token_buf *t, mpq_t x)
{
    if (!next_token (f, t)) return (false) ;
    if (mpq_set_str (x, t->s, 10) != 0) return (false) ;
    if (mpz_sgn (mpq_denref (x)) == 0) return (false) ;
    mpq_canonicalize (x) ;
    return (true) ;
}

//------------------------------------------------------------------------------
// the sections
//------------------------------------------------------------------------------

// A: read the triplets, sort them by column, make each column integral
static bool read_A (FILE *f, token_buf *t, pivot_log *L)
{
    int64_t m, N, nz ;
    if (!read_int (f, t, &m) || !read_int (f, t, &N) || !read_int (f, t, &nz)
        || m <= 0 || N <= 0 || nz < 0 || L->Ap != NULL) return (false) ;
    L->m = m ;
    L->N = N ;

    int64_t *ti = malloc (nz * sizeof (int64_t)) ;
    int64_t *tj = malloc (nz * sizeof (int64_t)) ;
    mpq_t *tx = malloc (nz * sizeof (mpq_t)) ;
    mpq_t *cx = malloc (nz * sizeof (mpq_t)) ;
    int64_t *next = calloc (N+1, sizeof (int64_t)) ;
    L->Ap = calloc (N+1, sizeof (int64_t)) ;
    L->Ai = malloc (nz * sizeof (int64_t)) ;
    L->Ax = malloc (nz * sizeof (mpz_t)) ;
    L->colscale = malloc (N * sizeof (mpz_t)) ;
    L->colbits = calloc (N, sizeof (int64_t)) ;
    bool ok = (ti && tj && tx && cx && next && L->Ap && L->Ai && L->Ax &&
        L->colscale && L->colbits) ;
    if (!ok)
    {
        free (ti) ; free (tj) ; free (tx) ; free (cx) ; free (next) ;
        free (L->Ax) ; L->Ax = NULL ;
        free (L->colscale) ; L->colscale = NULL ;
        return (false) ;
    }
    for (int64_t p = 0 ; p < nz ; p++) mpq_init (tx [p]) ;
    for (int64_t p = 0 ; p < nz ; p++) mpz_init (L->Ax [p]) ;
    for (int64_t j = 0 ; j < N ; j++) mpz_init (L->colscale [j]) ;
    L->nz = nz ;

    // the triplets
    for (int64_t p = 0 ; p < nz && ok ; p++)
    {
        ok = read_int (f, t, &ti [p]) && read_int (f, t, &tj [p]) &&
             read_mpq (f, t, tx [p]) &&
             ti [p] >= 0 && ti [p] < m && tj [p] >= 0 && tj [p] < N ;
    }

    if (ok)
    {
        // sort by column (counting sort); the mpq_t values are moved
        for (int64_t p = 0 ; p < nz ; p++) L->Ap [tj [p] + 1]++ ;
        for (int64_t j = 0 ; j < N ; j++) L->Ap [j+1] += L->Ap [j] ;
        for (int64_t j = 0 ; j < N ; j++) next [j] = L->Ap [j] ;
        for (int64_t p = 0 ; p < nz ; p++)
        {
            int64_t q = next [tj [p]]++ ;
            L->Ai [q] = ti [p] ;
            memcpy (cx [q], tx [p], sizeof (mpq_t)) ;
        }
        memcpy (tx, cx, nz * sizeof (mpq_t)) ;

        // column j times s_j = lcm of its denominators
        for (int64_t j = 0 ; j < N ; j++)
        {
            mpz_set_ui (L->colscale [j], 1) ;
            for (int64_t p = L->Ap [j] ; p < L->Ap [j+1] ; p++)
            {
                mpz_lcm (L->colscale [j], L->colscale [j],
                    mpq_denref (tx [p])) ;
            }
            for (int64_t p = L->Ap [j] ; p < L->Ap [j+1] ; p++)
            {
                // A(i,j) s_j = num * (s_j / den)
                mpz_divexact (L->Ax [p], L->colscale [j], mpq_denref (tx [p]));
                mpz_mul (L->Ax [p], L->Ax [p], mpq_numref (tx [p])) ;
                if (mpz_sgn (L->Ax [p]) != 0)
                {
                    L->colbits [j] += mpz_sizeinbase (L->Ax [p], 2) ;
                }
            }
        }
    }

    for (int64_t p = 0 ; p < nz ; p++) mpq_clear (tx [p]) ;
    free (ti) ; free (tj) ; free (tx) ; free (cx) ; free (next) ;
    return (ok) ;
}

// A sparse vector section: "n nz", then nz lines "index value".  With
// x_handle != NULL it is read into a dense rational vector of size n, which
// must equal n_expected; with x_handle == NULL it is skipped without parsing
// the values (bounds may be -inf or inf).
static bool read_vector (FILE *f, token_buf *t, mpq_t **x_handle,
    int64_t n_expected)
{
    int64_t n, nz ;
    if (!read_int (f, t, &n) || !read_int (f, t, &nz) || nz < 0) return (false);
    if (x_handle == NULL)
    {
        for (int64_t p = 0 ; p < 2*nz ; p++)
        {
            if (!next_token (f, t)) return (false) ;
        }
        return (true) ;
    }
    if (n != n_expected || *x_handle != NULL) return (false) ;
    mpq_t *x = malloc (n * sizeof (mpq_t)) ;
    if (x == NULL) return (false) ;
    for (int64_t i = 0 ; i < n ; i++) mpq_init (x [i]) ;
    (*x_handle) = x ;
    bool ok = true ;
    for (int64_t p = 0 ; p < nz && ok ; p++)
    {
        int64_t i ;
        ok = read_int (f, t, &i) && i >= 0 && i < n && read_mpq (f, t, x [i]) ;
    }
    return (ok) ;
}

bool lp_read (const char *filename, pivot_log *L)
{
    memset (L, 0, sizeof (pivot_log)) ;
    FILE *f = fopen (filename, "r") ;
    if (f == NULL)
    {
        perror (filename) ;
        return (false) ;
    }
    token_buf t ;
    t.size = 256 ;
    t.s = malloc (t.size) ;
    int64_t cap = 0 ;
    bool ok = (t.s != NULL) ;

    while (ok && next_token (f, &t))
    {
        const char *s = t.s ;
        bool after_A = (L->Ap != NULL) ;
        if (strcmp (s, "A") == 0)
        {
            ok = read_A (f, &t, L) ;
        }
        else if (strcmp (s, "b") == 0)
        {
            ok = after_A && read_vector (f, &t, &L->b, L->m) ;
        }
        else if (strcmp (s, "c") == 0)
        {
            ok = after_A && read_vector (f, &t, &L->c, L->N) ;
        }
        else if (strcmp (s, "y") == 0)
        {
            ok = after_A && read_vector (f, &t, &L->y, L->m) ;
        }
        else if (strcmp (s, "l") == 0 || strcmp (s, "u") == 0 ||
                 strcmp (s, "x") == 0 || strcmp (s, "d") == 0)
        {
            ok = read_vector (f, &t, NULL, 0) ;
        }
        else if (strcmp (s, "z") == 0)
        {
            ok = next_token (f, &t) ;               // the objective value
        }
        else if (strcmp (s, "B") == 0)
        {
            ok = after_A && L->basis0 == NULL ;
            if (ok) L->basis0 = malloc (L->m * sizeof (int64_t)) ;
            ok = ok && (L->basis0 != NULL) ;
            for (int64_t k = 0 ; k < L->m && ok ; k++)
            {
                ok = read_int (f, &t, &L->basis0 [k]) &&
                     L->basis0 [k] >= 0 && L->basis0 [k] < L->N ;
            }
        }
        else if (strcmp (s, "U") == 0)
        {
            if (L->nupdates == cap)
            {
                cap = (cap == 0) ? 1024 : 2*cap ;
                int64_t *lv = realloc (L->leaving, cap * sizeof (int64_t)) ;
                if (lv != NULL) L->leaving = lv ;
                int64_t *en = realloc (L->entering, cap * sizeof (int64_t)) ;
                if (en != NULL) L->entering = en ;
                ok = (lv != NULL && en != NULL) ;
            }
            int64_t k = L->nupdates ;
            ok = ok && read_int (f, &t, &L->leaving [k]) &&
                 read_int (f, &t, &L->entering [k]) ;
            L->nupdates++ ;
        }
        else
        {
            printf ("unexpected token \"%s\" in %s\n", s, filename) ;
            ok = false ;
        }
    }

    ok = ok && L->Ap && L->b && L->c && L->basis0 ;
    if (!ok) printf ("could not read the pivot log %s\n", filename) ;
    free (t.s) ;
    fclose (f) ;
    return (ok) ;
}

//------------------------------------------------------------------------------
// the basis
//------------------------------------------------------------------------------

bool lp_basis_start (const pivot_log *L, int64_t **basis, int64_t **slot_of)
{
    (*basis) = malloc (L->m * sizeof (int64_t)) ;
    (*slot_of) = malloc (L->N * sizeof (int64_t)) ;
    if (*basis == NULL || *slot_of == NULL) return (false) ;
    for (int64_t j = 0 ; j < L->N ; j++) (*slot_of) [j] = -1 ;
    for (int64_t k = 0 ; k < L->m ; k++)
    {
        (*basis) [k] = L->basis0 [k] ;
        (*slot_of) [L->basis0 [k]] = k ;
    }
    return (true) ;
}

int64_t lp_update_slot (const pivot_log *L, const int64_t *slot_of, int64_t k)
{
    int64_t lv = L->leaving [k], en = L->entering [k] ;
    bool ok = (lv >= 0 && lv < L->N && en >= 0 && en < L->N &&
               slot_of [lv] >= 0 && slot_of [en] < 0) ;
    if (!ok)
    {
        printf ("update %" PRId64 " (%" PRId64 " leaves, %" PRId64
            " enters) does not match the basis\n", k+1, lv, en) ;
        return (-1) ;
    }
    return (slot_of [lv]) ;
}

void lp_basis_apply (const pivot_log *L, int64_t *basis, int64_t *slot_of,
    int64_t k)
{
    int64_t lv = L->leaving [k], en = L->entering [k] ;
    int64_t slot = slot_of [lv] ;
    basis [slot] = en ;
    slot_of [en] = slot ;
    slot_of [lv] = -1 ;
}

//------------------------------------------------------------------------------
// SPEX matrices
//------------------------------------------------------------------------------

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&B, option) ;         \
}

SPEX_info lp_basis_matrix (SPEX_matrix *B_handle, const pivot_log *L,
    const int64_t *basis, const SPEX_options option)
{
    SPEX_matrix B = NULL ;
    int64_t m = L->m, nz = 0 ;
    for (int64_t k = 0 ; k < m ; k++)
    {
        nz += L->Ap [basis [k] + 1] - L->Ap [basis [k]] ;
    }
    SPEX_TRY (SPEX_matrix_allocate (&B, SPEX_CSC, SPEX_MPZ, m, m, nz, false,
        true, option)) ;
    int64_t q = 0 ;
    for (int64_t k = 0 ; k < m ; k++)
    {
        B->p [k] = q ;
        int64_t j = basis [k] ;
        for (int64_t p = L->Ap [j] ; p < L->Ap [j+1] ; p++)
        {
            B->i [q] = L->Ai [p] ;
            SPEX_TRY (SPEX_mpz_set (B->x.mpz [q], L->Ax [p])) ;
            q++ ;
        }
    }
    B->p [m] = q ;
    (*B_handle) = B ;
    return (SPEX_OK) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&vk, option) ;        \
}

SPEX_info lp_column (SPEX_matrix *vk_handle, const pivot_log *L, int64_t j,
    const SPEX_options option)
{
    SPEX_matrix vk = NULL ;
    int64_t nz = L->Ap [j+1] - L->Ap [j] ;
    SPEX_TRY (SPEX_matrix_allocate (&vk, SPEX_DYNAMIC_CSC, SPEX_MPZ, L->m, 1,
        0, false, true, option)) ;
    SPEX_TRY (SPEX_vector_realloc (vk->v [0], nz, option)) ;
    for (int64_t p = 0 ; p < nz ; p++)
    {
        vk->v [0]->i [p] = L->Ai [L->Ap [j] + p] ;
        SPEX_TRY (SPEX_mpz_set (vk->v [0]->x [p], L->Ax [L->Ap [j] + p])) ;
    }
    vk->v [0]->nz = nz ;
    (*vk_handle) = vk ;
    return (SPEX_OK) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&v, option) ;         \
}

// v = a dense rational vector, as an m-by-1 dense MPZ matrix made integral
// by SPEX_matrix_copy (which records the factor in the result's scale)
static SPEX_info dense_rhs (SPEX_matrix *rhs_handle, mpq_t *x, int64_t m,
    const SPEX_options option)
{
    SPEX_matrix v = NULL ;
    SPEX_TRY (SPEX_matrix_allocate (&v, SPEX_DENSE, SPEX_MPQ, m, 1, m, false,
        true, option)) ;
    for (int64_t i = 0 ; i < m ; i++)
    {
        SPEX_TRY (SPEX_mpq_set (v->x.mpq [i], x [i])) ;
    }
    SPEX_TRY (SPEX_matrix_copy (rhs_handle, SPEX_DENSE, SPEX_MPZ, v, option)) ;
    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

SPEX_info lp_rhs_b (SPEX_matrix *b_handle, const pivot_log *L,
    const SPEX_options option)
{
    return (dense_rhs (b_handle, L->b, L->m, option)) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&v, option) ;         \
}

SPEX_info lp_rhs_cB (SPEX_matrix *c_handle, const pivot_log *L,
    const int64_t *basis, const SPEX_options option)
{
    // D c_B: slot k holds c [basis [k]] * colscale [basis [k]]
    SPEX_matrix v = NULL ;
    int64_t m = L->m ;
    SPEX_TRY (SPEX_matrix_allocate (&v, SPEX_DENSE, SPEX_MPQ, m, 1, m, false,
        true, option)) ;
    for (int64_t k = 0 ; k < m ; k++)
    {
        int64_t j = basis [k] ;
        SPEX_TRY (SPEX_mpq_set_z (v->x.mpq [k], L->colscale [j])) ;
        SPEX_TRY (SPEX_mpq_mul (v->x.mpq [k], v->x.mpq [k], L->c [j])) ;
    }
    SPEX_TRY (SPEX_matrix_copy (c_handle, SPEX_DENSE, SPEX_MPZ, v, option)) ;
    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                          \
{                                               \
    SPEX_matrix_free (&B, option) ;             \
    SPEX_symbolic_analysis_free (&S, option) ;  \
}

SPEX_info lp_factorize (SPEX_factorization *F_handle, double *t_factor,
    const pivot_log *L, const int64_t *basis, const SPEX_options option)
{
    SPEX_matrix B = NULL ;
    SPEX_symbolic_analysis S = NULL ;
    SPEX_TRY (lp_basis_matrix (&B, L, basis, option)) ;
    double t0 = SUITESPARSE_TIME ;
    SPEX_TRY (SPEX_lu_analyze (&S, B, option)) ;
    SPEX_TRY (SPEX_lu_factorize (F_handle, B, S, option)) ;
    (*t_factor) = SUITESPARSE_TIME - t0 ;
    // the updatable form (columns of L, rows of U), as SDIPPS uses
    SPEX_TRY (SPEX_factorization_convert (*F_handle, true, option)) ;
    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// solves and checks
//------------------------------------------------------------------------------

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&cB, option) ;        \
}

SPEX_info lp_solves (SPEX_matrix *x, SPEX_matrix *y, double *t_solve,
    double *t_tsolve, bool do_tsolve, SPEX_factorization F,
    const SPEX_matrix b, const pivot_log *L, const int64_t *basis,
    const SPEX_options option)
{
    SPEX_matrix cB = NULL ;

    double t0 = SUITESPARSE_TIME ;
    SPEX_TRY (SPEX_update_solve (x, F, b, option)) ;
    (*t_solve) = SUITESPARSE_TIME - t0 ;

    (*t_tsolve) = 0 ;
    if (do_tsolve)
    {
        SPEX_TRY (lp_rhs_cB (&cB, L, basis, option)) ;
        t0 = SUITESPARSE_TIME ;
        SPEX_TRY (SPEX_update_tsolve (y, F, cB, option)) ;
        (*t_tsolve) = SUITESPARSE_TIME - t0 ;
    }

    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&B, option) ;         \
    SPEX_matrix_free (&BT, option) ;        \
    SPEX_matrix_free (&cB, option) ;        \
}

SPEX_info lp_check (bool *ok, const SPEX_matrix x, const SPEX_matrix y,
    const SPEX_matrix b, const pivot_log *L, const int64_t *basis,
    bool final_basis, const SPEX_options option)
{
    SPEX_matrix B = NULL, BT = NULL, cB = NULL ;
    SPEX_TRY (lp_basis_matrix (&B, L, basis, option)) ;

    // (B D) x = b, exactly
    bool ok_x = (spex_demo_check_solution (B, x, b, option) == SPEX_OK) ;
    printf ("    solve %s", ok_x ? "exact" : "WRONG") ;

    // (B D)^T y = D c_B, exactly
    bool ok_y = true ;
    if (y != NULL)
    {
        SPEX_TRY (SPEX_transpose (&BT, B, option)) ;
        SPEX_TRY (lp_rhs_cB (&cB, L, basis, option)) ;
        ok_y = (spex_demo_check_solution (BT, y, cB, option) == SPEX_OK) ;
        printf (", transposed solve %s", ok_y ? "exact" : "WRONG") ;
    }

    // the final basis: y must be QSopt_ex's dual solution
    bool ok_qs = true ;
    if (y != NULL && final_basis && L->y != NULL)
    {
        for (int64_t i = 0 ; i < L->m && ok_qs ; i++)
        {
            ok_qs = (mpq_equal (y->x.mpq [i], L->y [i]) != 0) ;
        }
        printf (", y %s QSopt_ex's dual solution",
            ok_qs ? "equals" : "DIFFERS FROM") ;
    }
    printf ("\n") ;

    (*ok) = ok_x && ok_y && ok_qs ;
    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// metrics
//------------------------------------------------------------------------------

// number of bits of |x| (0 for x = 0); mpz_sizeinbase does not allocate
static int64_t bits_of (const mpz_t x)
{
    return ((mpz_sgn (x) == 0) ? 0 : (int64_t) mpz_sizeinbase (x, 2)) ;
}

int64_t lp_vector_nnz (const SPEX_vector v)
{
    int64_t count = 0 ;
    for (int64_t p = 0 ; p < v->nz ; p++)
    {
        if (mpz_sgn (v->x [p]) != 0) count++ ;
    }
    return (count) ;
}

// add v's entries and pending scale to M: its total bits, and the largest
// off-diagonal entry in *max_entry; the pivot entry is measured with rhos
static void vector_metrics (const SPEX_vector v, lp_metrics *M,
    int64_t *max_entry)
{
    for (int64_t p = 0 ; p < v->nz ; p++)
    {
        int64_t bits = bits_of (v->x [p]) ;
        M->bits += bits ;
        if (p > 0 && bits > *max_entry) *max_entry = bits ;
    }
    if (mpq_cmp_ui (v->scale, 1, 1) != 0)
    {
        int64_t bits = bits_of (mpq_numref (v->scale)) +
                       bits_of (mpq_denref (v->scale)) ;
        M->bits += bits ;
        if (bits > M->max_bits_scale) M->max_bits_scale = bits ;
    }
}

void lp_factor_metrics (SPEX_factorization F, lp_metrics *M)
{
    memset (M, 0, sizeof (lp_metrics)) ;
    for (int64_t j = 0 ; j < F->L->n ; j++)
    {
        M->nnz_L += lp_vector_nnz (F->L->v [j]) ;
        M->nnz_U += lp_vector_nnz (F->U->v [j]) ;
        vector_metrics (F->L->v [j], M, &M->max_bits_L) ;
        vector_metrics (F->U->v [j], M, &M->max_bits_U) ;
        int64_t bits = bits_of (F->rhos->x.mpz [j]) ;
        M->bits += bits ;
        if (bits > M->max_bits_rho) M->max_bits_rho = bits ;
    }
}

int64_t lp_basis_bits (const pivot_log *L, const int64_t *basis)
{
    int64_t bits = 0 ;
    for (int64_t k = 0 ; k < L->m ; k++) bits += L->colbits [basis [k]] ;
    return (bits) ;
}
