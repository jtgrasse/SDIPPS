//------------------------------------------------------------------------------
// SPEX_SDIPPS/Test/sdipps_test.c: test SDIPPS column replacement
//------------------------------------------------------------------------------

/* Reads a matrix, factorizes it with SPEX LU, then applies a sequence of
 * column replacements with SDIPPS_LU_col_rep, checking the factorization
 * after every one.
 *
 * Usage:
 *
 *      sdipps_test matrix_file update_file [-c] [-g | -l] [-e] [-n N]
 *
 *          -c      prefer column pushes (default: prefer diagonal pushes)
 *          -g      eager backtracking scaling (SDIPPS_SCALING_EAGER)
 *          -l      lazy backtracking scaling (SDIPPS_SCALING_LAZY)
 *          -e      early exit
 *          -n N    check every frame against a dense REF LU only if n <= N
 *                  (default 200)
 *
 * matrix_file is in SPEX's triplet format (see ../ExampleMats).
 * update_file holds the number of updates, then for each update:
 *
 *      col nz          the column of A to replace (1-based), and the number
 *                      of nonzeros in the entering column
 *      row value       nz lines, one per nonzero (row is 1-based)
 *
 * The updates are applied in sequence, each to the factorization the previous
 * one produced.  After each update:
 *
 *  1. solve: SPEX_update_solve on the updated F solves A x = b exactly
 *     (checked in rational arithmetic by spex_demo_check_solution);
 *  2. tsolve: SPEX_update_tsolve on the updated F solves A^T x = b exactly
 *     (checked the same way, with A^T from SPEX_transpose);
 *  3. frames (n <= N): every pivot and every entry of L and U equals the
 *     REF LU of P A Q computed from scratch without pivoting, under F's own
 *     P and Q.  Pending scales are applied when reading F.
 *
 * An update that makes A singular must return SPEX_SINGULAR.  The test then
 * confirms that the would-be matrix is singular, refactorizes A (which is
 * unchanged), and goes on.  The exit status is 0 if every check passed.
 */

#include "spex_demos.h"
#include "SDIPPS.h"

//------------------------------------------------------------------------------
// read_update: read one column replacement
//------------------------------------------------------------------------------

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&vk, option) ;        \
}

static SPEX_info read_update
(
    SPEX_matrix *vk_handle,     // output: n-by-1 dynamic CSC MPZ
    int64_t *col,               // output: column to replace, 0-based
    FILE *file,
    int64_t n,
    const SPEX_options option
)
{
    SPEX_matrix vk = NULL ;
    int64_t c, nz ;
    if (fscanf (file, "%" SCNd64 " %" SCNd64, &c, &nz) != 2 ||
        c < 1 || c > n || nz < 0 || nz > n)
    {
        printf ("bad update header\n") ;
        return (SPEX_INCORRECT_INPUT) ;
    }
    (*col) = c - 1 ;

    // an empty n-by-1 dynamic CSC MPZ matrix, with room for nz entries
    SPEX_TRY (SPEX_matrix_allocate (&vk, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, 1, 0,
        false, true, option)) ;
    SPEX_TRY (SPEX_vector_realloc (vk->v [0], nz, option)) ;
    for (int64_t p = 0 ; p < nz ; p++)
    {
        int64_t row ;
        SPEX_TRY (SPEX_gmp_fscanf (file, "%" SCNd64 " %Zd", &row,
            vk->v [0]->x [p])) ;
        if (row < 1 || row > n)
        {
            printf ("bad row index %" PRId64 "\n", row) ;
            SPEX_TRY (SPEX_INCORRECT_INPUT) ;
        }
        vk->v [0]->i [p] = row - 1 ;
    }
    vk->v [0]->nz = nz ;
    (*vk_handle) = vk ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// factorize: F = SPEX LU of A (A is dynamic CSC MPZ), in updatable format
//------------------------------------------------------------------------------

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&A_csc, option) ;     \
    SPEX_symbolic_analysis_free (&S, option) ; \
}

static SPEX_info factorize
(
    SPEX_factorization *F_handle,
    const SPEX_matrix A,
    const SPEX_options option
)
{
    SPEX_matrix A_csc = NULL ;
    SPEX_symbolic_analysis S = NULL ;
    SPEX_TRY (SPEX_matrix_copy (&A_csc, SPEX_CSC, SPEX_MPZ, A, option)) ;
    SPEX_TRY (SPEX_lu_analyze (&S, A_csc, option)) ;
    SPEX_info info = SPEX_lu_factorize (F_handle, A_csc, S, option) ;
    if (info == SPEX_OK)
    {
        // L and U as columns/rows (SPEX_DYNAMIC_CSC), as SDIPPS uses them
        info = SPEX_factorization_convert (*F_handle, true, option) ;
    }
    FREE_WORKSPACE ;
    return (info) ;     // SPEX_SINGULAR is returned, not caught
}

//------------------------------------------------------------------------------
// check_solve: F solves A x = b (or A^T x = b) exactly
//------------------------------------------------------------------------------

// transpose = false: x from SPEX_update_solve, checked against A x = b.
// transpose = true:  x from SPEX_update_tsolve, checked against A^T x = b.

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&A_csc, option) ;     \
    SPEX_matrix_free (&AT, option) ;        \
    SPEX_matrix_free (&b, option) ;         \
    SPEX_matrix_free (&x, option) ;         \
}

static SPEX_info check_solve
(
    bool *ok,
    SPEX_factorization F,
    const SPEX_matrix A,        // dynamic CSC MPZ
    bool transpose,
    const SPEX_options option
)
{
    SPEX_matrix A_csc = NULL, AT = NULL, b = NULL, x = NULL ;
    int64_t n = A->n ;
    (*ok) = false ;

    // b = (1, 2, ..., n)
    SPEX_TRY (SPEX_matrix_allocate (&b, SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option)) ;
    for (int64_t i = 0 ; i < n ; i++)
    {
        SPEX_TRY (SPEX_mpz_set_si (b->x.mpz [i], i+1)) ;
    }

    // check the product with x in rational arithmetic (needs A as CSC)
    SPEX_TRY (SPEX_matrix_copy (&A_csc, SPEX_CSC, SPEX_MPZ, A, option)) ;
    if (transpose)
    {
        SPEX_TRY (SPEX_update_tsolve (&x, F, b, option)) ;
        SPEX_TRY (SPEX_transpose (&AT, A_csc, option)) ;
        (*ok) = (spex_demo_check_solution (AT, x, b, option) == SPEX_OK) ;
    }
    else
    {
        SPEX_TRY (SPEX_update_solve (&x, F, b, option)) ;
        (*ok) = (spex_demo_check_solution (A_csc, x, b, option) == SPEX_OK) ;
    }

    FREE_WORKSPACE ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// check_frames: F equals the REF LU of P A Q without pivoting
//------------------------------------------------------------------------------

// Plain GMP is used here; this is test code.  F is read as true values
// (x[p] * scale) into a dense frame matrix Fd by position: L(i,j) at
// Fd[i][j] for i >= j and U(j,i) at Fd[j][i].  B = P A Q is then reduced by
// integer-preserving Gaussian elimination without pivoting; frame k of B is
// row and column k of B after k steps, and must match Fd and rhos.

static mpz_t *dense_alloc (int64_t n)
{
    mpz_t *M = malloc (n * n * sizeof (mpz_t)) ;
    for (int64_t t = 0 ; t < n*n ; t++) mpz_init (M [t]) ;
    return (M) ;
}

static void dense_free (mpz_t *M, int64_t n)
{
    if (M == NULL) return ;
    for (int64_t t = 0 ; t < n*n ; t++) mpz_clear (M [t]) ;
    free (M) ;
}

// value = v->x[p] * v->scale; false if that is not an integer
static bool true_value (mpz_t value, const SPEX_vector v, int64_t p)
{
    if (!mpz_divisible_p (v->x [p], mpq_denref (v->scale))) return (false) ;
    mpz_divexact (value, v->x [p], mpq_denref (v->scale)) ;
    mpz_mul (value, value, mpq_numref (v->scale)) ;
    return (true) ;
}

static bool check_frames
(
    SPEX_factorization F,
    const SPEX_matrix A         // dynamic CSC MPZ, scale 1
)
{
    int64_t n = A->n ;
    int64_t *P = F->P_perm, *Pinv = F->Pinv_perm ;
    int64_t *Q = F->Q_perm, *Qinv = F->Qinv_perm ;
    mpz_t *rho = F->rhos->x.mpz ;
    mpz_t *Fd = dense_alloc (n) ;
    mpz_t *B  = dense_alloc (n) ;
    mpz_t value, prev ;
    mpz_init (value) ;
    mpz_init_set_ui (prev, 1) ;
    bool ok = true ;
    #define FD(i,j) Fd [(i)*n + (j)]
    #define BD(i,j) B  [(i)*n + (j)]

    // Fd from F, by position
    for (int64_t j = 0 ; j < n && ok ; j++)
    {
        SPEX_vector Lj = F->L->v [j], Uj = F->U->v [j] ;
        if (Lj->i [0] != P [j] || Uj->i [0] != Q [j])
        {
            printf ("  frame %" PRId64 ": pivot entry is not first\n", j) ;
            ok = false ;
        }
        for (int64_t p = 0 ; p < Lj->nz && ok ; p++)
        {
            ok = true_value (value, Lj, p) ;
            if (p == 0) ok = ok && (mpz_cmp (value, rho [j]) == 0) ;
            else mpz_set (FD (Pinv [Lj->i [p]], j), value) ;
        }
        for (int64_t p = 0 ; p < Uj->nz && ok ; p++)
        {
            ok = true_value (value, Uj, p) ;
            if (p == 0) ok = ok && (mpz_cmp (value, rho [j]) == 0) ;
            else mpz_set (FD (j, Qinv [Uj->i [p]]), value) ;
        }
        if (!ok) printf ("  frame %" PRId64 ": bad pivot or scale\n", j) ;
    }

    // B = P A Q
    for (int64_t c = 0 ; c < n && ok ; c++)
    {
        SPEX_vector Ac = A->v [c] ;
        for (int64_t p = 0 ; p < Ac->nz ; p++)
        {
            mpz_set (BD (Pinv [Ac->i [p]], Qinv [c]), Ac->x [p]) ;
        }
    }

    // integer-preserving Gaussian elimination, comparing frame k at step k
    for (int64_t k = 0 ; k < n && ok ; k++)
    {
        if (mpz_sgn (BD (k,k)) == 0 || mpz_cmp (BD (k,k), rho [k]) != 0)
        {
            printf ("  frame %" PRId64 ": wrong pivot\n", k) ;
            ok = false ;
            break ;
        }
        for (int64_t i = k+1 ; i < n && ok ; i++)
        {
            if (mpz_cmp (BD (i,k), FD (i,k)) != 0 ||
                mpz_cmp (BD (k,i), FD (k,i)) != 0)
            {
                printf ("  frame %" PRId64 ": wrong entry at position %"
                    PRId64 "\n", k, i) ;
                ok = false ;
            }
        }
        for (int64_t i = k+1 ; i < n ; i++)
        {
            for (int64_t j = k+1 ; j < n ; j++)
            {
                // B(i,j) = (B(k,k) B(i,j) - B(i,k) B(k,j)) / prev
                mpz_mul (BD (i,j), BD (i,j), BD (k,k)) ;
                mpz_submul (BD (i,j), BD (i,k), BD (k,j)) ;
                mpz_divexact (BD (i,j), BD (i,j), prev) ;
            }
        }
        mpz_set (prev, BD (k,k)) ;
    }

    #undef FD
    #undef BD
    dense_free (Fd, n) ;
    dense_free (B, n) ;
    mpz_clear (value) ;
    mpz_clear (prev) ;
    return (ok) ;
}

//------------------------------------------------------------------------------
// check_singular: true if A with column col replaced by vk is singular
//------------------------------------------------------------------------------

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    SPEX_matrix_free (&A2, option) ;        \
    SPEX_matrix_free (&vk2, option) ;       \
    SPEX_factorization_free (&F2, option) ; \
}

static SPEX_info check_singular
(
    bool *singular,
    const SPEX_matrix A,
    const SPEX_matrix vk,
    int64_t col,
    const SPEX_options option
)
{
    SPEX_matrix A2 = NULL, vk2 = NULL ;
    SPEX_factorization F2 = NULL ;
    SPEX_TRY (SPEX_matrix_copy (&A2, SPEX_DYNAMIC_CSC, SPEX_MPZ, A, option)) ;
    SPEX_TRY (SPEX_matrix_copy (&vk2, SPEX_DYNAMIC_CSC, SPEX_MPZ, vk, option));
    SPEX_TRY (SPEX_update_matrix_colrep (A2, vk2, col, option)) ;
    SPEX_info info = factorize (&F2, A2, option) ;
    (*singular) = (info == SPEX_SINGULAR) ;
    FREE_WORKSPACE ;
    return ((info == SPEX_SINGULAR) ? SPEX_OK : info) ;
}

//------------------------------------------------------------------------------
// main
//------------------------------------------------------------------------------

static const char *scaling_name (SDIPPS_scaling scaling)
{
    if (scaling == SDIPPS_SCALING_EAGER) return ("eager") ;
    if (scaling == SDIPPS_SCALING_LAZY) return ("lazy") ;
    return ("regular") ;
}

static int64_t nnz_LU (SPEX_factorization F)
{
    int64_t nz = 0 ;
    for (int64_t j = 0 ; j < F->L->n ; j++)
    {
        nz += F->L->v [j]->nz + F->U->v [j]->nz - 1 ;  // pivot counted once
    }
    return (nz) ;
}

#undef  FREE_WORKSPACE
#define FREE_WORKSPACE                      \
{                                           \
    if (mat_file != NULL) fclose (mat_file) ; \
    if (upd_file != NULL) fclose (upd_file) ; \
    SPEX_matrix_free (&A_csc, option) ;     \
    SPEX_matrix_free (&A, option) ;         \
    SPEX_matrix_free (&vk, option) ;        \
    SPEX_factorization_free (&F, option) ;  \
    SPEX_FREE (option) ;                    \
    SPEX_finalize ( ) ;                     \
}

int main (int argc, char *argv [])
{
    FILE *mat_file = NULL, *upd_file = NULL ;
    SPEX_matrix A_csc = NULL, A = NULL, vk = NULL ;
    SPEX_factorization F = NULL ;
    SPEX_options option = NULL ;

    //--------------------------------------------------------------------------
    // command line
    //--------------------------------------------------------------------------

    SDIPPS_options opt = SDIPPS_DEFAULT_OPTIONS ;
    int64_t max_check = 200 ;
    if (argc < 3)
    {
        printf ("usage: %s matrix_file update_file [-c] [-l] [-e] [-n N]\n",
            argv [0]) ;
        return (1) ;
    }
    for (int a = 3 ; a < argc ; a++)
    {
        if (strcmp (argv [a], "-c") == 0)
        {
            opt.push_choice = SDIPPS_PREFER_COLUMN ;
        }
        else if (strcmp (argv [a], "-g") == 0)
        {
            opt.scaling = SDIPPS_SCALING_EAGER ;
        }
        else if (strcmp (argv [a], "-l") == 0)
        {
            opt.scaling = SDIPPS_SCALING_LAZY ;
        }
        else if (strcmp (argv [a], "-e") == 0)
        {
            opt.early_exit = true ;
        }
        else if (strcmp (argv [a], "-n") == 0 && a+1 < argc)
        {
            max_check = atoll (argv [++a]) ;
        }
        else
        {
            printf ("unknown option %s\n", argv [a]) ;
            return (1) ;
        }
    }

    SPEX_TRY (SPEX_initialize ( )) ;
    SPEX_TRY (SPEX_create_default_options (&option)) ;

    //--------------------------------------------------------------------------
    // read A and factorize it
    //--------------------------------------------------------------------------

    mat_file = fopen (argv [1], "r") ;
    upd_file = fopen (argv [2], "r") ;
    if (mat_file == NULL || upd_file == NULL)
    {
        perror ("cannot open the input files") ;
        FREE_WORKSPACE ;
        return (1) ;
    }
    SPEX_TRY (spex_demo_tripread (&A_csc, mat_file, SPEX_MPZ, option)) ;
    SPEX_TRY (SPEX_matrix_copy (&A, SPEX_DYNAMIC_CSC, SPEX_MPZ, A_csc,
        option)) ;
    SPEX_matrix_free (&A_csc, option) ;
    int64_t n = A->n ;
    if (A->m != n)
    {
        printf ("A must be square\n") ;
        FREE_WORKSPACE ;
        return (1) ;
    }

    SPEX_TRY (factorize (&F, A, option)) ;
    bool check = (n <= max_check) ;

    printf ("SDIPPS test: %s, n = %" PRId64 ", nnz(L+U) = %" PRId64 "\n",
        argv [1], n, nnz_LU (F)) ;
    printf ("options: prefer %s pushes, %s scaling, early exit %s\n",
        (opt.push_choice == SDIPPS_PREFER_COLUMN) ? "column" : "diagonal",
        scaling_name (opt.scaling), opt.early_exit ? "on" : "off") ;
    printf ("frame check: %s\n\n", check ? "on" : "off (n too large)") ;
    printf ("  #   col     xi -> landing   diag  col    nnz(L+U)  time (s)"
        "   solve  tsolve  frames\n") ;

    //--------------------------------------------------------------------------
    // the updates, in sequence
    //--------------------------------------------------------------------------

    int64_t nupdates ;
    if (fscanf (upd_file, "%" SCNd64, &nupdates) != 1 || nupdates < 0)
    {
        printf ("bad update file\n") ;
        FREE_WORKSPACE ;
        return (1) ;
    }

    int nfail = 0 ;
    for (int64_t u = 1 ; u <= nupdates ; u++)
    {
        int64_t col ;
        SPEX_TRY (read_update (&vk, &col, upd_file, n, option)) ;

        SDIPPS_stats stats ;
        double t0 = SUITESPARSE_TIME ;
        SPEX_info status = SDIPPS_LU_col_rep (F, vk, col, &opt, &stats,
            option) ;
        double t = SUITESPARSE_TIME - t0 ;

        if (status == SPEX_SINGULAR)
        {
            // A is unchanged; confirm the exchange is singular, refactorize
            bool singular ;
            SPEX_TRY (check_singular (&singular, A, vk, col, option)) ;
            printf ("%3" PRId64 " %5" PRId64 "   singular: %s\n", u, col+1,
                singular ? "confirmed" : "WRONG, A would be nonsingular") ;
            if (!singular) nfail++ ;
            SPEX_factorization_free (&F, option) ;
            SPEX_TRY (factorize (&F, A, option)) ;
        }
        else
        {
            SPEX_TRY (status) ;
            SPEX_TRY (SPEX_update_matrix_colrep (A, vk, col, option)) ;

            // The SPEX solves apply pending scales of the factor they use
            // for forward substitution (L for solve, U for tsolve) in place,
            // so the second solve sees fewer of them.  Alternate the order
            // so that both meet the pending scales SDIPPS leaves.
            bool solve_ok, tsolve_ok, frames_ok = true ;
            bool tsolve_first = (u % 2 == 0) ;
            if (tsolve_first)
            {
                SPEX_TRY (check_solve (&tsolve_ok, F, A, true, option)) ;
            }
            SPEX_TRY (check_solve (&solve_ok, F, A, false, option)) ;
            if (!tsolve_first)
            {
                SPEX_TRY (check_solve (&tsolve_ok, F, A, true, option)) ;
            }
            if (check) frames_ok = check_frames (F, A) ;
            if (!solve_ok || !tsolve_ok || !frames_ok) nfail++ ;

            printf ("%3" PRId64 " %5" PRId64 "   %5" PRId64 " -> %5" PRId64
                "   %5" PRId64 " %4" PRId64 "  %10" PRId64 "  %8.4f"
                "    %-5s  %-6s  %s\n",
                u, col+1, stats.start+1, stats.landing+1, stats.n_diagonal,
                stats.n_column, nnz_LU (F), t, solve_ok ? "ok" : "FAIL",
                tsolve_ok ? "ok" : "FAIL",
                check ? (frames_ok ? "ok" : "FAIL") : "-") ;
        }
        SPEX_matrix_free (&vk, option) ;
    }

    printf ("\n%s (%d failure%s)\n", (nfail == 0) ? "all checks passed" :
        "TEST FAILED", nfail, (nfail == 1) ? "" : "s") ;
    FREE_WORKSPACE ;
    return ((nfail == 0) ? 0 : 1) ;
}
