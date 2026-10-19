//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/SDIPPS_Reduce: reduce a pivot-over-pivot ratio
//------------------------------------------------------------------------------

// Purpose: ratio = num / den in lowest terms.  Called once for each new
// pivot-over-pivot ratio with eager or lazy scaling: the backtracking ratio
// lambda = rho~_{k-1} / rho_k at the start of a run of pure-scaling pushes,
// and rho~_landing / rho_landing after an early exit.  Eager scaling applies
// the reduced ratio to the entries (so the reduction is its only saving);
// lazy scaling multiplies it into the vectors' pending scales.
//
// This is the place to try other reduction methods.  Whatever the method,
// the result must be in lowest terms: eager scaling divides each entry by q
// before multiplying by p, which is exact only because gcd (p, q) = 1 (x p / q
// is an integer, so q divides x) [Chen 2022, Sec. 5.4]; and lazy scaling
// multiplies it into a SPEX_vector's pending scale (an mpq_t), which GMP
// requires to be canonical.  For example, a method could first test whether
// den divides num (or num divides den) with one exact division, and only
// fall back to the GCD when it does not.
//
// Current method: GMP's GCD (mpq_canonicalize).

#include "sdipps_internal.h"

SPEX_info SDIPPS_Reduce
(
    mpq_t ratio,                // output: num / den, in lowest terms
    const mpz_t num,
    const mpz_t den             // nonzero
)
{
    SPEX_info info ;
    if (sdipps_is_zero (den))
    {
        return (SPEX_INCORRECT_INPUT) ;
    }
    SPEX_MPQ_SET_NUM (ratio, num) ;
    SPEX_MPQ_SET_DEN (ratio, den) ;
    SPEX_MPQ_CANONICALIZE (ratio) ;     // divide by gcd (num, den); den > 0
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_reduce: SDIPPS_Reduce, counted in W->stats
//------------------------------------------------------------------------------

SPEX_info sdipps_reduce
(
    sdipps_work *W,
    mpq_t ratio,
    const mpz_t num,
    const mpz_t den
)
{
    SPEX_info info ;
    int64_t bits_before = sdipps_bits (num) + sdipps_bits (den) ;
    SPEX_CHECK (SDIPPS_Reduce (ratio, num, den)) ;
    int64_t bits_saved = bits_before - sdipps_bits (SPEX_MPQ_NUM (ratio))
        - sdipps_bits (SPEX_MPQ_DEN (ratio)) ;
    W->stats.n_reduce++ ;
    if (bits_saved > 0) W->stats.n_reduced++ ;  // a common factor >= 2
    W->stats.reduce_bits_saved += bits_saved ;
    return (SPEX_OK) ;
}
