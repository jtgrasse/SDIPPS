//------------------------------------------------------------------------------
// SPEX_SDIPPS/Source/sdipps_vector: SPEX_vector utilities for SDIPPS
//------------------------------------------------------------------------------

// A column of L or row of U is a SPEX_vector: x[0] is the pivot entry, the
// rest are the off-diagonal entries in any order, keyed by original index.
// Its true entries are x[p] * scale (SPEX's pending scale).

#include "sdipps_internal.h"

//------------------------------------------------------------------------------
// sdipps_apply_scale: make x[p] the true values and set the scale to 1
//------------------------------------------------------------------------------

// As in SPEX, x[p] / den(scale) is an integer, so dividing first keeps the
// numbers small.  The pivot is taken from the given true pivot.

SPEX_info sdipps_apply_scale
(
    SPEX_vector v,
    const mpz_t pivot           // the true pivot of v's frame
)
{
    SPEX_info info ;
    int r ;
    SPEX_MPQ_CMP_UI (&r, v->scale, 1, 1) ;
    if (r == 0)
    {
        return (SPEX_OK) ;      // no pending scale
    }

    for (int64_t p = 1 ; p < v->nz ; p++)
    {
        SPEX_MPZ_DIVEXACT (v->x [p], v->x [p], SPEX_MPQ_DEN (v->scale)) ;
        SPEX_MPZ_MUL (v->x [p], v->x [p], SPEX_MPQ_NUM (v->scale)) ;
    }
    SPEX_MPZ_SET (v->x [0], pivot) ;
    SPEX_MPQ_SET_UI (v->scale, 1, 1) ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_set_pivot: x[0] = pivot / scale, so that x[0] * scale == pivot
//------------------------------------------------------------------------------

SPEX_info sdipps_set_pivot
(
    SPEX_vector v,
    const mpz_t pivot
)
{
    SPEX_info info ;
    SPEX_MPZ_MUL (v->x [0], pivot, SPEX_MPQ_DEN (v->scale)) ;
    SPEX_MPZ_DIVEXACT (v->x [0], v->x [0], SPEX_MPQ_NUM (v->scale)) ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_reserve: make room for `need` entries
//------------------------------------------------------------------------------

SPEX_info sdipps_reserve
(
    SPEX_vector v,
    int64_t need,
    const SPEX_options option
)
{
    SPEX_info info ;
    if (v->nzmax < need)
    {
        SPEX_CHECK (SPEX_vector_realloc (v, need, option)) ;
    }
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_start_vector: v <- (pivot_index : pivot), scale 1
//------------------------------------------------------------------------------

SPEX_info sdipps_start_vector
(
    SPEX_vector v,
    int64_t pivot_index,
    const mpz_t pivot,
    int64_t need,               // room for this many entries, pivot included
    const SPEX_options option
)
{
    SPEX_info info ;
    SPEX_CHECK (sdipps_reserve (v, need, option)) ;
    v->i [0] = pivot_index ;
    SPEX_MPZ_SET (v->x [0], pivot) ;
    v->nz = 1 ;
    SPEX_MPQ_SET_UI (v->scale, 1, 1) ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_append: append (index, value) if value != 0
//------------------------------------------------------------------------------

// Zeros are never stored.  v must have scale 1 and room for one more entry.

SPEX_info sdipps_append
(
    SPEX_vector v,
    int64_t index,
    const mpz_t value
)
{
    SPEX_info info ;
    if (sdipps_is_zero (value))
    {
        return (SPEX_OK) ;
    }
    ASSERT (v->nz < v->nzmax) ;
    v->i [v->nz] = index ;
    SPEX_MPZ_SET (v->x [v->nz], value) ;
    v->nz++ ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// sdipps_remove_entry: remove the off-diagonal entry at original index `index`
//------------------------------------------------------------------------------

// The last entry takes its place; the pivot stays first.

SPEX_info sdipps_remove_entry
(
    SPEX_vector v,
    int64_t index
)
{
    SPEX_info info ;
    for (int64_t p = 1 ; p < v->nz ; p++)
    {
        if (v->i [p] == index)
        {
            int64_t last = v->nz - 1 ;
            v->i [p] = v->i [last] ;
            SPEX_MPZ_SWAP (v->x [p], v->x [last]) ;
            v->nz-- ;
            break ;
        }
    }
    return (SPEX_OK) ;
}
