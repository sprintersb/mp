#ifndef MP_CHEBY_H
#define MP_CHEBY_H

/* A Short Reminder on Properties of Tn, Chebyshev Polynomials of the 1st Kind:

    Degree: n

    T_n(x) = 2x T_{n-1}(x) - T_{n-2}(x)
    
           = cos (n acos x),                x in [-1, 1]
           = cosh (n arcosh x),             x > 1
           = (-1)^n cosh (n arcosh (-1)),   x < -1

           = ((x + w)^n + (x - w)^n) / 2,   w = sqrt(x^2-1)
    
    T_n (cos x) = cos (n x)

    T_n is the polynomial of degree n with maximal leading coefficient
    such that  | T_n(x) | <= 1  for all  x  in  [-1,1]

    n zeros in (-1,1) are of order 1 at: cos ((2k+1)/2n * pi),  k = 0 ... n-1.

    n + 1 extrema in [-1,1] at:  cos (k/n * pi), k = 0 ... n
    with  | T_n(x_k) | = 1  and  T_n(x_k)  have alternating signs.
    
    DGL:  (1-x^2) y" - x y' + n^2 y = 0

    Orthogonal:  W.r.t. integral over [-1,1] with kernel 1 / sqrt(1-x^2).
*/

#include "mp-ratio.h"
#include "mp-poly.h"
#include "mp-matrix.h"

namespace Cheby
{
    // Chebyshev polynomial of th 1st kind of degree n in Z[x].
    const Poly<Ratio>& Tn (int);

    // Clear the Tn cache.
    // Invalidates all references returned by the function above.
    void Tn_clear();

    // Transform from basis { T_0, T_1, ... T_n }  to  { 1, x, x^2, ... x^n}.
    // It's a (1+n) x (1+n) matrix that:
    //  * It's upper triangular, all diagonal elements != 0.
    //  * All entries are integers.
    FMatrix<Float> Tn_matrix (int);

    // Change of basis T_n <-> x^n.  DIM must be at least 1 + degree
    // of the polynomial or -1 (defaults to 1 + degree).
    Poly<Float> basis_Tn_to_xn (const Poly<Float>&, int dim = -1);
    Poly<Float> basis_xn_to_Tn (const Poly<Float>&, int dim = -1);

    FVector<Float> basis_Tn_to_xn (const FVector<Float>&);
    FVector<Float> basis_xn_to_Tn (const FVector<Float>&);

    // The k-th zero of T_n lies in (-1, 1) and are in increasing order:
    // x_k < x_j  <=>  k < j.
    template<class F> F Tn_zero (int n , int k); // 0 <= k < n.

    // The k-th maximum (by absolute value) lies in [-1,1].
    // The maxima at -1 and 1 are due to the interval boundaries;
    // the other n-1 maxima are "proper" maxima.  In any case,
    // we have  | T_n(x_k) | = 1  for the maxima and.
    // x_k < x_j  <=>  k < j.
    template<class F> F Tn_maximum (int n , int k); // 0 <= k <= n.

    extern template double Tn_zero (int, int);
    extern template Float  Tn_zero (int, int);

    extern template double Tn_maximum (int, int);
    extern template Float  Tn_maximum (int, int);

}; // ::Cheby

#endif // MP_CHEBY_H
