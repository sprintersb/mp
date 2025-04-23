#ifndef MP_CUBIC_SPLINE_H
#define MP_CUBIC_SPLINE_H

#include "mp-poly.h"

#include <vector>

extern void test_cubic_spline();

template<class F>
class CubicSpline
{
private:
    CubicSpline () = default;

public:
    typedef F scalar_type;
    typedef Poly<F> polynomial_type;

    int n_points_ = 0, n_polys_ = 0;
    std::vector<F> xs_, ys_;
    std::vector<polynomial_type> polys_;
    // The linear functions that smoothly connect to the spline outside
    // of the range given by the x's.
    Poly<F> poly_l_, poly_r_;

    // Construct a cubic spline that smoothly (to 2nd order) connects
    // the points (x_i, y_i).  The number of x-coordinates and y-coordinates
    // must be the same, and the x's must be strictly increasing.
    CubicSpline (const std::vector<F>& xs, const std::vector<F>& ys);

    // Evaluate the spline at x.  If x is outside the range given by the x's,
    // then use a linear function that smoothly connects to the spline.
    F operator () (const F& x) const;

    // Get the index of the polynomial that's responsible for x-value x.
    // Return -1 if x is smaller than the smallest of x's.
    // Return -2 if x is larger than the largest of x's.
    int index (const F& x) const;

    // A (const) reference to the polynomial at the specified index.
    // i may be -1 or -2 as returned by .index(), in which case a linear
    // polynomial is referenced that smoothly continues outside the range
    // of the very spline.
    Poly<F>& operator [] (int i);
    const Poly<F>& operator [] (int i) const;

    // A (const) reference to the polynomial that is responsible for mapping x.
    // If x is outside the range as specified by the x's, then return the
    // linear continuation.
    Poly<F>& poly_at (const F& x);
    const Poly<F>& poly_at (const F& x) const;

    // The 1st derivative.  In general, this is no more a cubic spline.
    CubicSpline<F> D() const;
};

extern template class CubicSpline<double>;

#endif // MP_CUBIC_SPLINE_H
