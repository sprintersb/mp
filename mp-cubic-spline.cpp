#include "mp-cubic-spline.h"

#include "mp-matrix.h"
#include "diagnostic.h"

#include <cassert>

template<class F>
CubicSpline<F>::CubicSpline (const std::vector<F>& xs,
                             const std::vector<F>& ys)
    : n_points_((int) xs.size()), n_polys_(n_points_-1), xs_(xs), ys_(ys)
{
    assert (xs.size() == ys.size());
    assert (xs.size() >= 2);

    /* We number the polynomials from  0  to  n_polys-1, and FVector v[] holds
       the coefficients of the polynomials:
            v[0] = poly[0].d
            v[1] = poly[0].c
            v[2] = poly[0].b
            v[3] = poly[0].a
            v[4] = poly[1].d ...
       where the polynomials have the form
            poly  = a*x^3 + b*x^2 + c*x + d,
            poly' = 3a*x^2 + 2b*x + c,
            poly" = 6a*x + 2b.
       The constraints on the coefficients are (with k = i+1):
       (1) poly[i](x_i)  = y_i             # n_polys eqs, i in [0, n_polys-1].
       (2) poly[i](x_k)  = y_k             # n_polys eqs, i in [0, n_polys-1].
       (3) poly[i]'(x_k) = poly[k]'(x_k)   # n_polys-1 eqs, i in [0, n_polys-2].
       (4) poly[i]"(x_k) = poly[k]"(x_k)   # n_polys-1 eqs, i in [0, n_polys-2].
       (5) poly[0]"(x_0) = 0               # 1 eq.
       (6) poly[n_polys-1]"(x_n_polys) = 0 # 1 eq.

       This are 4*n_polys equations for the 4*n_polys unknown coefficients.
       (1) and (2) state that the polynomials interpolate over the respective
       intervals given by adjacent x's.  (3) states that the transition at the
       junction points is smooth of 1st order at least, and (4) states that the
       transitions are smooth of 2nd order at least.  (5) states that left of
       the leftmost x, the graph can be continued by a line and that such a
       continuation can be made smoothly of 2nd order.  (6) states a similar
       condition like (5) for the area right of the rightmost x.  */

    auto m = FMatrix<F>::id (4*n_polys_, FMatrix<F>::element_0);
    auto b = FVector<F>::make (4*n_polys_);

    for (int j = 0; j < n_polys_; ++j)
    {
        int i = 4*j, r;
        int k = 4*(j+1);
        const F& x0 = xs[j+0];
        const F& x1 = xs[j+1];
        assert (x0 < x1);

        // (1) poly[i](x_i) = y_i       # n_polys eqs, i in [0, n_polys-1].
        r = j + 0*n_polys_;
        b[r] = ys[j+0];
        m[r][i + 0] = F{1};
        m[r][i + 1] = x0;
        m[r][i + 2] = x0 * x0;
        m[r][i + 3] = x0 * x0 * x0;

        // (2) poly[i](x_k) = y_k       # n_polys eqs, i in [0, n_polys-1].
        r = j + 1*n_polys_;
        b[r] = ys[j+1];
        m[r][i + 0] = F{1};
        m[r][i + 1] = x1;
        m[r][i + 2] = x1 * x1;
        m[r][i + 3] = x1 * x1 * x1;

        if (j != n_polys_ -1)
        {
            // (3) poly[i]'(x_k) = Poly[k]'(x_k) # n_polys-1 eqs,
            //     i in [0, n_polys-2].
            r = j + 2*n_polys_;
            b[r] = F{0};
            //m[r][i + 0] = F{0};
            m[r][i + 1] = F{1};
            m[r][i + 2] = F{2} * x1;
            m[r][i + 3] = F{3} * x1 * x1;
            //m[r][k + 0] = - m[r][i + 0];
            m[r][k + 1] = - m[r][i + 1];
            m[r][k + 2] = - m[r][i + 2];
            m[r][k + 3] = - m[r][i + 3];
            
            // (4) poly[i]"(x_k) = Poly[k]"(x_k) # n_polys-1 eqs,
            //     i in [0, n_polys-2].
            r = j + 3*n_polys_;
            b[r] = F{0};
            //m[r][i + 0] = F{0};
            //m[r][i + 1] = F{0};
            m[r][i + 2] = F{2};
            m[r][i + 3] = F{6} * x1;
            //m[r][k + 0] = - m[r][i + 0];
            //m[r][k + 1] = - m[r][i + 1];
            m[r][k + 2] = - m[r][i + 2];
            m[r][k + 3] = - m[r][i + 3];
        }
        else
        {
            // (5) poly[0]"(x_0) = 0                    # 1 eq.
            r = j + 2*n_polys_;
            b[r] = F{0};
            //m[r][0 + 0] = F{0};
            //m[r][0 + 1] = F{0};
            m[r][0 + 2] = F{2};
            m[r][0 + 3] = F{6} * xs[0];
            
            // (6) poly[n_polys-1]"(x_n_polys) = 0      # 1 eq.
            r = j + 3*n_polys_;
            b[r] = F{0};
            //m[r][i + 0] = F{0};
            //m[r][i + 1] = F{0};
            m[r][i + 2] = F{2};
            m[r][i + 3] = F{6} * x1;
        }
    }

    //cout << b << endl;
    //cout << m;

    FMatrix<F>::push_verbose (0);
    FVector<F> v = m.solve (b, Algo::QR_vectors);
    FMatrix<F>::pop_verbose();

    //cout << v << endl;

    polys_.resize (n_polys_);

    for (int j = 0; j < n_polys_; ++j)
    {
        int i = 4*j;
        Poly<F> p { v[i+0], v[i+1], v[i+2], v[i+3] };
        polys_[j] = p;
    }

    const F& x0 = xs_[0];
    const F& y0 = ys_[0];
    F a0 = polys_[0].D() (x0);
    poly_l_ = Poly<F> { y0 - a0*x0, a0 };

    const F& x1 = xs_[n_points_ - 1];
    const F& y1 = ys_[n_points_ - 1];
    F a1 = polys_[n_polys_ - 1].D() (x1);
    poly_r_ = Poly<F> { y1 - a1*x1, a1 };
}


template<class F>
int CubicSpline<F>::index (const F& x) const
{
    assert (xs_.size() == ys_.size());
    assert (xs_.size() >= 2);
    assert (n_points_ == (int) xs_.size());
    assert (n_polys_ == (int) polys_.size());
    assert (n_polys_ == n_points_ - 1);

    int lo = 0, hi = n_polys_ - 1;

    if (x < xs_[lo])
        return -1;
    else if (x > xs_[hi+1])
        return -2;

    while (lo != hi)
    {
        assert (lo < hi && lo >= 0 && hi < n_polys_);
        int mid = (lo + hi) / 2;
        assert (mid >= 0 && mid < n_polys_);

        if (x < xs_[mid])
            hi = mid - 1;
        else if (x > xs_[mid + 1])
            lo = mid + 1;
        else
            return mid;
    }

    return lo;
}

template<class F>
F CubicSpline<F>::operator () (const F& x) const
{
    return poly_at(x) (x);
}

template<class F>
Poly<F>& CubicSpline<F>::operator [] (int i)
{
    assert (i >= -2 && i < n_polys_);
    return i == -1 ? poly_l_ : i == -2 ? poly_r_ : polys_[i];
}

template<class F>
const Poly<F>& CubicSpline<F>::operator [] (int i) const
{
    assert (i >= -2 && i < n_polys_);
    return i == -1 ? poly_l_ : i == -2 ? poly_r_ : polys_[i];
}

template<class F>
Poly<F>& CubicSpline<F>::poly_at (const F& x)
{
    return this->operator[] (index(x));
}

template<class F>
const Poly<F>& CubicSpline<F>::poly_at (const F& x) const
{
    return this->operator[] (index(x));
}

template<class F>
CubicSpline<F> CubicSpline<F>::D () const
{
    CubicSpline d;
    d.n_points_ = n_points_;
    d.n_polys_ = n_polys_;

    d.poly_l_ = poly_l_.D();
    d.poly_r_ = poly_r_.D();
    d.polys_.resize (n_polys_);
    for (int i = 0; i < n_polys_; ++i)
        d.polys_[i] = polys_[i].D();

    d.xs_ = xs_;
    d.ys_.resize (n_points_);
    d.ys_[0] = d.polys_[0] (xs_[0]);
    for (int i = 0; i < n_polys_; ++i)
        d.ys_[i+1] = d.polys_[i] (xs_[i+1]);

    return d;
}

///////////////////////////////////////////////////////////////////////
// Explicit instanciate extern templates for double.

template class CubicSpline<double>;

///////////////////////////////////////////////////////////////////////

#include <iostream>

using std::cout;
using std::endl;

void test_cubic_spline()
{
    info (0, "%s", __PRETTY_FUNCTION__);
    using CS = CubicSpline<double>;
    using Vs = std::vector<CS::scalar_type>;

    Vs xs = Vs{1, 2, 3, 4};
    Vs ys = Vs{1, 3, 2, 5};
    
    auto s = CS (xs, ys);
    const CS::scalar_type& x0 = s.xs_[0];
    const CS::scalar_type& x1 = s.xs_[s.n_points_ - 1];

    for (CS::scalar_type x = x0 - 1; x <= x1 + 1; x += 0.5)
        printf ("index(%.2f) = %d\n", x, s.index(x));

    for (CS::scalar_type x = x0 - 1; x <= x1 + 1; x += 0.1)
        out ("%.3f -> %+.3f\n", x, s(x));

    auto ds = s.D();
    for (CS::scalar_type x = x0 - 1; x <= x1 + 1; x += 0.1)
        printf ("D %.3f -> %+.3f\n", x, ds(x));

    auto d2s = ds.D();
    for (CS::scalar_type x = x0 - 1; x <= x1 + 1; x += 0.1)
        printf ("DD %.3f -> %+.3f\n", x, d2s(x));

    auto d3s = d2s.D();
    for (CS::scalar_type x = x0 - 1; x <= x1 + 1; x += 0.1)
        printf ("DDD %.3f -> %+.3f\n", x, d3s(x));

    for (int j = 0; j < s.n_polys_; ++j)
    {
        const CS::polynomial_type& p = s.polys_[j];

        printf ("p[%d] = ", j); cout << p << endl;
        printf ("p[%d](x0) = ", j); cout << p(xs[j+0]) << endl;
        printf ("p[%d](x1) = ", j); cout << p(xs[j+1]) << endl;
        printf ("  p[%d]'(x0) = ", j); cout << p.D()(xs[j+0]) << endl;
        printf ("  p[%d]'(x1) = ", j); cout << p.D()(xs[j+1]) << endl;
        printf ("    p[%d]\"(x0) = ", j); cout << p.D().D()(xs[j+0]) << endl;
        printf ("    p[%d]\"(x1) = ", j); cout << p.D().D()(xs[j+1]) << endl;
    }
}
