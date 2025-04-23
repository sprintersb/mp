#include "mp-vector.h"
#include "mp-matrix.h"

/* The algebra stuff herein should be independent of the very scalar type,
 * however we have
 *
 *    #include "mp-float.h"
 *
 * down the line to explicitly instanciate for double and for Float.
 */

#include <cstdio>
#include <cstring>
#include <cassert>

#include "diagnostic.h"

namespace
{
    template<class F> bool _isnumber (const F&) { return true; }
} // anon

template<class F>
FVector<F>::~FVector<F>()
{
    free();
}

template<class F>
void FVector<F>::alloc (int n)
{
    if (n == dim_)
        return;

    if (dim_)
        free();

    assert (v_ == nullptr);

    dim_ = n;
    v_ = n
        ? new element_type[n]
        : nullptr;
}


template<class F>
void FVector<F>::free()
{
    // delete[] will call destructors for all the array elements.
    delete[] v_;
    v_ = nullptr;
    dim_ = 0;
}

template<class F>
FVector<F>::FVector (const V& w)
{
    alloc (w.dim_);
    for (int i = 0; i < dim_; i++)
        v_[i] = w[i];
}

template<class F>
FVector<F>::FVector (V&& w)
{
    std::memcpy (this, &w, sizeof (V));
    w.dim_ = 0;
    w.v_ = nullptr;
}


template<class F>
FVector<F>::FVector (std::initializer_list<F> fs)
{
//info (0, "V init_list %d\n", (int) fs.size());
    alloc (fs.size());
    const F* f = fs.begin();
    for (int i = 0; i < dim_; i++)
        v_[i] = f[i];
}


template<class F>
auto FVector<F>::operator = (const V& w) -> V&
{
    alloc (w.dim_);

    for (int i = 0; i < dim_; i++)
        v_[i] = w[i];

    return *this;
}

template<class F>
auto FVector<F>::operator = (V&& w) -> V&
{
    assert (this != &w);
    free();
    std::memcpy (this, &w, sizeof (V));
    w.dim_ = 0;
    w.v_ = nullptr;

    return *this;
}

template<class F>
FVector<F>::FVector (const std::vector<F>& w)
{
    alloc ((int) w.size());

    for (int i = 0; i < dim_; ++i)
        v_[i] = w[i];
}

template<class F>
FVector<F>::FVector (std::vector<F>&& w)
{
    alloc ((int) w.size());

    for (int i = 0; i < dim_; ++i)
        v_[i] = std::move (w[i]);
}

template<class F>
auto FVector<F>::operator = (const std::vector<F>& w) -> V&
{
    alloc ((int) w.size());

    for (int i = 0; i < dim_; ++i)
        v_[i] = w[i];

    return *this;
}

template<class F>
auto FVector<F>::operator = (std::vector<F>&& w) -> V&
{
    alloc ((int) w.size());

    for (int i = 0; i < dim_; ++i)
        v_[i] = std::move (w[i]);

    return *this;
}

template<class F>
FVector<F>::operator std::vector<F>() const
{
    std::vector<F> w;
    for (int i = 0; i < dim_; ++i)
        w.push_back (v_[i]);

    return w;
}

template<class F>
bool FVector<F>::vector_p () const
{
    for (int i = 0; i < dim_; i++)
        if (! _isnumber<F> (v_[i]))
            return false;
    return true;
}

template<class F>
auto FVector<F>::make (int n, const F& elt) -> V
{
    V z { Alloc (n) };

    for (int i = 0; i < z.dim_; i++)
        z.v_[i] = elt;
    return z;
}

template<class F>
auto FVector<F>::e_k (int dim, int k, const F& k_value) -> V
{
    V w { Alloc (dim) };

    for (int i = 0; i < dim; i++)
        w[i] = (i == k) ? k_value : element_0;

    return w;
}

template<class F>
auto FVector<F>::random_range (int dim, const F& x0, const F& x1) -> V
{
    V z { Alloc (dim) };

    for (int i = 0; i < dim; i++)
        z[i] = ::random_range (x0, x1);

    return z;
}


template<class F>
F& FVector<F>::operator[] (int n)
{
    if (!is_index (n))
        fatal ("FVector %p: bad index %d/%d for operator[]", this, n, dim_);

    return v_[n];
}

template<class F>
const F& FVector<F>::operator[] (int n) const
{
    if (!is_index (n))
        fatal ("FVector %p: bad index %d/%d for operator[] const", this, n, dim_);

    return v_[n];
}


template<class F>
auto FVector<F>::operator - () const -> V
{
    V z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.v_[i] = -v_[i];

    return z;
}

template<class F>
void FVector<F>::operator -= (const V &w)
{
    if (dim_ != w.dim_)
        fatal ("FVector %p -= %p: bad [%d] -= [%d]", this, &w, dim_, w.dim_);

    for (int i = 0; i < dim_; i++)
        v_[i] -= w.v_[i];
}

template<class F>
void FVector<F>::operator += (const V &w)
{
    if (dim_ != w.dim_)
        fatal ("FVector %p += %p: bad [%d] += [%d]", this, &w, dim_, w.dim_);

    for (int i = 0; i < dim_; i++)
        v_[i] += w.v_[i];
}

template<class F>
void FVector<F>::operator *= (const F& f)
{
    for (int i = 0; i < dim_; i++)
        v_[i] *= f;
}

template<class F>
void FVector<F>::operator /= (const F& f)
{
    for (int i = 0; i < dim_; i++)
        v_[i] /= f;
}

template<class F>
void FVector<F>::sub_scaled (const F& f, const FVector<F>& w)
{
    if (dim_ != w.dim_)
        fatal ("FVector %p.sub_scaled (%f, %p): bad [%d] != [%d]",
               this, (double) f, &w, dim_, w.dim_);

    for (int i = 0; i < dim_; i++)
        v_[i] -= f * w[i];
}

template<class F>
auto FVector<F>::operator - (const V& w) const -> V
{
    if (dim_ != w.dim_)
        fatal ("FVector %p - %p: bad [%d] - [%d]", this, &w, dim_, w.dim_);

    V z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.v_[i] = v_[i] - w.v_[i];

    return z;
}

template<class F>
auto FVector<F>::operator + (const V& w) const -> V
{
    if (dim_ != w.dim_)
        fatal ("FVector %p + %p: bad [%d] + [%d]", this, &w, dim_, w.dim_);

    V z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.v_[i] = v_[i] + w.v_[i];

    return z;
}

template<class F>
auto FVector<F>::operator * (const F& f) const -> V
{
    V w { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        w.v_[i] = v_[i] * f;

    return w;
}

template<class F>
FVector<F> operator * (const F& f, const FVector<F>& w)
{
    FVector<F> z { Alloc (w.dim_) };

    for (int i = 0; i < w.dim_; i++)
        z.v_[i] = f * w.v_[i];

    return z;
}

template<class F>
auto FVector<F>::operator / (const F& f) const -> V
{
    V z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.v_[i] = v_[i] / f;

    return z;
}

template<class F>
F FVector<F>::operator * (const V& w) const
{
    if (dim_ != w.dim_)
        fatal ("FVector %p * %p: bad [%d] * [%d]", this, &w, dim_, w.dim_);

    F s {0};
    for (int i = 0; i < dim_; i++)
        s += v_[i] * w.v_[i];
    return s;
}

template<class F>
F FVector<F>::min() const
{
    F m { -element_inf };
    for (int i = 0; i < dim_; i++)
        m = fmin (m, v_[i]);
    return m;
}

template<class F>
F FVector<F>::max() const
{
    F m { -element_inf };
    for (int i = 0; i < dim_; i++)
        m = fmax (m, v_[i]);
    return m;
}

template<class F>
F FVector<F>::height() const
{
    F m { -element_inf };
    for (int i = 0; i < dim_; i++)
        m = fmax (m, fabs (v_[i]));
    return m;
}


template<class F>
auto FVector<F>::dyade () const -> M
{
    M z { Alloc (dim_) };
    V row;

    for (int i = 0; i < dim_; i++)
    {
        row.alloc (dim_);
        for (int c = 0; c < dim_; c++)
        {
            row[c] = (i > c)
                ? z[c][i]
                : v_[i] * v_[c];
        }
        z[i] = std::move (row);
    }

    return z;
}


template<class F>
auto FVector<F>::dyade (const V& w) const -> M
{
    if (this == &w)
        return dyade();

    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z[i] = v_[i] * w;

    return z;
}

template<class F>
F FVector<F>::abs2 () const
{
    return (*this) * (*this);
}

template<class F>
F FVector<F>::abs () const
{
    return sqrt (abs2());
}

template<class F>
auto FVector<F>::normed () const -> V
{
    return (*this) * (element_1 / this->abs());
}

// Project the vector onto  y  and return the length of the resulting vector
// in units of the length of y.  */
template<class F>
auto FVector<F>::project_on (const V& y) const -> F
{
    return (*this) * y / y.abs2();
}


/*
    Map vector  x  by reflection  H  to a given direction  y:  Hx = a * y
    with |a*y| = |x|.  The plane of reflection runs through the origin and
    is orthogonal to  v = x - Hx.

    Return vector  v  with  |v| = 1.

    The linear mapping  H  has the following properties.  Denoting the
    dyadic product (a  dim x dim  matrix) as  D(v,v):

    * H = Id - 2 D(v,v) / <v,v>

      For the matter of brevity, let |y| = 1, then Hx = a * y where  a  is
      a scalar satisfying  |a| = |x|, then

          v = (x - a*y) / |x - a*y|

    * H has one eigenvalue -1 and dim-1 eigenvalues 1.  The eigenvector
      for -1 is  v, and the other eigenvectors span the orthogonal
      complement of  v.

    * H^2 = Id, H is symmteric and orthogonal: H = H^T, H^T * H = Id.


    For numerical stability we chose  a  in such a way that  a*<x,y> <= 0 to
    make the denominator as big as possible.
*/

template<class F>
auto FVector<F>::reflect_to (const V& y) const -> V
{
    if (dim_ != y.dim_)
        fatal ("FVector %p.reflect_to(%p): [%d].reflect_to[%d]", this, &y,
               dim_, y.dim_);

    const V& x = *this;

    F xy { x * y };
    F aa { sqrt ((x*x) / (y*y)) };
    if (xy > element_0)
        aa = -aa;

    V v { x - aa * y };
    assert (v.abs2() >= (x + aa * y).abs2());

    return v.normed();
}

/*  A reflection matrix that reflects a vector to direction W.
    "Reflection" means that one eigenvalue is -1 and the others are 1.
    The result is orthogonal and symmetric.

    In practice, we don't want to use Householder matrices because the very
    transform can be performed with less overhead.  In particular,
    Dyad(v,v)*x = <v,x> * v  where the right side needs only 2*N multiplies
    whereas the left side needs 2*N^2.  */

template<class F>
auto FVector<F>::householder (const V& w) const -> M
{
    if (dim_ != w.dim_)
        fatal ("FVector %p.householder(%p): [%d].householder[%d]", this, &w,
               dim_, w.dim_);

    V v { reflect_to (w) };

    return M::id (dim_) - F{2} * v.dyade();
}


// Multiply tails of 2 vectors starting at k.
// Empty vector for k = dim_ will return 0.
template<class F>
F FVector<F>::tail_mul_tail (int k, const V &w) const
{
    assert (dim_ == w.dim_ && k >= 0 && k <= dim_);

    F s { element_0 };
    for (int i = k; i < dim_; i++)
        s += v_[i] * w[i];

    return s;
}


// Multiply tail starting at k with w.
template<class F>
F FVector<F>::tail_mul (int k, const V &w) const
{
    assert (dim_ == k + w.dim_);

    F s { element_0 };
    for (int i = 0; i < w.dim_; i++)
        s += v_[k + i] * w[i];

    return s;
}

// Multiply with the tail of column c, tail starting at [k][c].
template<class F>
F FVector<F>::mul_column_tail (int k, const M& m, int col) const
{
    int d = m.dim() - k;
    if (dim_ != d)
        fatal ("FVector %p mul_column_tail (%d, %p, %d): bad [%d] * [%d]-%d",
               this, k, &m, col, dim_, m.dim(), k);

    F s { element_0 };

    for (int i = 0; i < d; i++)
        s += v_[i] * m[k + i][col];

    return s;
}


template<class F>
F FVector<F>::mul_column (const M& m, int col) const
{
    if (dim_ != m.dim())
        fatal ("FVector %p mul_column (%p, %d): bad [%d] * [%d]", this, &m, col,
               dim_, m.dim());

    return mul_column_tail (0, m, col);
}


template<class F>
std::ostream& operator << (std::ostream& ost, const FVector<F>& w)
{
    ost << " (";
    for (int i = 0; i < w.dim_; i++)
        ost << " " << w[i];
    ost << " ) ";
    return ost;
}


template<class F> const F FVector<F>::element_0 { 0 };
template<class F> const F FVector<F>::element_1 { 1 };
template<> const double FVector<double>::element_nan { std::nan("") };
template<> const Float  FVector<Float>::element_nan {  };
template<> const double FVector<double>::element_inf { __builtin_huge_val() };
template<> const Float  FVector<Float>::element_inf { __builtin_huge_val() };

namespace
{
    template<> bool _isnumber (const Float& x) { return !x.inf_p() && !x.nan_p(); }
    template<> bool _isnumber (const double& x) { return !std::isinf(x) && !std::isnan(x); }
} // anon



// Explicit instanciate Float.

template const Float FVector<Float>::element_nan;
template FVector<Float> operator * (const Float&, const FVector<Float>&);
template std::ostream& operator << (std::ostream&, const FVector<Float>&);

template class FVector<Float>;

// Explicit instanciate double.

template const double FVector<double>::element_nan;
template FVector<double> operator * (const double&, const FVector<double>&);
template std::ostream& operator << (std::ostream&, const FVector<double>&);

template class FVector<double>;

