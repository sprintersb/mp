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

const char *const Algo::name[] =
{
    [QR_slow] = "QR.Slow",
    [QR_fast] = "QR.Fast",
    [QR_vectors] = "QR.Fastest",
};


template<class F>
void FMatrix<F>::apply (void (*func)(M&, int, int, void*), void* data)
{
    for (int i = 0; i < dim_; i++)
        for (int c = 0; c < r_[i].dim_; c++)
            func (*this, i, c, data);
}


template<class F>
FMatrix<F>::~FMatrix<F>()
{
    free();
}

template<class F>
void FMatrix<F>::alloc (int n)
{
    if (n == dim_)
        return;

    if (dim_)
        free();

    assert (r_ == nullptr);

    dim_ = n;
    r_ = n
        ? new row_type[n]
        : nullptr;
}


template<class F>
void FMatrix<F>::alloc (int dim, int v_dim)
{
//info (0, "M alloc(%d,%d)\n", dim, v_dim);
    if (dim < 0)
        fatal ("FMatrix.alloc(%d,%d))", dim, v_dim);

    alloc (dim);
    if (v_dim >= 0)
        for (int i = 0; i < dim_; i++)
            r_[i].alloc (v_dim);
}

template<class F>
void FMatrix<F>::free()
{
    // delete[] will call destructors for all array elements.
    delete[] r_;
    r_ = nullptr;
    dim_ = 0;
}

template<class F>
FMatrix<F>::FMatrix (const M& m)
{
    alloc (m.dim_);
    for (int i = 0; i < dim_; i++)
        r_[i] = m[i];
}

template<class F>
FMatrix<F>::FMatrix (M&& m)
{
    free();
    std::memcpy (this, &m, sizeof (M));
    m.dim_ = 0;
    m.r_ = nullptr;
}


template<class F>
FMatrix<F>::FMatrix (std::initializer_list<V> vs)
{
//info (0, "M init_list %d\n", (int) ms.size());
    alloc (vs.size());
    const V* v = vs.begin();
    for (int i = 0; i < dim_; i++)
        r_[i] = v[i];
    if (n_cols() < 0)
        fatal ("FMatrix <initializer_list>[%d]: row[%d] = %d, row[0] = %d", dim_,
               -n_cols(), r_[-n_cols()].dim_, r_[0].dim_);
}

template<class F>
auto FMatrix<F>::operator = (const M& m) -> M&
{
    alloc (m.dim_);

    for (int i = 0; i < dim_; i++)
        r_[i] = m[i];

    return *this;
}

template<class F>
auto FMatrix<F>::operator = (M&& m) -> M&
{
    free();
    std::memcpy (this, &m, sizeof (M));
    m.dim_ = 0;
    m.r_ = nullptr;

    return *this;
}


template<class F>
int FMatrix<F>::dim () const
{
    return dim_;
}

template<class F>
bool FMatrix<F>::matrix_p () const
{
    for (int i = 0; i < dim_; i++)
        if (! r_[i].vector_p ())
            return false;
    return true;
}

template<class F>
int FMatrix<F>::n_cols () const
{
    if (dim_ == 0)
        return 0;

    int cc = r_[0].dim_;

    for (int i = 1; i < dim_; i++)
        if (r_[i].dim_ != cc)
            return -i;
    return cc;
}


template<class F>
auto FMatrix<F>::make (int dim, int cols, const F& elt) -> M
{
    M z { Alloc (dim) };

    if (cols > 0)
    {
        V row { V::make (cols, elt) };
        for (int i = 0; i < dim; i++)
            z.r_[i] = row;
    }
    return z;
}

template<class F>
auto FMatrix<F>::id (int r, const F& diag) -> M
{
    M m { Alloc (r) };
    //info (10, "::id diag = ...");
    //std::cout << "... " << diag << std::endl;
    for (int i = 0; i < r; i++)
        m[i] = V::e_k (r, i, diag);
    return m;
}


template<class F>
auto FMatrix<F>::random_range (int dim, int cols, const F& x0, const F& x1) -> M
{
    M z { Alloc (dim) };

    for (int i = 0; i < dim; i++)
        z[i] = V::random_range (cols, x0, x1);

    return z;
}

template<class F>
auto FMatrix<F>::operator[] (int r) -> M::row_type&
{
    if (!is_index (r))
        fatal ("FMatrix %p: bad index %d/%d for operator[]", this, r, dim_);

    return r_[r];
}


template<class F>
auto FMatrix<F>::operator[] (int r) const -> const M::row_type&
{
    if (!is_index (r))
        fatal ("FMatrix %p: bad index %d/%d for operator[] const", this, r, dim_);

    return r_[r];
}

template<class F>
auto FMatrix<F>::column (int c, const char *caller) const -> V
{
    V w { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
    {
        if (r_[i].dim_ <= c)
            fatal ("FMatrix%s %p column(%d): row[%d]: [%d] <= %d", caller,
                   this, c, i, r_[i].dim_ , c);
        w.v_[i] = r_[i].v_[c];
    }

    return w;
}

template<class F>
auto FMatrix<F>::T() const -> M
{
    int z_rows = n_cols();
    if (z_rows < 1)
        fatal ("FMatrix::T() %p row[%d].dim = %d, row[0].dim = %d", this,
               -z_rows, r_[-z_rows].dim_, r_[0].dim_);

    M z { Alloc (z_rows) };

    for (int i = 0; i < z.dim_; i++)
        z.r_[i] = column (i, "::T()");

    return z;
}

template<class F>
auto FMatrix<F>::operator - () const -> M
{
    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.r_[i] = - r_[i];

    return z;
}


template<class F>
void FMatrix<F>::operator -= (const M &m)
{
    int s_cols = n_cols();
    int m_cols = m.n_cols();
    if (dim_ != m.dim_ || s_cols < 1 || s_cols != m_cols)
        fatal ("FMatrix %p -= %p: bad [%d][%d] -= [%d][%d]", this, &m,
               dim_, s_cols, m.dim_, m_cols);

    for (int i = 0; i < dim_; i++)
        r_[i] -= m.r_[i];
}


template<class F>
void FMatrix<F>::operator += (const M &m)
{
    int s_cols = n_cols();
    int m_cols = m.n_cols();
    if (dim_ != m.dim_ || s_cols < 1 || s_cols != m_cols)
        fatal ("FMatrix %p += %p: bad [%d][%d] += [%d][%d]", this, &m,
               dim_, s_cols, m.dim_, m_cols);

    for (int i = 0; i < dim_; i++)
        r_[i] += m.r_[i];
}


template<class F>
auto FMatrix<F>::operator - (const M &m) const -> M
{
    int s_cols = n_cols();
    int m_cols = m.n_cols();
    if (dim_ != m.dim_ || s_cols < 1 || s_cols != m_cols)
        fatal ("FMatrix %p - %p: bad [%d][%d] - [%d][%d]", this, &m,
               dim_, s_cols, m.dim_, m_cols);

    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.r_[i] = r_[i] - m.r_[i];

    return z;
}


template<class F>
auto FMatrix<F>::operator + (const M &m) const -> M
{
    int s_cols = n_cols();
    int m_cols = m.n_cols();
    if (dim_ != m.dim_ || s_cols < 1 || s_cols != m_cols)
        fatal ("FMatrix %p + %p: bad [%d][%d] + [%d][%d]", this, &m,
               dim_, s_cols, m.dim_, m_cols);

    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.r_[i] = r_[i] + m.r_[i];

    return z;
}


template<class F>
auto FMatrix<F>::operator * (const M &m) const -> M
{
    int m_cols = m.n_cols();
    int s_cols = n_cols();
    if (s_cols < 1 || m_cols < 1 || s_cols != m.dim_)
        fatal ("FMatrix %p * %p: bad [%d][%d] * [%d][%d]", this, &m,
               dim_, s_cols, m.dim_, m_cols);

    // dim_ of the result remains the same, but we must use a fresh intermediate
    // matrix due to potential early-clobber.

    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
    {
        V& row = z[i];
        row.alloc (m_cols);

        for (int c = 0; c < m_cols; c++)
        {
            if (r_[i].dim_ != m.dim_)
                fatal ("FMatrix %p *= %p: bad row[%d] * col[%d] = [%d] * [%d]",
                       this, &m, i, c, r_[i].dim_, m.dim_);

            row.v_[c] = r_[i].mul_column (m, c);
        }
    }

    return z;
}

template<class F>
void FMatrix<F>::operator *= (const M &m)
{
    *this = (*this) * m;
}


template<class F>
auto FMatrix<F>::operator * (const FVector<F> &v) const -> V
{
    int s_cols = n_cols();
    if (s_cols < 1 || s_cols != v.dim_)
        fatal ("FMatrix %p * %p: bad [%d][%d] * [%d]", this, &v,
               dim_, s_cols, v.dim_);

    V z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.v_[i] = r_[i] * v;

    return z;
}


template<class F>
auto FMatrix<F>::operator / (const F& f) const -> M
{
    int s_cols = n_cols();
    if (s_cols < 1)
        fatal ("FMatrix %p * %p: bad [%d][%d] / F", this, &f,
               dim_, s_cols);

    M z { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
        z.r_[i] = r_[i] / f;

    return z;
}

template<class F>
FMatrix<F> operator * (const F& f, const FMatrix<F>& m)
{
    int m_cols = m.n_cols();
    if (m_cols < 1)
        fatal ("FMatrix %p * %p: bad F * [%d][%d]", &f, &m, m.dim_, m_cols);

    FMatrix<F> z { Alloc (m.dim_) };

    for (int i = 0; i < m.dim_; i++)
        z.r_[i] = f * m.r_[i];

    return z;
}


template<class F>
std::ostream& operator << (std::ostream& out, const FMatrix<F>& m)
{
    for (int i = 0; i < m.dim_; i++)
        out << "( " << m[i] << " )" << std::endl;
    return out;
}


// In the square  d x d  minor  S  starting at  [k][k], d = dim - k:
// Get a vector  v  to be used in a Householder reflection that maps
// the first row of  S  to the direction of the 1st unit vector in R^d.
// The Householder transform for  v  will be Id(d) - 2*Dyade(v,v).
// See FVector::reflect_to() for a detailed description.

template<class F>
auto FMatrix<F>::reflect_to_e_k (int k) const -> V
{
    if (k >= dim_)
        fatal ("FMatrix %p.reflect_to(%d): %d >= %d", this, k, k, dim_);

    const M& m = *this;
    int d = dim_ - k;
    // x * e_1  where  x  is the vector below m[k][k] and  e_1  is the
    // 1st unit vector in  R^d, d=dim-k.
    F x_e1 { m[k][k] };

    F aa { element_0 };
    // Compute the norm of x.
    for (int i = 0; i < d; i++)
        aa += m[k+i][k] * m[k+i][k];
    aa = sqrt (aa);

    // Make sure the following division is with denominator as big as it gets.
    if (x_e1 > element_0)
        aa = -aa;
    // Compute v = x - aa * e_k
    V v { Alloc (d) };

    v[0] = m[k][k] - aa;  // aa * e_k has only 1 non-0 component, the 1st.
    for (int i = 1; i < d; i++)
        v[i] = m[k+i][k];

    return v.normed();
}


// Householder from  d x d  Minor starting at  [k][k],  d = dim_ - k.
template<class F>
auto FMatrix<F>::householder (int k) const -> M
{
    int d = dim_ - k;
    V v { this->reflect_to_e_k (k) };

    return M::id (d) - F{2} * v.dyade(v);
}

// Multiply from the right with a smaller symmetric matrix M that starts
// at [k][k] and fits into the lower right corner.  M is extended implicitly
// to a  dim_ x dim_  matrix as follows:  West and North of M are zeros.
// North-West is a identity matrix of dimension k x k.

template<class F>
void FMatrix<F>::mul_minor_right (int k, const M& m)
{
    int d = m.dim_;
    assert (k == dim_ - d);

    // The left k columns 0...k-1 are unaltered.
    // The  dim_ x d  matrix resulting from the remaining columns is multiplied
    // as usual to get a  dim_ x d = (dim_ x d) * (d x d)  matrix that replaces
    // the right columns.  In addition we know that M is symmetric.

    V row { Alloc (dim_) };

    for (int i = 0; i < dim_; i++)
    {
        //info (0, "%d", i);
        const V& ri = r_[i];
        for (int c = 0; c < k; c++)
            row[c] = ri[c];
        //info (0, "%d", i);
        for (int c = k; c < dim_; c++)
        {
            //info (0, "%d,%d", i, c);
            row[c] = ri.tail_mul (k, m[c-k]);
        }
        r_[i] = row;
    }
}


// Multiply from the left with a smaller symmetric matrix M that starts
// at [k][k] and fits into the lower right corner.  M is extended implicitly
// to  a  dim_ x dim_  matrix as follows:  West and North of M are zeros.
// North-West is a identity matrix of dimension k x k.  In addition, we
// may assume that the first k columns of *this are zero below the diagonal.
// This implies that the 1st k columns of *this and the 1st k rows of *this
// are unaltered and all what's happening is in the lower right  d x d  part.

template<class F>
void FMatrix<F>::mul_minor_left (int k, const M& m)
{
    int d = m.dim_;
    assert (k == dim_ - d);

    // The top k rows 0...k-1 are unaltered, same for the left k columns.
    // The remaining  d x d  matrix in the lower right corner is multiplied
    // as usual.

    M A { Alloc (d, d) };

    for (int i = 0; i < d; i++)
    {
        const V& row = m[i];
        for (int c = 0; c < d; c++)
            A[i][c] = row.mul_column_tail (k, *this, c+k);
    }

    for (int i = 0; i < d; i++)
        for (int c = 0; c < d; c++)
            (*this)[i+k][c+k] = A[i][c];
}


// Like  void mul_minor_left (int k, const FMatrix<F>& H)
// with  H = Id(d) - 2 * Dyade(v,v)  with  d = dim(v) = dim - k.
// Notice that if  c  is a column vector, then  Dyade(v,v) * c = <v,c> * v
// is the resulting column vector.

template<class F>
void FMatrix<F>::mul_householder_minor_left (int k, const V& v)
{
    int d = v.dim_;
    assert (k == dim_ - d && k >= 0);

    // The top k rows 0...k-1 are unaltered, same for the left k columns
    // because the assertion is that left of  H  are only zeroes.
    // The remaining  d x d  matrix in the lower right corner is multiplied
    // from the left with H as usual.

    M A { Alloc(d) };
    F s, m2{-2};

    for (int c = 0; c < d; c++)
    {
        s = m2 * v.mul_column_tail (k, *this, c+k);
        //info (0, "mul_column_tail [%d].(%d, [%d] %d)\n", v.dim_, k, dim_, c+k);

        A[c] = s * v; // Will be transposed implicitly below.
    }

    for (int i = 0; i < d; i++)
    {
        V& row = r_[i+k];
        for (int c = 0; c < d; c++)
            row[c+k] += A[c][i];
    }
}


// Like  void mul_minor_right (int k, const M& H)
// with  H = Id(d) - 2 * Dyade(v,v)  with  d = dim(v) = dim - k.
// Notice that if  r  is a row vector, then  r * Dyade(v,v) = <r,v> * v
// is the resulting row vector.

template<class F>
void FMatrix<F>::mul_householder_minor_right (int k, const V& v)
{
    int d = v.dim_;
    assert (k == dim_ - d && k >= 0);

    // The left k columns 0...k-1 are unaltered.  The remaining  dim x d  matrix
    // in the right is multiplied from the right with H as usual.

    F s, m2{-2};

    for (int i = 0; i < dim_; i++)
    {
        V& row = r_[i];
        s = m2 * row.tail_mul (k, v);
        //info (0, "tail_mul [%d].(%d [%d])\n", row.dim_, k, v.dim_);

        for (int c = 0; c < d; c++)
            row[c+k] += s * v[c];
    }
}


template<class F>
void FMatrix<F>::QR_householder_step (int k, M& Q, M& R) const
{
    assert (k >= 0 && k < dim_ - 1);

    if (k == 0)
    {
        //M h { column(0).householder (V::e_k (dim_, 0)) };
        M h { this->householder (0) };

        R = h * (*this);
        Q = h;
    }
    else
    {
        M h { R.householder (k) };

        R.mul_minor_left (k, h);
        Q.mul_minor_right (k, h);
    }

/*    std::cout << "QR: H." << k << "\n" << h << std::endl;
    std::cout << "QR: H." << k << "^2\n" << (h*h) << std::endl;*/
}


template<class F>
void FMatrix<F>::QR_householder_step_fast (int k, M& Q, M& R) const
{
    assert (k >= 0 && k < dim_ - 1);

    if (k == 0)
        R = (*this);

    V v { R.reflect_to_e_k (k) };
    R.mul_householder_minor_left (k, v);

    if (k == 0)
        Q = M::id (dim_) - F{2} * v.dyade();
    else
        Q.mul_householder_minor_right (k, v);
}


// Q ist NOT a matrix but only contains vectors of DIFFERENT dimensions
// that represent Householder reflections of dimension: dim ... 2.
// Q[k] has dimension dim(Q) - k, 0 <= k <= dim_ - 2.

template<class F>
void FMatrix<F>::QR_householder_step_reflection_only (int k, M& Q, M& R) const
{
    assert (k >= 0 && k < dim_ - 1);

    if (k == 0)
    {
        Q.alloc (dim_ - 1);
        R = (*this);
    }

    Q[k] = R.reflect_to_e_k (k);
    R.mul_householder_minor_left  (k, Q[k]);
}


template<class F>
void FMatrix<F>::QR_decompose (M& Q, M& R, Algo::How how) const
{
    if (dim_ == 1)
    {
        Q = M::id(1);
        R = *this;
        return;
    }

    info (verbose, "QR %s: ", Algo::name[how]);
    for (int k = 0; k < dim_ - 1; k++)
    {
        const char *cr = k && k % 30 == 0 ? "\n" : "";
        if (verbose) out ("%s% 4d", cr, dim_ - k);

        switch (how)
        {
            default:
                assert (0);

            case Algo::QR_slow:
                QR_householder_step (k, Q, R);
                break;

            case Algo::QR_fast:
                QR_householder_step_fast (k, Q, R);
                break;

            case Algo::QR_vectors:
                QR_householder_step_reflection_only (k, Q, R);
                break;
        }

        if (how != Algo::QR_vectors)
        {
            assert (Q.dim_ == dim_);
        }
        else
        {
            assert (Q.dim_ == dim_ - 1);
            for (int i = 0; i <= k; i++)
                assert (Q[i].dim_ == dim_ - i);
        }

        assert (R.dim_ == dim_);

        for (int i = k + 1; i < dim_; i++)
        {
            double f = (double) R[i][k];
            if (fabs (f) > 1e-10)
                warning ("weak elimination: R[%d][%d] = %.3e", i, k, f);

            R[i][k] = element_0;
        }

#if 0
        info (0, "QR: Q.%d", k);
        std::cout << Q << std::endl;
        //std::cout << "QR: Q." << k << " * Q." << k << "^T\n" << (Q*Q.T()) << std::endl;
#endif // 0
        //std::cout << "QR: Delta." << k << " = A - Q*R\n" << ((*this) - Q*R) << std::endl;
    }
    if (verbose) out (".\n");
#if 0
    info (0, "QR: R");
    std::cout << R << std::endl;
#endif
    //if (fast == 2)
    //    std::cout << "Q.vectors =\n" << Q << std::endl;
}

// Solve  A*x = b  where A is an upper triangle (below diagonal is assumed 0).
template<class F>
auto FMatrix<F>::solve_triangle (const V& b) const -> V
{
    int s_cols = n_cols();
    if (s_cols < 1 || s_cols != dim_ || b.dim_ != dim_)
        fatal ("FMatrix %p.solve %p: bad [%d][%d].solve ([%d])", this, &b,
               dim_, s_cols, b.dim_);

    F aii_xi;
    V x { Alloc (dim_) };

    info (verbose, "solve_triangle: ");
    for (int i = dim_ - 1; i >= 0; i--)
    {
        int cc = dim_ - 1 - i;
        const char *cr = cc && cc % 30 == 0 ? "\n" : "";
        if (verbose) out ("%s% 4d", cr, dim_ - i);

        aii_xi = b[i] - r_[i].tail_mul_tail (i + 1, x);
        x[i] = aii_xi / r_[i][i];
    }
    if (verbose) out (".\n");
    return x;
}


template<class F>
auto FMatrix<F>::transpose_mul (const V& w) const -> V
{
    int s_cols = n_cols();
    if (dim_ != w.dim_ || s_cols < 1)
        fatal ("FMatrix %p.transpose_mul %p: bad [%d][%d]^T * [%d]", this, &w,
               dim_, s_cols, w.dim_);

    V z { Alloc (s_cols) };

    for (int c = 0; c < s_cols; c++)
        z[c] = w.mul_column (*this, c);

    return z;
}

// Q = *this  is  NOT an actual matrix but an array of vectors:
// dim_Q[i] = dim_Q + 1 - i;  dim_Q = dim_w - 1.
// Each row of Q represents a Householder-Matrix of dim = dim_w
// H(i) = Id(dim_w) - 2 * Q[i].dyade(Q[i])  where  Q[i] is extended
// with i zeros at the start to yield dimension dim_w.
// Each H(i) is orthogonal and symmetric.
// Q = H(0) * H(1) * H(2) * ... * H(dim_w-1)  is the orthogonal matrix of
// dim_w as produced by a QR decomposition.
//
// This routine returns Q^t * w = H(dim_w-1) * ... * H(1) * H(0) * w.
// Due to the structure of the H(i), each multipliction affects only
// the  dim_w - i  tail components  w[i], ... w[dim_w - 1]  of  w.
//
// Again, we make use of  w * dyade(v,v) = <w,v> * v.
template<class F>
auto FMatrix<F>::householder_vectors_mul (const V& w) const -> V
{
    const M& Q = *this;

    if (Q.dim_ != w.dim_ - 1)
        fatal ("FMatrix %p.householder_vectors_mul %p: bad [%d]+1 != [%d]",
               &Q, &w, Q.dim_, w.dim_);
    F f;
    V z { w };

    info (verbose, "HH_vec_mul: ");
    for (int i = 0; i < w.dim_ - 1; i++)
    {
        const char *cr = i && i % 30 == 0 ? "\n" : "";
        if (verbose) out ("%s% 4d", cr, z.dim_ - i);

        if (w.dim_ - i != Q[i].dim_)
            fatal ("FMatrix %p[%d].householder_vectors_mul %p: bad (%d,[%d]) * [%d]",
                   &Q, i, &w, i, w.dim_, Q[i].dim_);

        f = F{2} * z.tail_mul (i, Q[i]);

        for (int k = i; k < w.dim_; k++)
            z[k] -= f * Q[i][k - i];
    }
    if (verbose) out (".\n");

    return z;
}

// Solve M*z = w  as follows:
// 1) Compute the QR-decomposition of M:  M = Q*R  with an orthogonal matrix Q
//    and an upper triangular matrix R.
// 2) M*z = w  <=>  Q*R*z = w  <=>  R*z = Q^t*w.
// 3) Solve the upper triangular system  R*z = w';  w' = Q^t*w.

template<class F>
auto FMatrix<F>::solve (const V& w, Algo::How how) const -> V
{
    info (verbose, "Solve.%s [%d]\n", Algo::name[how], w.dim_);

    int s_cols = n_cols();
    if (dim_ != w.dim_ || s_cols < 1)
        fatal ("FMatrix %p.solve %p: bad [%d][%d] * [%d]", this, &w,
               dim_, s_cols, w.dim_);

    M Q, R;
    QR_decompose (Q, R, how);

    V Qt_w;
    if (how != Algo::QR_vectors)
        Qt_w = Q.transpose_mul (w);
    else
        Qt_w = Q.householder_vectors_mul (w);

    V z = R.solve_triangle (Qt_w);

    info (verbose, "Done Solve.%s [%d]\n", Algo::name[how], w.dim_);

    if (verbose)
    {
        Float::stat.scale = 0;
        V delta = (*this) * z - w;
        info (verbose, "Done Solve.%s [%d] delta = %e\n", Algo::name[how], w.dim_,
              (double) delta.abs());
        Float::stat.scale = 1;
    }
    return z;
}



// Orthogonalization on the rows.  The µ's are the optional coefficients
// as used in the projections.  Pass µ from the LLL-algorithm to avoid
// their re-computation.  µ is a square matrix with the same mumber of
// rows like *this.

template<class F>
auto FMatrix<F>::gram_schmidt (bool normalize, M *mu) const -> M
{
    M m (*this);
    m.gram_schmidt (normalize, 0, mu);
    return m;
}

template<class F>
void FMatrix<F>::gram_schmidt (bool normalize, int from_row, M *mu)
{
    M& m = *this;
    for (int r = from_row; r < dim_; ++r)
    {
#if 0
        V row { m[r] };
        for (int v = 0; v < r; ++v)
        {
            F proj = m[r].project_on (m[v]);
            row -= proj * m[v];
            if (mu)
                (*mu)[r][v] = std::move (proj);
        }

        m[r] = std::move (row);
#else
        // Operating iteratively on the current row is the *modified*
        // Gram-Schmidt which has better numerical stability than the
        // original version of the algorithm.

        for (int v = 0; v < r; ++v)
        {
            F proj = m[r].project_on (m[v]);
            m[r] -= proj * m[v];
            if (mu)
                (*mu)[r][v] = std::move (proj);
        }
#endif // 0

        if (normalize)
            m[r] = std::move (m[r].normed());
    }
}


#include <stack>

namespace {
    template<class F>
    struct FStack
    {
        std::stack<int> sp;
    };

    FStack<Float> vstack;
    //FStack<double> vstack;
}; // ::anon

template<class F>
void FMatrix<F>::push_verbose (int v)
{
    vstack.sp.push (verbose);
    verbose = v;
}

template<class F>
int FMatrix<F>::pop_verbose()
{
    if (vstack.sp.size())
    {
        verbose = vstack.sp.top();
        vstack.sp.pop();
    }
    return verbose;
}

template<class F> const F FMatrix<F>::element_0 { 0.0 };
template<class F> const F FMatrix<F>::element_1 { 1.0 };
template<> const double FMatrix<double>::element_nan { std::nan("") };
template<> const Float  FMatrix<Float>::element_nan { };


#include "mp-float.h"


int LLL::verbose = 0;

// Used during LLL-algorithm to avoid re-computations of Gram-Schmidt
// orthogonalization.  Namings follow the LLL-algorithm on Wikipedia.  */
class LLL::LLLData
{
    friend LLL;

    // Only rows up to and including .row_valid_ of .b_star_ and .mu_ are valid.
    int row_valid_ = -1;

    // .b_ is the original basis that is to be reduced.
    // .b_star_ is a matrix we get from Gram-Schmidt orthogonalization
    //      from .b_ and that is repreatedly updated as .b_ evolves.
    // .mu_ is a square matrix that contains scalars generated during that
    //      Gram-Schmidt process.  Only entries below the main diagonal are
    //      valid resp. used.
    M b_, b_star_, mu_;

    LLLData (const M& b)
        : b_(b), b_star_(Alloc (b.dim())), mu_(Alloc (b.dim(), b.dim())) {}

    void swap_b_rows (int r1, int r2)
    {
        std::swap (b_[r1], b_[r2]);
        validate_row (r1, false);
        validate_row (r2, false);
    }

    V& b_row (int r)
    {
        return b_[r];
    }

    const V& b_star_row (int r)
    {
        recompute_rows_upto (r);
        return b_star_[r];
    }

    const F& mu (int r, int c)
    {
        assert (r > c);
        recompute_rows_upto (r);
        return mu_[r][c];
    }

    void validate_row (int r, bool valid)
    {
        row_valid_ = valid
            ? std::max (row_valid_, r)
            : std::min (row_valid_, r - 1);
    }

    void recompute_rows_upto (int row)
    {
        M& m = b_star_;

        for (int r = 1 + row_valid_; r <= row; ++r)
        {
            // Operating iteratively on the current row is the *modified*
            // Gram-Schmidt which has better numerical stability than the
            // original algorithm.

            m[r] = b_[r];

            for (int v = 0; v < r; ++v)
            {
                F& mu_rv = mu_[r][v];
                mu_rv = m[r].project_on (m[v]);
                m[r].sub_scaled (mu_rv, m[v]);
            }
        }

        validate_row (row, true);
    }
};


auto LLL::reduce (const M& b, double delta) -> M
{
    if (delta <= 0.25 || delta >= 1.0)
        fatal ("delta = %f not in (1/4, 1)", delta);

    LLLData lll (b);

    for (int k = 1; k < b.dim(); )
    {
        if (verbose) out ("LLL k = %d: ", k);

        for (int j = k - 1; j >= 0; --j)
        {
            const Float& mu_kj = lll.mu (k, j);
            if (mu_kj.abscmp (0.5) > 0)
            {
                if (verbose) printf (" %d-%d", k, j);
                lll.b_row (k).sub_scaled (mu_kj.round(), lll.b_row (j));
                lll.validate_row (k, false);
            }
        }

        Float bk  = lll.b_star_row(k).abs2();
        Float bk1 = lll.b_star_row(k-1).abs2();
        const Float& mm  = lll.mu (k, k-1);

        if (bk >= (Float{delta} - mm*mm) * bk1)
            ++k;
        else
        {
            if (verbose) printf (" %d<->%d", k, k-1);
            lll.swap_b_rows (k, k-1);
            k = std::max (k-1, 1);
        }

        if (verbose) out ("\n");
    }

    return std::move (lll.b_);
}


auto LLL::find_relation (const V& v, double LLL_delta, const F& zoom) -> V
{
    const int& n = v.dim();

    M m { M::make (n, 1 + n, F{0}) };

    for (int i = 0; i < n; ++i)
    {
        m[i][i] = F{1};
        m[i][n] = zoom * v[i];
    }

    M r = LLL::reduce (m, LLL_delta);

    return r[0];
}

auto LLL::find_power_relation (const F& x, int deg, double delta, const F& zoom) -> V
{
    V v { Alloc (1 + deg) };
    F xk { 1 };

    for (int k = 0; k <= deg; ++k)
    {
        v[k] = xk;
        xk *= x;
    }

    return LLL::find_relation (v, delta, zoom);
}


// Explicit instanciate Float.

template const Float FMatrix<Float>::element_nan;
template FMatrix<Float> operator * (const Float&, const FMatrix<Float>&);
template std::ostream& operator << (std::ostream&, const FMatrix<Float>&);
template void FMatrix<Float>::push_verbose (int);
template int FMatrix<Float>::pop_verbose();
template int FMatrix<Float>::verbose;
template<> int FMatrix<Float>::verbose;

template class FMatrix<Float>;


// Explicit instanciate double.

template const double FMatrix<double>::element_nan;
template FMatrix<double> operator * (const double&, const FMatrix<double>&);
template std::ostream& operator << (std::ostream&, const FMatrix<double>&);
template void FMatrix<double>::push_verbose (int);
template int FMatrix<double>::pop_verbose();
template int FMatrix<double>::verbose;
template<> int FMatrix<double>::verbose;

template class FMatrix<double>;
