#ifndef MP_MATRIX_H
#define MP_MATRIX_H

#include <iostream>
#include <initializer_list>
#include <cassert>

#include "mp-vector.h"

/*
 * Some mimimalist functions for vector and matrix for scalars F.
 */

class Algo
{
public:
    enum How { QR_slow, QR_fast, QR_vectors };
    static const char *const name[];
};


// A Matrix of FVector's as rows.
template<class F>
class FMatrix
{
public:
    typedef FVector<F> row_type;
    //typedef FVector<F> V;
    //using M = typedef FMatrix<F> M;
    using V = FVector<F>;
    using M = FMatrix<F>;
    typedef F element_type;

    static const F element_0;
    static const F element_1;
    static const F element_nan;

    void free();
    void alloc (int);
    void alloc (int, int);

    static int verbose;
    static void push_verbose (int);
    static int pop_verbose();

private:
    int dim_ = 0;
    row_type *r_ = nullptr;

    bool is_index (int n) const
    {
        return n >= 0 && n < dim_;
    }

public:

    ~FMatrix();
    FMatrix() {}

    FMatrix (const M&);
    FMatrix (M&&);
    FMatrix (Alloc a) { alloc (a.x1, a.x2); }
    FMatrix (std::initializer_list<V>);
    FMatrix& operator = (const M&);
    FMatrix& operator = (M&&);

    int dim() const;
    int n_cols() const; // < 0 if row[-r] has other dim than row[0].
    bool matrix_p() const;

    static M make (int dim, int cols = 0, const F& = element_nan);
    static M make (int dim, const V&);
    static M id (int, const F& = FMatrix<F>::element_1);
    static M random_range (int dim, int cols, const F&, const F&);

    V& operator[] (int row);
    const V& operator[] (int row) const;

    FVector<F> column (int, const char *caller = "") const;
    FMatrix<F> T() const;
    void operator -= (const FMatrix<F>&);
    void operator += (const FMatrix<F>&);
    void operator *= (const FMatrix<F>&);
    void operator *= (const F&);
    void operator /= (const F&);
    FMatrix<F> operator - () const;
    FMatrix<F> operator - (const FMatrix<F>&) const;
    FMatrix<F> operator + (const FMatrix<F>&) const;
    FMatrix<F> operator * (const FMatrix<F>&) const;
    FVector<F> operator * (const FVector<F>&) const;
    FMatrix<F> operator * (const F&) const;
    FMatrix<F> operator / (const F&) const;

    void apply (void (*)(M&, int, int, void*), void*);
    void apply (void (*)(const M&, int, int, void*), void*) const;

    void QR_decompose (M& Q, M& R, Algo::How) const;
    V reflect_to_e_k (int) const;
    M householder (int) const;
    void mul_minor_right (int, const M&);
    void mul_minor_left (int, const M&);
    void mul_householder_minor_left (int, const V&);
    void mul_householder_minor_right (int, const V&);
    void QR_householder_step (int, M&, M&) const;
    void QR_householder_step_fast (int, M&, M&) const;
    void QR_householder_step_reflection_only (int, M&, M&) const;

    V solve (const V&, Algo::How = Algo::QR_vectors) const;
    V solve_triangle (const V&) const;
    V transpose_mul (const V&) const;
    V householder_vectors_mul (const V&) const;

    // Orthogonalization on the rows.
    void gram_schmidt (bool normalize, int from_row, M* mu = nullptr);
    M gram_schmidt (bool normalize, M* mu = nullptr) const;

    template<class S>
    friend FMatrix<S> operator * (const S&, const FMatrix<S>&);

    template<class S>
    friend std::ostream& operator << (std::ostream&, const FMatrix<S>&);
};

template<class S>
FMatrix<S> operator * (const S&, const FMatrix<S>&);


#include "mp-float.h"

// At the end of mp-float.cpp:

extern template class FMatrix<double>;

extern template FMatrix<double> operator * (const double&, const FMatrix<double>&);
//extern template std::ostream& operator << (std::ostream&, const FMatrix<double>&);

extern template class FMatrix<Float>;

extern template FMatrix<Float> operator * (const Float&, const FMatrix<Float>&);
//extern template std::ostream& operator << (std::ostream&, const FMatrix<Float>&);

class LLL
{
    class LLLData;

    using F = Float;
    using V = FVector<F>;
    using M = FMatrix<F>;

    LLL() = delete;

public:
    static int verbose;

    // The very LLL algorithm.  DELTA in (1/4, 1) where the reduction is
    // better if DELTA is greater; but that also slows down the algorithm.
    // A typical choice is DELTA = 0.75.
    static M reduce (const M&, double delta);

    // Find integers a[] such that the relation  a[] * v[] = 0  holds
    // approximately.  DELTA is passed to the LLL-algorithm.  ZOOM is a scaling
    // factor for v[].  The bigger ZOOM, the more likely it is to find the
    // relation (providing it exists), but the algorithm will be slower and
    // it requires for higher Float precision.
    //
    // The returned vector has one excess element:  If the last element
    // is (close to) zero, then there is a relation.  Otherwise, there is no
    // relation, or it could not be found with the provided precision and
    // values of DELTA and ZOOM.
    //
    // A = max |a[i]|, log = log_10 and be P the precision in decimal digits.
    // Then:
    // * The precision should be big enough so it can descrimitante
    //   between  A*ZOOM  and  1.0 + A*ZOOM, i.e.  P > log A + log ZOOM.
    // * ZOOM should satisfy  ZOOM > A^dim.
    static V find_relation (const V& v, double delta, const F& zoom);

    // Find a relation amongst 1, x, x^2, ... x^deg using LLL::find_relation.
    static V find_power_relation (const F&, int deg, double delta, const F& zoom);
};

#endif // MP_MATRIX_H
