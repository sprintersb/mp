#ifndef MP_VECTOR_H
#define MP_VECTOR_H

#include <iostream>
#include <initializer_list>
#include <cassert>
#include <vector>

struct Alloc
{
    int x1, x2;
    Alloc() = delete;
    Alloc (int x1) : x1(x1), x2(-1) {};
    Alloc (int x1, int x2) : x1(x1), x2(x2) {};
};

template<class F> class FMatrix;

template<class F>
class FVector
{
    friend class FMatrix<F>;
    void alloc (int);
    void free();
private:
    int dim_ = 0;
    F* v_ = nullptr;

    bool is_index (int n) const
    {
        return n >= 0 && n < dim_;
    }

public:
    typedef F element_type;
    using V = FVector<F>;
    using M = FMatrix<F>;

    static const F element_0;
    static const F element_1;
    static const F element_nan;
    static const F element_inf;

    ~FVector();
    FVector() {}

    FVector (const FVector&);
    FVector (FVector&&);
    explicit FVector (Alloc a) { assert (a.x2 == -1); alloc (a.x1); }
    FVector (std::initializer_list<F>);
    FVector& operator= (const FVector&);
    FVector& operator= (FVector&&);

    FVector (const std::vector<F>&);
    FVector (std::vector<F>&&);
    FVector& operator= (const std::vector<F>&);
    FVector& operator= (std::vector<F>&&);
    explicit operator std::vector<F>() const;

    int dim() const { return dim_; }
    bool vector_p() const;
    static V make (int dim, const F& = element_nan);
    static V e_k (int dim, int k, const F& k_value = element_1);
    static V random_range (int dim, const F&, const F&);

    F& operator[] (int);
    const F& operator[] (int) const;

    F operator * (const FVector<F>&) const;
    void operator -= (const FVector<F>&);
    void operator += (const FVector<F>&);
    void operator *= (const F&);
    void operator /= (const F&);
    void sub_scaled (const F&, const FVector<F>&);
    FVector<F> operator - () const;
    FVector<F> operator - (const FVector<F>&) const;
    FVector<F> operator + (const FVector<F>&) const;
    FVector<F> operator * (const F&) const;
    FVector<F> operator / (const F&) const;
    F abs2 () const;
    F abs () const;
    V normed () const;
    F min () const;
    F max () const;
    F height () const;

    V reflect_to (const FVector&) const;
    M dyade () const;
    M dyade (const V&) const;
    M householder (const V&) const;
    F tail_mul (int, const V&) const;
    F tail_mul_tail (int, const V&) const;
    F mul_column (const M&, int) const;
    F mul_column_tail (int, const M&, int) const;
    F project_on (const V&) const;

    template<class S>
    friend FVector<S> operator * (const S&, const FVector<S>&);

    template<class S>
    friend std::ostream& operator << (std::ostream&, const FVector<S>&);
};

template<class S>
FVector<S> operator * (const S&, const FVector<S>&);

#include "mp-float.h"

// At the end of mp-float.cpp:

extern template class FVector<double>;

extern template FVector<double> operator * (const double&, const FVector<double>&);
//extern template std::ostream& operator << (std::ostream&, const FVector<double>&);

extern template class FVector<Float>;

extern template FVector<Float> operator * (const Float&, const FVector<Float>&);
//extern template std::ostream& operator << (std::ostream&, const FVector<Float>&);

#endif // MP_VECTOR_H
