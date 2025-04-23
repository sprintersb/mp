#define MP_POLY_CPP

#include "diagnostic.h"

#include "mp-poly.h"
#include "mp-int.h"
#include "mp-ratio.h"
#include "mp-float.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cinttypes>
#include <ostream>

namespace
{
    using _complex = std::complex<double>;

    template<class A> bool _is_zero (const A& a) { return a == 0; }
    template<> bool _is_zero<Float> (const Float& a) { return a.is_zero(); }
    template<> bool _is_zero<Ratio> (const Ratio& a) { return a.is_zero(); }

    template<class A,class I> A _make (I i);// { return A { i }; }
    template<class A> A _make (int i);// { return A { i }; }

    template<> double _make (int i) { return (double) i; }
    template<> Float _make (int i) { return Float::from_int (i); }
    template<> Ratio _make (int i) { return i ? Ratio { i } : 0_Ratio; }
    template<> Float _make (Ratio i) { return (Float) i; }
    template<> _complex _make (double i) { return (_complex) i; }
    template<> Poly<Float> _make (Ratio i) { return Poly<Float> { (Float) i }; }

    template<> RationalFunction<Float>
    _make (Float f)
    {
        return RationalFunction<Float> { Poly<Float>{f}, Poly<Float>{1} };
    }
    template<> RationalFunction<Ratio>
    _make (Ratio f)
    {
        return RationalFunction<Ratio> { Poly<Float>{f}, Poly<Float>{1} };
    }

    template<class A> bool _isnan (const A&) { return false; }
    template<> bool _isnan (const Float& a) { return a.nan_p(); }
    template<> bool _isnan (const double& a) { return std::isnan (a); }

    template<class A> bool _isinf (const A&) { return false; }
    template<> bool _isinf (const Float& a) { return a.inf_p(); }
    template<> bool _isinf (const double& a) { return std::isinf (a); }

    template<class A> bool _isnumber (const A &a) { return !_isnan<A>(a) && !_isinf<A>(a); }

    template<class A> A _max (const A& a, const A& b) { return std::max (a, b); }
    template<class R, class A> R _abs (const A& a) { return abs (a); }

    template<class X> X _fma (const X& x, const X& y, const X& z)
    {
        return x * y + z;
    }
    template<> Float _fma (const Float& x, const Float& y, const Float& z)
    {
        return x.fma (y, z);
    }
    template<> double _fma (const double& x, const double& y, const double& z)
    {
        return std::fma (x, y, z);
    }
}; // ::anon


#include "mp-poly-impl.h"


// Explicit instanciate some non-trivial evaluations.
#define _poly_eval(A, X)                                                \
    template auto Poly<A>::operator () (const X&) const -> result_type<X>

_poly_eval (double, std::complex<double>);
_poly_eval (Float, RationalFunction<Float>);
_poly_eval (Ratio, RationalFunction<Ratio>);
_poly_eval (Ratio, Poly<Float>);

template Poly<Ratio>::operator Poly<Float>() const;
template Poly<Ratio>::operator Poly<double>() const;
template Poly<Float>::operator Poly<double>() const;
template Poly<double>::operator Poly<Float>() const;

template Poly<Float>::operator Poly<Ratio>() const;

#define MP_POLY_PRINT_H
#include "mp-poly-print.h"
#undef MP_POLY_PRINT_H

Poly<Ratio> Ratio::operator + (const Poly<Ratio>& q) const { return q + (*this); }
Poly<Ratio> Ratio::operator - (const Poly<Ratio>& q) const { return -q + (*this); }
Poly<Ratio> Ratio::operator * (const Poly<Ratio>& q) const { return q * (*this); }

Poly<Float> Float::operator + (const Poly<Float>& q) const { return q + (*this); }
Poly<Float> Float::operator - (const Poly<Float>& q) const { return -q + (*this); }
Poly<Float> Float::operator * (const Poly<Float>& q) const { return q * (*this); }


template Float Poly<Ratio>::operator () (const Float&) const;
//template Poly<double> Poly<Ratio>::operator () (const Poly<double>&) const;

template std::ostream& operator << (std::ostream&, const Poly<Ratio>&);
template std::ostream& operator << (std::ostream&, const Poly<Float>&);
template std::ostream& operator << (std::ostream&, const Poly<double>&);

//template<class A> const A Poly<A>::A_One { 1 };
//template<class A> const A Poly<A>::A_Zero { 0 };
template<class A> int Poly<A>::Debug { 0 };

template class Poly<Ratio>;
template class Poly<Float>;
template class Poly<double>;

template class RationalFunction<Ratio>;
template class RationalFunction<Float>;

#undef MP_POLY_CPP
