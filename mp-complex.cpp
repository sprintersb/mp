#include "mp-complex.h"
#include "diagnostic.h"

#include <utility> // std::move
#include <ostream>
#include <cstdlib>
#include <type_traits> // std::is_same

#include <mpfr.h>
#include <mpc.h>

namespace
{
    // For subtlety of () in decltype introduced by macro mpc_realref(), see
    //    Significance of parentheses in decltype((c)) ?
    //    https://stackoverflow.com/q/14115744/1556746
    mpc_t c;
    static_assert (sizeof (mpc_t) == 2 * sizeof (mpfr_t), "");
    static_assert (std::is_same<decltype(&mpc_realref(c)), mpfr_t*>::value, "");
    static_assert (std::is_same<decltype(&mpc_imagref(c)), mpfr_t*>::value, "");
}


namespace
{
    void _mimic (mpc_ptr c, const Complex& z)
    {
        __builtin_memcpy (mpc_realref (c), z.real().mpfr(), sizeof (mpfr_t));
        __builtin_memcpy (mpc_imagref (c), z.imag().mpfr(), sizeof (mpfr_t));
    }

    Complex& _mimic (Complex& z, mpc_srcptr c)
    {
        __builtin_memcpy (z.real().mpfr(), mpc_realref (c), sizeof (mpfr_t));
        __builtin_memcpy (z.imag().mpfr(), mpc_imagref (c), sizeof (mpfr_t));
        return z;
    }
};

#define MIMIC(_mpc, _cplx) \
    mpc_t _mpc; _mimic (_mpc, _cplx)

Complex::Complex (Float&& x) : x_(std::move(x)), y_(0.0) {};

Complex::Complex (const Float& x) : x_(x), y_(0.0) {};
Complex::Complex (const Float& x, const Float& y) : x_(x), y_(y) {};

Complex::Complex (double x) : x_(x), y_(0.0) {};
Complex::Complex (double x, double y) : x_(x), y_(y) {};

Float& Complex::real() { return x_; }
Float& Complex::imag() { return y_; }
const Float& Complex::real() const { return x_; }
const Float& Complex::imag() const { return y_; }


//////////////////////////////////////////////////////////////////

void Complex::operator -= (const Float& f) { x_ -= f; }
void Complex::operator += (const Float& f) { x_ += f; }
void Complex::operator *= (const Float& f) { x_ *= f; y_ *= f; }
void Complex::operator /= (const Float& f) { x_ /= f; y_ /= f; }

void Complex::operator -= (const Complex& z)
{
    x_ -= z.x_;
    y_ -= z.y_;
}

void Complex::operator += (const Complex& z)
{
    x_ += z.x_;
    y_ += z.y_;
}

void Complex::operator *= (const Complex& z)
{
    *this = (*this) * z;
}

void Complex::operator /= (const Complex& z)
{
    *this = (*this) / z;
}

//////////////////////////////////////////////////////////////////
// (Complex, Complex) -> Complex

#define MK_OP(NAME, MPC_NAME)                       \
    Complex Complex::NAME (const Complex& z) const  \
    {                                               \
        Complex w;                                  \
        MIMIC (c, *this);                           \
        MIMIC (cz, z);                              \
        MIMIC (cw, w);                              \
        mpc_## MPC_NAME (cw, c, cz, MPFR_RNDN);     \
        return _mimic (w, cw);                      \
    }

MK_OP (operator -, sub)
MK_OP (operator +, add)
MK_OP (operator *, mul)
MK_OP (operator /, div)
MK_OP (pow, pow)

#undef MK_OP


//////////////////////////////////////////////////////////////////
// (Float, Complex) -> Complex

Complex operator + (const Float& f, const Complex& z)
{
    return Complex { f + z.real(), z.imag() };
}

Complex operator - (const Float& f, const Complex& z)
{
    return Complex { f - z.real(), -z.imag() };
}

Complex operator * (const Float& f, const Complex& z)
{
    return Complex { f * z.real(), f * z.imag() };
}

Complex operator / (const Float& f, const Complex& z)
{
    return f * z.inv();
}


//////////////////////////////////////////////////////////////////

bool Complex::zero_p () const
{
    return x_.zero_p() && y_.zero_p();
}

bool Complex::nan_p () const
{
    return x_.nan_p() || y_.nan_p();
}

bool Complex::operator == (const Complex& w) const
{
    return x_ == w.x_ && y_ == w.y_;
}

bool Complex::operator != (const Complex& w) const
{
    return x_ != w.x_ || y_ != w.y_;
}

bool Complex::operator == (const Float& f) const
{
    return x_ == f && y_.zero_p();
}

bool Complex::operator != (const Float& f) const
{
    return x_ != f || ! y_.zero_p();
}

bool Complex::inf_p() const
{
    return ! nan_p() && (x_.inf_p() || y_.inf_p());
}

bool Complex::number_p() const
{
    return ! nan_p() && ! inf_p();
}

int Complex::abscmp (const Complex& z) const
{
    MIMIC (c, *this);
    MIMIC (cz, z);
    return mpc_cmp_abs (c, cz);
}

//////////////////////////////////////////////////////////////////
// (Complex, Complex) -> Float

Float Complex::absdiff (const Complex& z) const
{
    return (*this - z).abs();
}

//////////////////////////////////////////////////////////////////
// Complex -> Float

#define MK_OP(NAME, MPC_NAME)                      \
    Float Complex::NAME() const                    \
    {                                              \
        Float f;                                   \
        MIMIC (c, *this);                          \
        mpc_## MPC_NAME (f.mpfr(), c, MPFR_RNDN);  \
        return f;                                  \
    }

MK_OP (arg, arg)
MK_OP (abs, abs)
MK_OP (norm, norm)

#undef MK_OP


//////////////////////////////////////////////////////////////////
// Complex -> Complex

Complex Complex::inv() const
{
    Float a = abs2();
    return Complex { x_ / a, -y_ / a };
}

#define MK_OP(NAME, MPC_NAME)                   \
    Complex Complex::NAME() const               \
    {                                           \
        Complex w;                              \
        MIMIC (c, *this);                       \
        MIMIC (cw, w);                          \
        mpc_## MPC_NAME (cw, c, MPFR_RNDN);     \
        return _mimic (w, cw);                  \
    }

MK_OP (operator -, neg)
MK_OP (operator ~, conj)
MK_OP (conj, conj)

MK_OP (exp, exp)
MK_OP (log, log)
MK_OP (sqrt, sqrt)

MK_OP (sin, sin)
MK_OP (cos, cos)
MK_OP (tan, tan)
MK_OP (asin, asin)
MK_OP (acos, acos)
MK_OP (atan, atan)

MK_OP (sinh, sinh)
MK_OP (cosh, cosh)
MK_OP (tanh, tanh)
MK_OP (asinh, asinh)
MK_OP (acosh, acosh)
MK_OP (atanh, atanh)
#undef MK_OK

Complex Complex::ldexp (int ex) const
{
    return Complex { x_.ldexp (ex), y_.ldexp (ex) };
}

Complex Complex::operator << (int ex) const
{
    return ldexp (ex);
}

Complex Complex::operator >> (int ex) const
{
    return ldexp (-ex);
}

Complex Complex::powi (long n) const
{
    Complex w;
    MIMIC (c, *this);
    MIMIC (cw, w);
    mpc_pow_si (cw, c, n, MPFR_RNDN);
    return _mimic (w, cw);
}

//////////////////////////////////////////////////////////////////

bool Complex::out_with_precision = true;
const char* Complex::out_format = "%RNe";

Complex Complex::from_str (const char*, int)
{
    fatal ("todo");
    return Complex();
}

Complex operator "" _C (const char *str, size_t)
{
    return Complex::from_str (str);
}

Complex operator "" _C (const char *str)
{
    return Complex::from_str (str);
}

char* Complex::to_str (const char *fmt) const
{
    char *re, *im, *str;
    mpfr_asprintf (&re, fmt ? fmt : Complex::out_format, x_.mpfr(), x_.mpfr());
    mpfr_asprintf (&im, fmt ? fmt : Complex::out_format, y_.mpfr(), y_.mpfr());
    mpfr_asprintf (&str, "(%s %s)", re, im);
    free (re);
    free (im);
    return str;
}

void Complex::print (const char* fmt) const
{
    mpfr_printf ("(");
    mpfr_printf (fmt ? fmt : Complex::out_format, x_.mpfr());
    mpfr_printf (" ");
    mpfr_printf (fmt ? fmt : Complex::out_format, y_.mpfr());
    mpfr_printf (")");
}

std::ostream& Complex::print (std::ostream& ost, const char* fmt) const
{
    bool with_prec = !fmt && Complex::out_with_precision;
    if (with_prec)
        ost << "[" << x_.get_precision() << "," << y_.get_precision() << "]=";

    char *str = to_str();
    ost << str;
    std::free (str);

    return ost;
}

std::ostream& operator << (std::ostream &ost, const Complex& z)
{
    return z.print (ost);
}

////////////////////////////////////////////////////////////////////////////////

#include "mp-poly.h"

namespace
{
    template<class A> bool _is_zero (const A& a) { return a == 0; }
    template<> bool _is_zero<Complex> (const Complex& a) { return a.is_zero(); }

    template<class A,class I> A _make (I i);// { return A { i }; }
    template<class A> A _make (int i);// { return A { i }; }

    template<> Complex _make (int i) { return Complex { Float::from_int(i) }; }
    template<> Float   _make (int i) { return Float::from_int (i); }

    template<class A> bool _isnan (const A&) { return false; }
    template<> bool _isnan (const Complex& a) { return a.nan_p(); }

    template<class A> A _max (const A& a, const A& b);
    template<> Float _max (const Float& x, const Float &y) { return max (x,y); }

    template<class R, class A> R _abs (const A& a) { return abs (a); }
    template<> Float _abs (const Complex& z) { return z.abs(); }
        
    template<class X> X _fma (const X& x, const X& y, const X& z)
    {
        return x * y + z;
    }

    template<> Complex _fma (const Complex& x, const Complex& y,
                             const Complex& z)
    {
        Complex w;
        MIMIC (cx, x);
        MIMIC (cy, y);
        MIMIC (cz, z);
        MIMIC (cw, w);
        mpc_fma (cw, cx, cy, cz, MPFR_RNDN);
        return _mimic (w, cw);
    }
    template<class A> bool _isnumber (const A &a)
    {
        return ! a.nan_p() && ! a.inf_p();
    }
} // ::anon


#include "mp-poly-impl.h"

namespace
{
/*    template<class A> bool _is_zero (const A& a) { return a == 0; }
    template<> bool _is_zero<Float> (const Float& a) { return a.is_zero(); }
    template<> bool _is_zero<Ratio> (const Ratio& a) { return a.is_zero(); }
*/
    template<class A,class I> A _make (I i);// { return A { i }; }
    template<class A> A _make (int i);// { return A { i }; }

    template<> Complex _make (Float f) { return Complex { f }; }
}

// Explicit instanciate some non-trivial evaluations.

extern template auto Poly<Float>::operator () (const Complex&) const -> Complex;
template auto Poly<Float>::operator () (const Complex&) const -> Complex;


template<>
std::ostream& Poly<Complex>::print (std::ostream& ost, Style,
                                    const char*, int /*prec2*/) const
{
    ost << "todo:Poly<Complex>";
    return ost;
}

/*template<> Complex Poly<Complex>::height() const
{
    error ("%s", "todo: Poly<Complex>::height()");
    }*/

template std::ostream& operator << (std::ostream&, const Poly<Complex>&);

//template<class A> const A Poly<A>::A_One { 1 };
//template<class A> const A Poly<A>::A_Zero { 0 };
//template<class A> int Poly<A>::Debug { 0 };

template class Poly<Complex>;
