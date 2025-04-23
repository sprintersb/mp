#ifndef MP_COMPLEX_H
#define MP_COMPLEX_H

#include "mp-float.h"

#include <iosfwd> // Just Forward decls for iostream.

class Int;
class Ratio;
template<class A> class Poly;

class Complex
{
    Float x_, y_;

public:
    static bool out_with_precision;
    static const char *out_format;
    static int Debug;

    ~Complex() = default;
    Complex() = default;
    Complex (double);
    Complex (double, double);
    Complex (Float&&);
    Complex (const Float&);
    Complex (const Float&, const Float&);

    Complex (const Complex&) = default;
    Complex (Complex&&) = default;
    Complex& operator= (const Complex&) = default;
    Complex& operator= (Complex&&) = default;
    //Complex& operator= (Float&&);

    //explicit operator double() const;

    Complex& set_precision_round (int bits, mpfr_rnd_t);
    mpfr_prec_t get_precision () const;
    static void SetPrecision (int digits, int radix);
    static int GetPrecision (int radix);
    static mpfr_prec_t SaturatePrecision (int);

    Complex conj () const; // Conjugation
    Complex operator ~ () const; // Conjugation
    Complex operator - () const;
    Complex operator + (const Complex&) const;
    Complex operator - (const Complex&) const;
    Complex operator * (const Complex&) const;
    Complex operator / (const Complex&) const;
    Complex inv () const; // 1 / Complex
    void operator -= (const Float&);
    void operator += (const Float&);
    void operator *= (const Float&);
    void operator /= (const Float&);
    void operator += (const Complex&);
    void operator -= (const Complex&);
    void operator *= (const Complex&);
    void operator /= (const Complex&);
    Complex operator << (int) const;
    Complex operator >> (int) const;

    bool operator == (const Complex&) const;
    bool operator != (const Complex&) const;
    bool operator == (const Float&) const;
    bool operator != (const Float&) const;

    bool nan_p () const;
    bool inf_p () const;
    bool number_p () const;
    bool zero_p () const;
    bool nonzero_p () const;

    static const Complex One;
    static const Complex Zero;
    bool is_zero() const { return zero_p(); }
    bool is_integer() const { return x_.is_integer() && y_.zero_p(); }

    Float abs () const;
    Float norm () const;
    Float abs2 () const { return norm(); }
    Float absdiff (const Complex& z) const;
    Float arg() const;
    Float& real();
    Float& imag();
    const Float& real() const;
    const Float& imag() const;

    Complex round () const;
    Complex trunc () const;
    Complex sqrt () const;
    Complex exp () const;
    Complex log () const;
    Complex sin () const;
    Complex cos () const;
    Complex tan () const;
    Complex asin () const;
    Complex acos () const;
    Complex atan () const;
    Complex sinh () const;
    Complex cosh () const;
    Complex tanh () const;
    Complex asinh () const;
    Complex acosh () const;
    Complex atanh () const;

    Complex ldexp (int) const;
    Complex pow (const Complex&) const;
    Complex root (unsigned long) const;
    int abscmp (const Complex&) const;

    Complex powi (long) const;

    Complex linear (const Complex&, const Complex&) const;

    static Complex from_str (const char*, int base = 0);
    char* to_str (const char* = nullptr) const; // malloc.
    void print (const char*) const;
    std::ostream& print (std::ostream&, const char* = nullptr) const;
    friend std::ostream& operator << (std::ostream&, const Complex&);
    void dump() const;

    Poly<Complex> operator + (const Poly<Complex> &q) const;
    Poly<Complex> operator - (const Poly<Complex> &q) const;
    Poly<Complex> operator * (const Poly<Complex> &q) const;

    static Complex make (int);
};

Complex operator + (const Float&, const Complex&);
Complex operator - (const Float&, const Complex&);
Complex operator * (const Float&, const Complex&);
Complex operator / (const Float&, const Complex&);

Complex operator "" _C (const char*);
Complex operator "" _C (const char[], size_t);

#include <type_traits>

namespace PolyHelp
{
#ifndef HAVE_POLYHELP
#define HAVE_POLYHELP
    template<typename T> struct real_type {};
    template<typename T> struct Prio {};
#endif

    template<> struct real_type<Complex> { using type = Float; };
    
    template<>
    struct Prio<Complex> { using type = std::integral_constant<int,10>; };
}

#endif // MP_COMPLEX_H
