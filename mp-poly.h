#ifndef MP_POLY_H
#define MP_POLY_H

class Int;
class Ratio;
class Float;
class Complex;
template<typename A> class Poly;
template<typename A> class RationalFunction;

#include <vector>
#include <cmath>
#include <iosfwd> // Just Forward decls for iostream.
#include <initializer_list>
#include <type_traits>

#include <stdio.h>
#include <complex>

namespace PolyHelp
{
#ifndef HAVE_POLYHELP
#define HAVE_POLYHELP
    template<typename T> struct real_type {};
    template<typename T> struct Prio {};
#endif

    template<> struct real_type<std::complex<double>> { using type = double; };

    template<> struct real_type<double> { using type = double; };
    template<> struct real_type<Float>  { using type = Float; };
    template<> struct real_type<Ratio>  { using type = Ratio; };
    template<> struct real_type<Int>    { using type = Int; };

#define _PPRIO(T, N)                                                    \
    template<>                                                          \
    struct Prio<T> { using type = std::integral_constant<int,N>; }

    _PPRIO (Int, 2);
    _PPRIO (Ratio, 3);
    _PPRIO (double, 4);
    _PPRIO (std::complex<double>, 5);
    _PPRIO (Float, 6);

    // The result type when we evaluate Poly<A> at X, used below
    // for Poly<A>::result_type<X>.
    template<typename A, typename X>
    struct Result
    {
        static constexpr bool cond =
            Prio<A>::type::value > Prio<X>::type::value;
        using type = typename std::conditional<cond, A, X>::type;
    };

    template<typename A, typename X>
    struct Result<A, Poly<X>>
    {
        using type = Poly<typename Result<A,X>::type>;
    };

    template<typename A, typename X>
    struct Result<A, RationalFunction<X>>
    {
        using type = RationalFunction<typename Result<A,X>::type>;
    };

} // PolyHelp

template<class A>
class Poly
{
    typedef std::vector<A> As;
public:
    using scalar_type = A;
    using scalar_real_type = typename PolyHelp::real_type<A>::type;
    struct Scalar
    {
        using type = A;
        using real_type = typename PolyHelp::real_type<A>::type;
        static A get (int);
        static A zero() { return get(0); };
        static A one() { return get(1); }
    };
private:
    int deg_ = 0;
    As a_;

public:
    static int Debug;

    template<typename X>
    using result_type = typename PolyHelp::Result<A,X>::type;

    Poly () : deg_(0), a_(As(1, Scalar::get(0))) {}
    Poly (const As&);
    Poly (As&&);
    Poly (std::initializer_list<A>);
    int deg() const { return deg_; }
    int set_deg ();
    A& at (int);
    const A& at (int) const;
    A& operator [] (int i) { return at(i); }
    const A& operator [] (int i) const { return at(i); }
    A operator () (const A&) const;
    Poly operator () (const Poly&) const;
    const As& a() const { return a_; }
    std::vector<A> as_vector() const { return a_; }
    // Will adjust degree if i > deg.
    void set (int i, const A&);
    void clear (int);
    template<class X>
    result_type<X> operator () (const X&) const;
    Poly operator - () const;
    Poly operator - (const Poly&) const;
    Poly operator + (const Poly&) const;
    Poly operator * (const Poly&) const;
    Poly& operator -= (const Poly& q) { return *this = (*this) - q; }
    Poly& operator += (const Poly& q) { return *this = (*this) + q; }
    Poly& operator *= (const Poly& q) { return *this = (*this) * q; }
    Poly& operator -= (const A& a) { at(0) -= a; return *this; }
    Poly& operator += (const A& a) { at(0) += a; return *this; }
    Poly& operator *= (const A&);
    Poly& operator /= (const A&);
    Poly operator - (const A& a) const { Poly q (*this); return q -= a; }
    Poly operator + (const A& a) const { Poly q (*this); return q += a; }
    Poly operator * (const A& a) const { Poly q (*this); return q *= a; }
    Poly operator / (const A& a) const { Poly q (*this); return q /= a; }

    Poly operator / (const Poly&) const;
    Poly operator % (const Poly&) const;
    bool operator == (const Poly&) const;
    bool operator != (const Poly& q) const { return ! (*this == q); }
    bool zero_p() const;
    bool polynomial_p() const;
    int ord_x() const;
    Poly pow (const Int&) const;
    Poly operator << (int) const;
    Poly operator >> (int) const;
    Poly& operator <<= (int);
    Poly& operator >>= (int);
    Poly D() const;
    Poly swap() const; // Swap highest and lowest: return x^deg(p) * p(1/x)
    scalar_real_type height() const;
    Poly crop (int lo, int hi, int stride = 1) const;
private:
    Poly msub (const Poly&, const A&, int) const;
public:
    Poly divmod (const Poly&, Poly& rem) const;

    template<class Y> operator Poly<Y>() const
    {
        Poly<Y> y;
        for (int i = deg_; i >= 0; --i)
            y.set (i, (Y) at(i));

        return y;
    }

    operator RationalFunction<A>() const
    {
        return RationalFunction<A> { *this, Poly { Scalar::get(1) } };
    }

    RationalFunction<A> operator - (const RationalFunction<A>&) const;
    RationalFunction<A> operator + (const RationalFunction<A>&) const;
    RationalFunction<A> operator * (const RationalFunction<A>&) const;
    RationalFunction<A> operator / (const RationalFunction<A>&) const;

    enum Style
    {
        None, Gnuplot, GnuplotDown, GnuplotHorner,
        HTML, HTMLDown, HTMLMinus, HTMLMinusDown, // Minus: --> "&minus;"
        TeX, TeXDown, TeXFrac, TeXFracDown, Text, TextDown,
        Python, PythonHorner, PythonHornerPow, CHorner, Desmos, DesmosDown,
        VHDLTab, FLT40tab,
        // Only lists following.  The string argument of print() is printed
        // as separator between the list elements.  If the string looks like
        // "L@R|S", then
        // L is printed left of either element,
        // R is printed right of either element, and
        // S is printed as seperator between them.
        // For example, to print F7_CONST_DEF like in libf7-const.def,
        // you can use "F7_CONST_DEF (X, @)\n|".
        List, ListDouble, ListDoubleDown,
        ListFloat, ListFloatHex,
        // Print flags, bytes and expo like "0, 0x80,0x00,... 0" for 1.
        // For ListF7Normalized and if the highest coefficient is normalized,
        // then skip the highest coefficient and set a special flag (8) for
        // the 2nd-highest coefficient.
        ListF7, ListF7Normalized,
        ListDoubleX64, ListFloatX32,
        ListVHDLFloatX32
    };
    // If VAR is like "L@R|S", then:
    // * L is printed left of either coefficient.
    // * R is printed right of either coefficient.
    // * If STYLE is a list: S is the separator printed between coefficients.
    // * If STYLE is a non-list: S is the polynomials's variable name.
    // VAR defaults to "F7_CONST_DEF (X, @)\n|" if style is ListF7[Normalized]
    // and "x", otherwise.  Default list separator is ", ".
    std::ostream& print (std::ostream&, Style = None, const char *var = nullptr,
                         int prec2 = 0) const;

    template<class X> friend Poly<X> poly_sin (int);
    template<class X> friend Poly<X> poly_asin (int);

    void print_VHDL_Table (std::ostream&, const char* var, const char* typ, int prec2 = 0, int n_coeff = -1) const;
private:
    std::ostream& print_horner (std::ostream&, Style, const char*, int, const char*, const char*) const;
    std::ostream& print_VHDLtab (std::ostream&, const char*, int prec2 = 0) const;
    std::ostream& print_FLT40tab (std::ostream&, const char* = nullptr) const;
};

// Just for a test Eval Poly<Float> at RatioalFunction<Float>
extern template auto Poly<Float>::operator () (const RationalFunction<Float>&) const -> Poly<Float>::result_type<RationalFunction<Float>>;


// For different code, use specialization
// template<> template<> Poly<XXX>::operator Poly<YYY>() const;


extern template Poly<Ratio>::operator Poly<Float>() const;
extern template Poly<Ratio>::operator Poly<double>() const;
extern template Poly<Float>::operator Poly<double>() const;
extern template Poly<double>::operator Poly<Float>() const;
extern template Poly<Float>::operator Poly<Ratio>() const;

extern template std::complex<double> Poly<double>::operator () (const std::complex<double>&) const;

template<class A> Poly<A> poly_sin (int);
template<class A> Poly<A> poly_cos (int);
template<class A> Poly<A> poly_cossqrt (int);
template<class A> Poly<A> poly_sin_x (int);
template<class A> Poly<A> poly_sinsqrt_sqrt (int);
template<class A> Poly<A> poly_asin (int);

template<class A> std::ostream& operator << (std::ostream&, const Poly<A>&);

extern template class Poly<Ratio>;
extern template class Poly<Float>;
extern template class Poly<double>;

template<class A>
class RationalFunction
{
public:
    typedef A scalar_type;
    typedef Poly<A> polynomial_type;
    typedef RationalFunction Self;
    Poly<A> p_;
    Poly<A> q_;
    A operator () (const A& x) const { return p_(x) / q_(x); };
    bool is_polynomial () const;
    Self normalize_q (int) const;
    Self normalize_p (int) const;
    Self normalize_q0() const { return normalize_q (0); };
    Self normalize_p0() const { return normalize_p (0); };
    Self normalize_q_deg() const { return normalize_q (q_.deg()); };
    Self normalize_p_deg() const { return normalize_p (p_.deg()); };
    Self at_1_over_x() const; // return Self (1/x)
    Self operator - () const;
    Self operator - (const Self&) const;
    Self operator + (const Self&) const;
    Self operator * (const Self&) const;
    Self operator / (const Self&) const;
    Self operator - (const Poly<A>&) const;
    Self operator + (const Poly<A>&) const;
    Self operator * (const Poly<A>&) const;
    Self operator / (const Poly<A>&) const;

    Self operator () (const Self&) const;
    Self operator () (const Poly<A>&) const;

    void print_VHDL_Tables (std::ostream&, const char* varP, const char *varQ,
                            const char *typ, int prec2 = 0) const;
};

extern template class RationalFunction<Float>;
extern template class RationalFunction<Ratio>;

template<class A>
Poly<A> poly_sin (int d)
{
    Poly<A> f;
    A ai { 1 };
    f.a_.resize (1 + d, A{0});
    for (int i = 0; i <= d; i++)
    {
        if (i % 2 != 0)
        {
            f.a_[i] = ai;
            ai /= A (- (i + 1) * (i + 2));
        }
    }
    f.set_deg();
    return f;
}

template<class A>
Poly<A> poly_cos (int d)
{
    return poly_sin<A> (1 + d).D();
}

template<class A>
Poly<A> poly_cossqrt (int d)
{
    return poly_cos<A> (2*d).crop (0, 2*d, 2);
}

template<class A>
Poly<A> poly_sin_x (int d)
{
    return poly_sin<A> (1 + d) >> 1;
}

template<class A>
Poly<A> poly_sinsqrt_sqrt (int d)
{
    return poly_sin<A> (2*d+1).crop (1, 2*d+1, 2);
}

template<class A>
Poly<A> poly_asin (int d)
{
    Poly<A> f;
    A ai { 1 };
    f.a_.resize (1 + d, A{0});
    for (int i = 0; i <= d; i++)
    {
        if (i % 2 != 0)
        {
            f.a_[i] = ai / A(i);
            ai *= A (i) / A (i + 1);
        }
    }
    f.set_deg();
    return f;
}

template<class A>
Poly<A> poly_asinsqrt_sqrt (int d)
{
    return poly_asin<A> (2*d+1).crop (1, 2*d+1, 2);
}

// asin(w)/w, w = sqrt(x/2)
template<class A>
Poly<A> poly_asinsqrt2_sqrt2 (int d)
{
    Poly<A> x2 { A(0), A(1) / A(2) };
    return poly_asinsqrt_sqrt<A> (d) (x2);
}

#endif // MP_POLY_H
