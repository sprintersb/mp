#ifndef MP_FLOAT_H
#define MP_FLOAT_H

/*
 * Purpose of this header module is to introduce features that allow to
 * use double and Float (with C++ wraps for mpfr_t) by means of the same
 * interfaces without the need to wrap a class around double.
 *
 * There are some new double functions and variables though, so that the
 * same interface can be maintained, e.g. to output these types.
 */

/* From gmp.h:

   The prototypes for gmp_vprintf etc are provided only if va_list is defined,
   via an application having included <stdarg.h>.  Usually va_list is a typedef
   so can't be tested directly, but C99 specifies that va_start is a macro.

   A similar mechanism is used in <mpfr.h>.  Thus, include <cstdarg> now. */

#include <cstdarg>

#include <mpfr.h>

#include <cmath>
#include <iosfwd> // Just Forward decls for iostream.

class Int;
class Ratio;
template<class A> class Poly;

/*  Float precision:
 *  Operations are performd using the MPFR default precision that can be set by
 *  means of:
 *      * Float::SetPrecision (int n_digits, int radix);
 *      * set_precision<Float> (int n_digits, int radix);
 *      * push_precision<Float> (int n_digits, int radix);
 *      * push_precision<Float> ();
 *      * pop_precision<Float> ();
 *
 *  The precision of an individual Float can be set by
 *      * f.set_precision (n_bits);
 *
 *  For all functions, the chosen precision is saturated to the range
 *  supported by MPFR which is [MPFR_PREC_MIN, MPFR_PREC_MAX].
 *
 *  The function templates are also instanciated for double as no-ops.
 *
 *  Either of  f.copy(x)  and  f.move(x)  copy / move the input unaltered
 *  to the target, i.e. it will have the same value and precision.
 *  f.move(x)  assumes that  x  has been initialized properly using mpfr_init*,
 *  and the "ownership" will be transferred to the Float, i.e.  you MUST
 *  NOT CALL  mpfr_clear(x)  after you called f.move(x).
 */

class Float
{
    mpfr_t x_;
    bool clear_p_;
    void clear();

public:
    class Stat;

    static Stat stat;
    static int Debug;
    static bool use_fma;

    struct Precision { mpfr_prec_t prec; };

    ~Float();
    Float();
    Float (Float::Precision);
    Float (double);
    Float (const Ratio&);
    Float (const Ratio&, int prec2);

    Float (const Float&);
    Float (Float&&);
    Float& operator= (const Float&);
    Float& operator= (Float&&);

    explicit operator Int () const; // Round to nearest, with ties to even.
    explicit operator Ratio () const; // Round to nearest int.

    static Float from_Int (const Int&);
    static Float from_int (int);
    // In order to initialize Float from mpfr_t use
    // Float().copy(const mpfr_t) or Float().move(mpfr_t).
    const mpfr_t& mpfr() const { return x_; }
    mpfr_t& mpfr() { return x_; }
    Float& move (mpfr_t);
    // Exact copy: Same precision, same value, no rounding.
    Float& copy (const mpfr_t);
    Float& copy (const Float& z) { return copy (z.mpfr()); }
    // Use these if the target had been initialized.
    void move_to (mpfr_t);
    void copy_to (mpfr_t) const;
    // Use these if the target is not initialized yet.
    void init_move_to (mpfr_t);
    void init_copy_to (mpfr_t) const;

    explicit operator double() const;

    Float& set_precision_round (int bits, mpfr_rnd_t);
    mpfr_prec_t get_precision () const;
    static void SetPrecision (int digits, int radix);
    static int GetPrecision (int radix);
    static mpfr_prec_t SaturatePrecision (mpfr_prec_t);

    struct Format
    {
        const char* str = default_str;
        bool with_precision = true;
        static constexpr const char *default_str = "%RNe";
        static Format push();
        static Format push (const char*, bool = false);
        static Format push (const Format&);
        static Format pop();
        Format() {}
        Format (const char *str, bool prec = false)
            : str(str), with_precision(prec) {}
        friend std::ostream& operator << (std::ostream&, const Format&);
    };
    static Format format;

    Float operator - () const;
    Float operator + () const;
    Float operator + (const Float&) const;
    Float operator - (const Float&) const;
    Float operator * (const Float&) const;
    Float operator / (const Float&) const;
    Float fma (const Float&, const Float&) const; // a * b + c
    Float inv (long = 1) const; // long / Float
    Float operator % (const Float&) const;
    void operator += (const Float&);
    void operator -= (const Float&);
    void operator *= (const Float&);
    void operator /= (const Float&);
    void operator %= (const Float&);
    void operator <<= (int);
    void operator >>= (int);
    Float operator << (int) const;
    Float operator >> (int) const;

    int signbit () const; // 0 or != 0
    int sgn () const; // < 0, == 0 or > 0
    int cmp (const Float&) const;
    int cmp (long) const;
    int cmp (int) const;
    int cmp (double) const;

#define MK_OP(OP)                          \
    bool operator OP (const Float&) const; \
    bool operator OP (int) const;          \
    bool operator OP (long) const;         \
    bool operator OP (double) const
    MK_OP (==);
    MK_OP (!=);
    MK_OP (<=);
    MK_OP (>=);
    MK_OP (>);
    MK_OP (<);
#undef MK_OP

    bool nan_p () const;
    bool inf_p () const;
    bool number_p () const;
    bool zero_p () const;
    bool nonzero_p () const;

    static Float urandom ();
    static Float urandom (const Float&, const Float&);
    static const Float One;
    static const Float Zero;
    static const Float Nan;
    static const Float Inf;
    static Float pi();
    static Float e();
    static Float ln2();
    static Float ln10();
    static Float euler();
    static Float catalan();
    static Float apery();
    bool is_zero() const { return zero_p(); }
    bool is_integer() const { return mpfr_integer_p (x_); }

    Float round () const;
    Float ceil () const;
    Float floor () const;
    Float trunc () const;
    Float abs () const;
    Float sqrt () const;
    Float cbrt () const;
    Float exp () const;
    Float exp2 () const;
    Float exp10 () const;
    Float expm1 () const;   // e^x - 1
    //Float exp2m1 () const;  // 2^x - 1
    //Float exp10m1 () const; // 10^x - 1
    Float log () const;
    Float log2 () const;
    Float log10 () const;
    Float log1p () const;   // ln (1+x)
    //Float log2p1 () const;  // ld (1+x)
    //Float log10p1 () const; // lg (1+x)
    Float sin () const;
    Float cos () const;
    Float tan () const;
    Float sec () const;
    Float csc () const;
    Float cot () const;
    Float asin () const;
    Float acos () const;
    Float atan () const;
    Float sinh () const;
    Float cosh () const;
    Float tanh () const;
    Float sech () const;
    Float csch () const;
    Float coth () const;
    Float asinh () const;
    Float acosh () const;
    Float atanh () const;
    Float zeta () const;
    Float gamma () const;
    static Float factorial (unsigned);

    Float ldexp (int) const;
    Float pow (const Float&) const;
    Float root (unsigned long) const;
    Float atan2 (const Float&) const;
    Float min (const Float&) const;
    Float max (const Float&) const;
    Float fmod (const Float&) const;
    int abscmp (const Float&) const;
    Float absdiff (const Float&) const;
    Float copysign (const Float&) const;
    Float setsign (int) const;
    Float nextabove () const;
    Float nextbelow () const;
    Float nextafter (const Float&) const;
    Float nextafter (bool above_p) const;

    // 0.5 <= |return| < 1  and  return * 2^exp =~ x.
    Float frexp (int *exp) const;

    // x =~ z * 2^return
    // Return = emin for 0, Inf and NaN (mpfr_get/set_emin()).
    int z_2exp (Int &z, mpfr_prec_t, mpfr_rnd_t = MPFR_RNDN) const;

    Float powi (long) const;

    Float linear (const Float&, const Float&) const;
    // Relative position [0,1] in interval [a,b].
    Float where_in (const Float&, const Float&) const;
    Float saturate (const Float&, const Float&) const;
    bool in_range (const Float&, const Float&) const;

    static Float from_str (const char*, int base = 0);
    char* to_str (const char* = nullptr) const; // malloc.
    void print (const char*) const;
    std::ostream& print (std::ostream&, const char* = nullptr) const;
    friend std::ostream& operator << (std::ostream&, const Float&);
    void dump() const;

    // Returned flags like in F7_CONST_DEF: 0=sign, 1=zero (not usd), 2=NaN,
    // 7=Inf.. bytes[0] is MSB, e.g. 1 = 0x80,0,0,... expo=0.
    int as_bytes (int n_bytes, int *bytes, int *expo, int prec=-1) const;

    // Return IEEE number as bytes[] with N_EXPO bits for the biased exponent,
    // and PREC bits for the *encoded* mantissa.
    // Returned bytes are little endian.
    int* as_IEEE (int n_bytes, int n_expo,
                  mpfr_prec_t = 0, mpfr_rnd_t = MPFR_RNDN) const;

    Poly<Float> operator + (const Poly<Float> &q) const;
    Poly<Float> operator - (const Poly<Float> &q) const;
    Poly<Float> operator * (const Poly<Float> &q) const;

    static Float make (int);
};

// X func (X, Y);

#define MK_OP(NAME)                                 \
    using std::NAME;                                \
    Float NAME (const Float&, const Float&);
    MK_OP (pow)
    MK_OP (atan2)
    MK_OP (fmin)
    MK_OP (fmax)
    MK_OP (fmod)
    MK_OP (copysign)
    MK_OP (nextafter)
#undef MK_OP

using std::ldexp;
Float ldexp (const Float&, int);

using std::frexp;
Float frexp (const Float&, int*);

Float powi (const Float&, long);
Float root (const Float&, unsigned long);
int abscmp (const Float&, const Float&);
Float min (const Float&, const Float&);
Float max (const Float&, const Float&);
Float fma (const Float&, const Float&, const Float&);

// X func (X);

Float abs (const Float&);

#define MK_OP(NAME)                                 \
    using std::NAME;                                \
    Float NAME (const Float&);
    MK_OP (fabs)
    MK_OP (round)
    MK_OP (ceil)
    MK_OP (floor)
    MK_OP (trunc)
    MK_OP (sqrt)
    MK_OP (cbrt)
    MK_OP (exp)
    MK_OP (exp2)
    MK_OP (expm1) // e^x - 1
//    MK_OP (exp10)
    MK_OP (log)
    MK_OP (log2)
    MK_OP (log10)
    MK_OP (log1p) // ln (1+x)
    MK_OP (sin)
    MK_OP (cos)
    MK_OP (tan)
    MK_OP (asin)
    MK_OP (acos)
    MK_OP (atan)
    MK_OP (sinh)
    MK_OP (cosh)
    MK_OP (tanh)
    MK_OP (asinh)
    MK_OP (acosh)
    MK_OP (atanh)
    MK_OP (nextafter)
#undef MK_OP

Float gamma (const Float&);
Float sec (const Float&);
Float csc (const Float&);
Float cot (const Float&);
Float sech (const Float&);
Float csch (const Float&);
Float coth (const Float&);
//Float log2p1 (const Float&);  // ld (1+x)
//Float log10p1 (const Float&); // lg (1+x)
//Float exp2m1 (const Float&);  // 2^x - 1
//Float exp10m1 (const Float&); // 10^x - 1

// int func (X);

#define MK_OP(NAME)                                 \
    using std::NAME;                                \
    int NAME (const Float&);
    MK_OP (isnan)
    MK_OP (isinf)
    MK_OP (isfinite)
#undef MK_OP


template<class F>
void set_precision (int digits, int radix);

template<class F>
void push_precision (int digits, int radix);

template<class F>
void push_precision ();

template<class F>
void pop_precision();

template<class F>
F random_range (const F&, const F&);


class Float::Stat
{
public:
    int scale;

#ifndef FLOAT_NO_STAT
#define UNSCALED(NAME)              \
    private: int n_## NAME;         \
    public: void NAME (int inc = 1) \
    {                               \
        n_## NAME += inc;           \
    }
#define SCALED(NAME)                \
    private: int n_## NAME;         \
    public: void NAME (int inc = 1) \
    {                               \
        n_## NAME += scale * inc;   \
    }
#else
#define UNSCALED(NAME)              \
    private: int n_## NAME;         \
    public: void NAME (int = 1) {}
#define SCALED(NAME)                \
    private: int n_## NAME;         \
    public: void NAME (int = 1) {}
#endif // FLOAT_NO_STAT
    SCALED (add)
    SCALED (sub)
    SCALED (mul)
    SCALED (div)
    SCALED (fmod)
    SCALED (sqr)
    SCALED (sqrt)
    SCALED (neg)
    SCALED (cmp)
    UNSCALED (new_copy)
    UNSCALED (new_move)
    UNSCALED (new_other)
    UNSCALED (new_nan)
    UNSCALED (new_double)
    UNSCALED (new_Float)
    UNSCALED (free)
    UNSCALED (copy)
    UNSCALED (move)
#undef SCALED
#undef UNSCALED

    Stat();
    void reset();
    void print () const;
    void print (bool /* reset after print? */);
};

Float operator "" _Float (const char[], size_t);
Float operator "" _Float (const char*);
Float operator "" _R (const char[], size_t);
Float operator "" _R (const char*);

std::ostream& operator << (std::ostream&, const Float::Stat&);

template<class X>
X where_in (const X& x, const X& a, const X& b)
{
    return (x - a) / (b - a);
}

#endif // MP_FLOAT_H
