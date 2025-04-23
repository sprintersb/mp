#include "mp-float.h"
#include "mp-ratio.h"
#include "mp-int.h"

#include "diagnostic.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <cmath> // std::max, or <algorithm> ?
#include <ostream>

int Float::Debug = 0;
bool Float::use_fma = true;

void Float::clear()
{
    if (clear_p_)
    {
        mpfr_clear (x_);
        clear_p_ = false;
    }
}

Float::~Float()
{
    stat.free();
    clear();
}

Float::Float()
{
    mpfr_init (x_);
    clear_p_ = true;
    stat.new_nan();
}

Float::Float (Float::Precision p) : clear_p_(true)
{
    stat.new_nan();
    mpfr_init2 (x_, SaturatePrecision (p.prec)); // NaN with binary precision.
}

Float::Float (double d) : Float()
{
    stat.new_nan(-1);
    stat.new_double();
    mpfr_set_d (x_, d, MPFR_RNDN);
}

Float::Float (const Ratio& q) : Float()
{
    stat.new_nan(-1);
    stat.new_Float();
    mpfr_set_q (x_, q.mpq(), MPFR_RNDN);
}

Float Float::from_int (int i)
{
    Float x;
    mpfr_set_si (x.mpfr(), i, MPFR_RNDN);
    return x;
}

Float Float::from_Int (const Int& z)
{
    return z.operator Float();
}

Int Int::from_Float (const Float& x, int rnd)
{
    Int z;
    mpfr_get_z (z.mpz(), x.mpfr(), (mpfr_rnd_t) rnd);
    return z;
}

mpfr_prec_t Float::SaturatePrecision (mpfr_prec_t p)
{
    mpfr_prec_t prec = (mpfr_prec_t) p;
    prec = std::max<mpfr_prec_t> (MPFR_PREC_MIN, prec);
    prec = std::min<mpfr_prec_t> (MPFR_PREC_MAX, prec);
    return prec;
}

Float::Float (const Ratio& q, int bits) : clear_p_(true)
{
    stat.new_Float();
    mpfr_init2 (x_, Float::SaturatePrecision (bits));
    mpfr_set_q (x_, q.mpq(), MPFR_RNDN);
}

Float::Float (const Float& f) : Float()
{
    stat.new_nan(-1);
    stat.new_copy();
    stat.copy();
    mpfr_set (x_, f.x_, MPFR_RNDN);
}

Float::Float (Float&& f)
{
    stat.new_move();
    stat.move();
    std::memcpy (this, &f, sizeof (f));
    f.clear_p_ = false;
}

Float& Float::operator = (const Float& f)
{
    assert (this != &f);
    assert (clear_p_);
    stat.copy();
    mpfr_set_prec (x_, mpfr_get_prec (f.x_)); // Mimic = (Float&&)
    mpfr_set (x_, f.x_, MPFR_RNDN);
    return *this;
}

Float& Float::operator = (Float&& f)
{
    assert (this != &f);
    stat.move();
    clear();
    std::memcpy (this, &f, sizeof (f));
    f.clear_p_ = false;
    return *this;
}

Float& Float::move (mpfr_t f)
{
    clear();
    std::memcpy (&x_, f, sizeof (mpfr_t));
    clear_p_ = true;
    return *this;
}

Float& Float::copy (const mpfr_t z)
{
    clear();
    mpfr_init2 (x_, mpfr_get_prec (z));
    mpfr_set (x_, z, MPFR_RNDN);
    clear_p_ = true;
    return *this;
}

void Float::move_to (mpfr_t f)
{
    mpfr_clear (f);
    init_move_to (f);
}

void Float::init_move_to (mpfr_t f)
{
    assert (clear_p_);
    std::memcpy (f, &x_, sizeof (mpfr_t));
    clear_p_ = false;
}

namespace
{
    // Must be a function so pointers are decaying.
    bool _mpfr_same (const mpfr_t a, const mpfr_t b)
    {
        return a == b;
    }
};

void Float::copy_to (mpfr_t f) const
{
    //static_assert (sizeof (f) == sizeof (void*), "pointer should decay");
    if (! _mpfr_same (x_, f))
    {
        mpfr_set_prec (f, get_precision());
        mpfr_set (f, x_, MPFR_RNDN);
    }
}

void Float::init_copy_to (mpfr_t f) const
{
    //static_assert (sizeof (f) == sizeof (void*), "pointer should decay");
    if (! _mpfr_same (x_, f))
    {
        mpfr_init2 (f, get_precision());
        mpfr_set (f, x_, MPFR_RNDN);
    }
}


Float& Float::set_precision_round (int bits, mpfr_rnd_t r)
{
    mpfr_prec_round (x_, Float::SaturatePrecision (bits), r);
    return *this;
}

mpfr_prec_t Float::get_precision () const
{
    return mpfr_get_prec (x_);
}


void Float::SetPrecision (int digits, int radix)
{
    digits = std::max (2, digits);
    radix  = std::max (2, radix);
    int bits = digits;
    if (radix != 2)
        bits = 1 + bits * std::log2 (radix);

    mpfr_set_default_prec (Float::SaturatePrecision (bits));
}

int Float::GetPrecision (int radix)
{
    int bits = (int) mpfr_get_default_prec();
    return radix == 2
        ? bits
        : 1 + bits / std::log2 (radix);
}

template<>
void set_precision<Float> (int digits, int radix)
{
    Float::SetPrecision (digits, radix);
}

#include <stack>

namespace {
    std::stack<int> _precisions;
    std::stack<Float::Format> _formats;
}

template<>
void push_precision<Float> (int digits, int radix)
{
    _precisions.push (Float::GetPrecision (2));
    Float::SetPrecision (digits, radix);
}

template<>
void push_precision<Float> ()
{
    _precisions.push (Float::GetPrecision (2));
}

template<>
void pop_precision<Float> ()
{
    if (_precisions.size())
    {
        Float::SetPrecision (_precisions.top(), 2);
        _precisions.pop();
    }
}

template<>
void set_precision<double> (int, int)
{
    // Empty.
}

Float::Format Float::Format::push()
{
    _formats.push (Float::format);
    return Float::format;
}

Float::Format Float::Format::push (const char *fmt, bool with_prec)
{
    return push (Format { fmt ? fmt : default_str, with_prec });
}

Float::Format Float::Format::push (const Format& f)
{
    Format old = push();
    Float::format = f;
    return old;
}

Float::Format Float::Format::pop()
{
    Format old = Float::format;

    if (! _formats.empty())
    {
        Float::format = _formats.top();
        _formats.pop();
    }
    else
        Float::format = Format{};

    return old;
}

std::ostream& operator << (std::ostream& ost, const Float::Format& f)
{
    Float::format = f;
    return ost;
}


Float Float::make (int i)
{
    Float f;
    mpfr_set_si (f.x_, (long) i, MPFR_RNDN);
    return f;
}

const Float Float::Zero { 0 };
const Float Float::One { 1 };
const Float Float::Inf = "inf"_R;
const Float Float::Nan = "nan"_R;

Float::operator double() const
{
    return mpfr_get_d (x_, MPFR_RNDN);
}

Int::operator Float() const
{
    Float f;
    mpfr_set_z (f.mpfr(), mpz(), MPFR_RNDN);
    return f;
}

Float::operator Int () const
{
    Int z;
    mpfr_get_z (z.mpz(), mpfr(), MPFR_RNDN);
    return z;
}

Float::operator Ratio () const
{
    return Ratio { operator Int() };
}

Float Float::operator - () const
{
    stat.neg();
    Float z;
    mpfr_neg (z.x_, x_, MPFR_RNDN);
    return z;
}

Float Float::operator + () const
{
    return *this;
}


#define MK_OP(NAME, OP)                             \
    Float Float::operator OP (const Float& y) const \
    {                                               \
        stat.NAME ();                               \
        Float z;                                    \
        mpfr_## NAME (z.x_, x_, y.x_, MPFR_RNDN);   \
        return z;                                   \
    }
    MK_OP (add, +)
    MK_OP (sub, -)
    MK_OP (div, /)
    MK_OP (fmod, %)
#undef MK_OP

Float Float::inv (long i) const
{
    stat.div();
    Float z;
    mpfr_si_div (z.x_, i, x_, MPFR_RNDN);
    return z;
}

void Float::operator *= (const Float& y)
{
    if (this == &y)
    {
        stat.sqr();
        mpfr_sqr (x_, x_, MPFR_RNDN);
    }
    else
    {
        stat.mul();
        mpfr_mul (x_, x_, y.x_, MPFR_RNDN);
    }
}

Float Float::operator * (const Float& y) const
{
    Float z;
    if (this == &y)
    {
        stat.sqr();
        mpfr_sqr (z.x_, x_, MPFR_RNDN);
    }
    else
    {
        stat.mul();
        mpfr_mul (z.x_, x_, y.x_, MPFR_RNDN);
    }
    return z;
}

Float Float::fma (const Float& a, const Float& b) const
{
    stat.mul();
    stat.add();
    if (Float::use_fma)
    {
        Float z;
        mpfr_fma (z.x_, x_, a.x_, b.x_, MPFR_RNDN);
        return z;
    }
    else
        return (*this) * a + b;
}

Float fma (const Float& x, const Float& a, const Float& b)
{
    return x.fma (a, b);
}


#define MK_OP(NAME, OP)                             \
    void Float::operator OP (const Float& y)        \
    {                                               \
        stat.NAME ();                               \
        mpfr_## NAME (x_, x_, y.x_, MPFR_RNDN);     \
    }
    MK_OP (add, +=)
    MK_OP (sub, -=)
    MK_OP (div, /=)
    MK_OP (fmod, %=)
#undef MK_OP

int Float::cmp (const Float& y) const
{
    stat.cmp();
    return mpfr_cmp (x_, y.x_);
}

int Float::cmp (int i) const
{
    stat.cmp();
    return mpfr_cmp_si (x_, (long) i);
}

int Float::cmp (long l) const
{
    stat.cmp();
    return mpfr_cmp_si (x_, l);
}

int Float::cmp (double d) const
{
    stat.cmp();
    return mpfr_cmp_d (x_, d);
}

int Float::sgn () const
{
    return mpfr_sgn (x_);
}

#define MK_OP(OP, NAME)                             \
    bool Float::operator OP (const Float& y) const  \
    {                                               \
        return mpfr_## NAME (x_, y.x_);             \
    }
    MK_OP (==, equal_p)
    MK_OP (!=, lessgreater_p)
    MK_OP (>,  greater_p)
    MK_OP (>=, greaterequal_p)
    MK_OP (<,  less_p)
    MK_OP (<=, lessequal_p)
#undef MK_OP

#define MK_OP(OP)                                   \
    bool Float::operator OP (int i) const           \
    {                                               \
        return mpfr_cmp_si (x_, (long) i) OP 0;     \
    }                                               \
    bool Float::operator OP (long l) const          \
    {                                               \
        return mpfr_cmp_si (x_, l) OP 0;            \
    }                                               \
    bool Float::operator OP (double d) const        \
    {                                               \
        return mpfr_cmp_d (x_, d) OP 0;             \
    }
    MK_OP (==)
    MK_OP (!=)
    MK_OP (>)
    MK_OP (<)
    MK_OP (>=)
    MK_OP (<=)
#undef MK_OP

#define MK_OP(METH, NAME)                           \
    bool Float::METH () const                       \
    {                                               \
        return mpfr_## NAME (x_);                   \
    }
    MK_OP (nan_p, nan_p)
    MK_OP (inf_p, inf_p)
    MK_OP (number_p, number_p)
    MK_OP (zero_p, zero_p)
    MK_OP (nonzero_p, regular_p)
#undef MK_OP

void Float::operator <<= (int ex)
{
    if (number_p())
        mpfr_mul_2si (x_, x_, ex, MPFR_RNDN);
}
void Float::operator >>= (int ex)
{
    if (number_p())
        mpfr_mul_2si (x_, x_, -ex, MPFR_RNDN);
}

Float Float::operator << (int ex) const
{
    return ldexp (ex);
}
Float Float::operator >> (int ex) const
{
    return ldexp (-ex);
}

Float Float::ldexp (int ex) const
{
    Float z { *this };
    if (z.nonzero_p ())
        mpfr_mul_2si (z.x_, z.x_, ex, MPFR_RNDN);
    return z;
}
Float ldexp (const Float& f, int ex)
{
    return f.ldexp (ex);
}

#define MK_OP(NAME, OP)                             \
    Float Float::NAME (const Float& y) const        \
    {                                               \
        Float z;                                    \
        mpfr_##OP (z.x_, x_, y.x_, MPFR_RNDN);      \
        return z;                                   \
    }                                               \
    Float NAME (const Float& x, const Float& y)     \
    {                                               \
        return x.NAME (y);                          \
    }
    MK_OP (pow, pow)
    MK_OP (min, min)
    MK_OP (max, max)
    MK_OP (atan2, atan2)
    MK_OP (fmod, fmod)
    MK_OP (copysign, copysign)
    MK_OP (absdiff, dim)
#undef MK_OP


// nextafter() has signature (F,F) while
// nexttoward has signature (F,long double).
// MPFR's mpfr_nexttoward(F,F)  deviates from that, so for now we only
// provide nextafter.

Float Float::nextafter (const Float& y) const
{
    Float z = *this;
    mpfr_nexttoward (z.x_, y.x_);
    return z;
}
Float nextafter (const Float& x, const Float& y)
{
    return x.nextafter (y);
}

Float Float::frexp (int *ex2) const
{
    mpfr_exp_t mexp;
    Float z;
    mpfr_frexp (&mexp, z.x_, x_, MPFR_RNDN);
    *ex2 = (int) mexp;
    return z;
}
Float frexp (const Float &x, int *ex2)
{
    return x.frexp (ex2);
}

Float Float::powi (long n) const
{
    Float z;
    mpfr_pow_si (z.x_, x_, n, MPFR_RNDN);
    return z;
}
Float powi (const Float& x, long n)
{
    return x.powi (n);
}

Float Float::root (unsigned long ex) const
{
    Float z;
    mpfr_rootn_ui (z.x_, x_, ex, MPFR_RNDN);
    return z;
}
Float root (const Float& f, unsigned long ex)
{
    return f.root (ex);
}

int Float::abscmp (const Float& y) const { return mpfr_cmpabs (x_, y.x_); }

int abscmp (const Float& x, const Float& y) { return x.abscmp (y); }

Float Float::setsign (int s) const
{
    Float z;
    mpfr_setsign (z.x_, x_, s, MPFR_RNDN);
    return z;
}

int Float::signbit() const
{
    return mpfr_signbit (mpfr());
}

Float Float::nextabove () const
{
    Float x = *this;
    mpfr_nextabove (x.x_);
    return x;
}

Float Float::nextbelow () const
{
    Float x = *this;
    mpfr_nextbelow (x.x_);
    return x;
}

Float Float::nextafter (bool above_p) const
{
    return above_p ? nextabove() : nextbelow();
}

Float Float::linear (const Float& x0, const Float& x1) const
{
    const Float& t = *this;
    return x0 + t * (x1 - x0);
}

Float Float::where_in (const Float& a, const Float& b) const
{
    return ((*this) - a) / (b - a);
}

Float Float::saturate (const Float& lo, const Float& hi) const
{
    return max(lo).min(hi);
}

bool Float::in_range (const Float& lo, const Float& hi) const
{
    return (*this) >= lo && (*this) <= hi;
}

Float Float::pi()  { Float f; mpfr_const_pi (f.x_, MPFR_RNDN); return f; } // 3.141...
Float Float::ln2() { Float f; mpfr_const_log2 (f.x_, MPFR_RNDN); return f; } // 0.693...
Float Float::ln10() { return Float::from_int(10).log(); } // 2.302...
Float Float::euler() { Float f; mpfr_const_euler (f.x_, MPFR_RNDN); return f; } // 0.577...
Float Float::apery() { Float f; mpfr_zeta_ui (f.x_, 3, MPFR_RNDN); return f; } // 1.202...
Float Float::catalan() { Float f; mpfr_const_catalan (f.x_, MPFR_RNDN); return f; } // 0.915...
Float Float::e()   { return Float::from_int(1).exp(); } // 2.718...


#define MK_OP(NAME)                                 \
    Float Float::NAME () const                      \
    {                                               \
        Float z;                                    \
        mpfr_## NAME (z.x_, x_);                    \
        return z;                                   \
    }                                               \
    Float NAME (const Float &y)                     \
    {                                               \
        return y.NAME();                            \
    }
    MK_OP (round)
    MK_OP (ceil)
    MK_OP (floor)
    MK_OP (trunc)
#undef MK_OP

#define MK_OP(NAME)                                 \
    Float Float::NAME () const                      \
    {                                               \
        Float z;                                    \
        mpfr_## NAME (z.x_, x_, MPFR_RNDN);         \
        return z;                                   \
    }                                               \
    Float NAME (const Float &y)                     \
    {                                               \
        return y.NAME();                            \
    }
    MK_OP (abs)
    MK_OP (sqrt)
    MK_OP (cbrt)
    MK_OP (exp)
    MK_OP (exp2)
    MK_OP (exp10)
    MK_OP (expm1)
    //MK_OP (exp2m1)
    //MK_OP (exp10m1)
    MK_OP (log)
    MK_OP (log2)
    MK_OP (log10)
    MK_OP (log1p)
    //MK_OP (log2p1)
    //MK_OP (log10p1)
    MK_OP (sin)
    MK_OP (cos)
    MK_OP (tan)
    MK_OP (gamma)
    MK_OP (sec)
    MK_OP (csc)
    MK_OP (cot)
    MK_OP (asin)
    MK_OP (acos)
    MK_OP (atan)
    MK_OP (sinh)
    MK_OP (cosh)
    MK_OP (tanh)
    MK_OP (sech)
    MK_OP (csch)
    MK_OP (coth)
    MK_OP (asinh)
    MK_OP (acosh)
    MK_OP (atanh)
    MK_OP (zeta)
#undef MK_OP

#define MK_OP(NAME, METH)                           \
    int NAME (const Float& f)                       \
    {                                               \
        return f.METH();                            \
    }
    MK_OP (isnan, nan_p)
    MK_OP (isinf, inf_p)
    MK_OP (isfinite, number_p)
#undef MK_OP

Float fmin (const Float& a, const Float& b) { return a.min (b); }
Float fmax (const Float& a, const Float& b) { return a.max (b); }
Float fabs (const Float& a) { return a.abs(); }

Float Float::factorial (unsigned u)
{
    Float f;
    mpfr_fac_ui (f.x_, (unsigned long) u, MPFR_RNDN);
    return f;
}

char* Float::to_str (const char *fmt) const
{
    char *str;
    mpfr_asprintf (&str, fmt ? fmt : Float::format.str, x_, x_);
    return str;
}

void Float::print (const char* fmt) const
{
    mpfr_printf (fmt ? fmt : Float::format.str, x_);
}

std::ostream& Float::print (std::ostream& ost, const char* fmt) const
{
    bool with_prec = !fmt && Float::format.with_precision;
    if (with_prec)
        ost << "[" << get_precision() << "]=";

    char *str = to_str();
    ost << str;
    std::free (str);

    return ost;
}

void Float::dump() const
{
    std::printf ("<@0x%p:clear_p_=%d", this, clear_p_);
    std::printf (":prec=%d", (int) x_[0]._mpfr_prec);
    std::printf (":sign=%d", (int) x_[0]._mpfr_sign);
    std::printf (":expo=%ld", (long) x_[0]._mpfr_exp);
    std::printf (":d=%p>", x_[0]._mpfr_d);
    out ("\n");
}

Float Float::from_str (const char *str, int base)
{
    Float f;

    char upper = '0' + (base ? base : 10);
    if (! str[0])
        { /* NaN */ }
    else if (str[0] >= '0' && ! str[1] && str[0] < upper)
        mpfr_set_ui (f.x_, (unsigned long) (str[0] - '0'), MPFR_RNDN);
    else if (str[0] == '-' && str[1] >= '0' && ! str[2] && str[1] < upper)
    {
        if (str[1] == '0')
            mpfr_set_d (f.x_, -0.0, MPFR_RNDN);
        else
            mpfr_set_si (f.x_, (long) ('0' - str[1]), MPFR_RNDN);
    }
    else
        mpfr_set_str (f.x_, str, base, MPFR_RNDN);

    return f;
}

// x =~ z * 2^return
// Return = emin for 0, Inf and NaN (mpfr_get/set_emin()).
int Float::z_2exp (Int &z, mpfr_prec_t prec, mpfr_rnd_t rnd) const
{
    mpfr_t op;
    mpfr_init2 (op, prec);
    mpfr_set (op, mpfr(), rnd);

    mpfr_exp_t ex = mpfr_get_z_2exp (z.mpz(), op);

    mpfr_clear (op);

    return (int) ex;
}

// Convert to IEEE representation with N_EXPO exponent bits
// PREC bits for the *encoded* mantissa.
// Returned bytes are little endian.
int* Float::as_IEEE (int n_bytes, int n_expo,
                     mpfr_prec_t prec, mpfr_rnd_t rnd) const
{
    if (prec == 0)
        prec = 8 * n_bytes - 1 - n_expo;

    if (n_expo < 2)
        error ("n_expo = %d", n_expo);
    if (8 * n_bytes < 2 + n_expo)
        error ("no space for mantissa: %d bytes, expo = %d bits",
               n_bytes, n_expo);

    const long expo_mask = (1L << n_expo) - 1;
    const long expo_bias = expo_mask >> 1;
    // Bit position of the implicit 1. of the mantissa.
    const int pos1 = 8 * n_bytes - n_expo - 1;
    // Normal and != 0 ?
    bool is_normal = false;

    // Mantissa bits, leading 1. already removed for normals.
    Int mant = 0_Int;

    // Biased exponent.
    long expo = 0;

    if (nan_p())
    {
        // expo = expo_max; mant != 0
        expo = expo_mask;
        mant.setbit (pos1 - 1);
    }
    else if (inf_p())
    {
        // expo = expo_max; mant = 0
        expo = expo_mask;
    }
    else if (! zero_p()) // ordinary number
    {
        // Biased expo is in 1 ... expo_mask-1.
        //    "    "    " 0 for subnormals.
        //    "    "    " expo_mask for Inf.
        expo = z_2exp (mant, 1 + prec, rnd);
        mant = mant.abs();
        //gmp_printf ("mant = 0x %Zx, expo = %d\n", mant.mpz(), expo);

        if (! mant.is_zero())
        {
            int msbit = mant.msbit();

            expo += expo_bias + msbit;

            if (expo >= expo_mask)
            {
                // Inf: expo = expo_max; mant = 0
                expo = expo_mask;
            }
            else if (expo <= 0)
            {
                // Subnormal.
                // The mantissa has less bits than anticipated.
                // Get expo and mant again, but with the now known
                // number of mantissa bits.
                int mant_msb = pos1 - 1 + expo;
                int mant_bits = 1 + mant_msb;
                int mant_lsb = pos1 - prec;

                if (mant_lsb < 0 && mant_msb >= 0)
                    expo = z_2exp (mant, mant_bits, rnd);
                else if (mant_lsb >= 0 && mant_lsb <= mant_msb)
                    expo = z_2exp (mant, 1 + mant_msb - mant_lsb, rnd);
                else
                {
                    expo = - expo_bias;
                    mant = 0_Int;
                }

                mant = mant.abs();
                msbit = mant.msbit();
                expo += expo_bias + msbit;

                assert (mant == 0 || expo <= 1);

                if (mant.is_zero())
                    expo = 0;
                else if (expo <= 0)
                {
                    // Remains subnormal
                    int p1 = pos1 - 1 + expo;
                    mant = p1 >= msbit
                        ? mant << (p1 - msbit)
                        : mant >> (msbit - p1);
                    expo = 0;
                }
                else
                {
                    // Rounding turned submormal into normal.
                    assert (expo == 1);
                    is_normal = true;
                }
            }
            else // Normal number
            {
                is_normal = true;
            }
        } // mant != 0
    }

    if (is_normal)
    {
        int msbit = mant.msbit();
        mant = mant.setbit (msbit, 0);
        mant = pos1 >= msbit
            ? mant << (pos1 - msbit)
            : mant >> (msbit - pos1);
    }

    assert (0 <= expo && expo <= expo_mask);
    assert (0 <= mant && mant < 0_Int .setbit (pos1));

    Int z = signbit();
    z = (z << n_expo) | expo;
    z = (z << pos1) | mant;

    int *bytes = new int[n_bytes];

    for (int i = 0; i < n_bytes; ++i)
    {
        bytes[i] = (int) (long) (z & 0xff);
        z >>= 8;
    }

    return bytes;
}

// Return flags as in F7_CONST_DEF,  MSB = bytes[0].
int Float::as_bytes (int n_bytes, int *bytes, int *expo, int prec) const
{
    if (prec < 0)
        prec = 8 * n_bytes;

    *expo = 0;
    __builtin_memset (bytes, 0, sizeof (int) * n_bytes);

    int sign = sgn() < 0;

    if (nan_p())
        return 1 << 2;
    else if (inf_p())
        return sign | (1 << 7);
    else if (zero_p())
        return 0;

    Int z;
    // Write *this as z * 2 ^ *expo.
    *expo = z_2exp (z, prec);
    z = z.abs();

    int ms = z.msbit();
    assert (ms == prec - 1);
    // Adjust decimal point to one *after* the MSB, like used by LibF7.
    *expo += ms;

    if (z.is_zero() || *expo < INT16_MIN)
        return 0;
    else if (*expo > INT16_MAX)
        return sign | (1 << 7);

    for (int i = n_bytes - 1; i >= 0; --i)
    {
        bytes[i] = (int) (long) (z & 0xff);
        z >>= 8;
    }

    return sign;
}

Float operator "" _Float (const char *str, size_t)
{
    return Float::from_str (str);
}

Float operator "" _Float (const char *str)
{
    return Float::from_str (str);
}

Float operator "" _R (const char *str, size_t)
{
    return Float::from_str (str);
}

Float operator "" _R (const char *str)
{
    return Float::from_str (str);
}

std::ostream& operator << (std::ostream &ost, const Float &f)
{
    return f.print (ost);
}


// Apéry's constant: zeta(3) with Riemann's zeta function.
//
//      zeta_3 = sum_{k=1}^oo 1 / k^3
//             = 2.5 * sum_{k=1}^oo (-1)^{k+1} / d_k;   d_k = binom(2k,k) * k^3
//
// For the denominators we have:  log_2 (d_k) ~ 2k + 2.5 * log_2 k > 2k
// i.e. for  k  bits of precision we need about  k/2  summands, and these
// will introduce an error of ~ 0.5 LSB each.  Hence, we have a total error
// of around  k/4  LSBs, and therefore bump the precision by  log_2 k  bits.

#if 0 // Superseded by apery -> mpfr_zeta_ui.

Float Float::zeta_3()
{
    int prec = Float::GetPrecision (2 /* radix */);
    push_precision<Float> (prec + 3 + (int) std::log2 (prec), 2 /* radix */);

    Float sum{0}, old_sum;

    for (int k = 1; ; ++k)
    {
        Int b { Int::binom (2*k, k) * k * k * k };
        Float f { b.operator Float() };

        sum += f.inv (k % 2 == 0 ? -1 : 1);

        if (old_sum == sum)
            break;
        old_sum = sum;
    }

    sum *= 2.5;

    pop_precision<Float>();
    sum.set_precision_round (prec, MPFR_RNDN);

    return sum;
}
#endif // 0

Float::Format Float::format{};


Float Float::urandom ()
{
    static gmp_randstate_t state;
    static bool init_done;

    if (!init_done)
    {
        gmp_randinit_mt (state);
        init_done = true;
    }

    Float x;
    mpfr_urandom (x.x_, state, MPFR_RNDN);
    return x;
}

Float Float::urandom (const Float& x0, const Float& x1)
{
    const Float *p0 = &x0;
    const Float *p1 = &x1;
    if (x0 > x1)
        std::swap (p0, p1);

    const Float& y0 = *p0;
    const Float& y1 = *p1;

    return y0 + Float::urandom() * (y1 - y0);
}

template<>
Float random_range (const Float& x0, const Float& x1)
{
    return Float::urandom (x0, x1);
}

#include <cstdlib>

template<>
double random_range (const double& x0, const double& x1)
{
    double f = (double) rand() / RAND_MAX;

    return x1 > x0
        ? x0 + f * (x1 - x0)
        : x1 + f * (x0 - x1);
}

Float::Stat Float::stat;

Float::Stat::Stat()
{
    scale = 1;
    reset();
}

void Float::Stat::reset ()
{
    int s = scale;
    std::memset (this, 0, sizeof (*this));
    scale = s;
}

void Float::Stat::print (bool reset_after)
{
    print();
    if (reset_after)
        reset();
}

void Float::Stat::print () const
{
    int n = n_new_copy + n_new_move + n_new_double + n_new_nan;
    double p1 = n ? 100. * n_new_nan / n : 0;
    double p2 = n ? 100. * n_new_double / n : 0;
    printf (": new  % 6d (nan %.2f%%, double %.2f%%), copy %d, move %d\n",
            n, p1, p2, n_copy, n_move);
    printf (": free % 6d\n", n_free);
    n = n_neg + n_add + n_sub + n_mul + n_sqr + n_div + n_fmod;
    printf (": +-*/%%: %d\n", n);
    n = n_mul + n_sqr;
    p1 = n ? 100. * n_sqr / n : 0;
    printf (": *: %d (square %d = %.2f%%)\n", n, n_sqr, p1);
}

std::ostream& operator << (std::ostream& out, const Float::Stat& s)
{
    s.print();
    return out;
}
