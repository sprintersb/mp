#define DIAGNOSTIC_NO_FORMAT_CHECK
#include "diagnostic.h"

#include "mp-int.h"

#include <cstring>
#include <cassert>
#include <cstdlib>
#include <ostream>

void Int::clear()
{
    if (clear_p_)
    {
        mpz_clear (z_);
        clear_p_ = false;
    }
}

Int::Int (long i) : clear_p_(true)
{
    mpz_init_set_si (z_, i);
}

Int& Int::copy (const mpz_t z)
{
    mpz_set (z_, z);
    return *this;
}

Int& Int::move (mpz_t z)
{
    clear();
    std::memcpy (&z_, z, sizeof (mpz_t));
    clear_p_ = true;
    return *this;
}

void Int::move_to (mpz_t z)
{
    assert (clear_p_);
    mpz_clear (z);
    init_move_to (z);
}

void Int::init_move_to (mpz_t z)
{
    assert (clear_p_);
    std::memcpy (z, &z_, sizeof (mpz_t));
    clear_p_ = false;
}

void Int::copy_to (mpz_t z) const
{
    mpz_set (z, z_);
}

void Int::init_copy_to (mpz_t z) const
{
    mpz_init (z);
    copy_to (z);
}

Int Int::from_str (const char *str, int base)
{
    Int z;
    return z.set_str (str, base);
}

Int& Int::set_str (const char *str, int base)
{
    mpz_set_str (z_, str, base);
    return *this;
}

Int operator "" _Int (const char *str, size_t)
{
    return Int::from_str (str);
}

Int operator "" _Int (const char *str)
{
    return Int::from_str (str);
}

Int operator "" _Z (const char *str, size_t)
{
    return Int::from_str (str);
}

Int operator "" _Z (const char *str)
{
    return Int::from_str (str);
}

Int::Int (const Int& z) : clear_p_(true)
{
    mpz_init_set (z_, z.z_);
}

Int::Int (Int&& z)
{
    std::memcpy (this, &z, sizeof (z));
    z.clear_p_ = false;
}

Int& Int::operator = (const Int& z)
{
    assert (this != &z);
    assert (clear_p_);
    mpz_set (z_, z.z_);
    return *this;
}

Int& Int::operator = (Int&& z)
{
    assert (this != &z);
    clear();
    std::memcpy (this, &z, sizeof (z));
    z.clear_p_ = false;
    return *this;
}

Int Int::make (int i)
{
    return Int ((long) i);
}

#include <climits>

Int::operator long() const
{
    if (mpz_fits_slong_p (z_))
        return mpz_get_si (z_);
    return mpz_cmp_si (z_, LONG_MAX) >= 0 ? LONG_MAX : LONG_MIN;
}

Int::operator double() const
{
    return mpz_get_d (z_);
}

char* Int::to_str (int base) const
{
    // By default GMP uses malloc, realloc and free.
    return mpz_get_str (nullptr /* alloc */, base, z_);
}

void Int::print (FILE *stream, int base) const
{
    // By default GMP uses malloc, realloc and free.
    char *str = to_str (base);
    fputs (str, stream);
    free (str);
}

Int Int::operator - () const
{
    Int w;
    mpz_neg (w.z_, z_);
    return w;
}

Int Int::operator ~ () const
{
    Int w;
    mpz_com (w.z_, z_);
    return w;
}

Int& Int::operator ++ () // pre
{
    mpz_add_ui (z_, z_, 1ul);
    return *this;
}

Int& Int::operator -- () // pre
{
    mpz_sub_ui (z_, z_, 1ul);
    return *this;
}

Int Int::operator ++ (int) // post
{
    Int i { *this };
    ++(*this);
    return i;
}

Int Int::operator -- (int) // post
{
    Int i { *this };
    --(*this);
    return i;
}

#define MAKE_OP(OP, NAME, DIV0TXT)              \
    void Int::operator OP ##= (const Int& z)    \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%Zd%s", z_, "" DIV0TXT);    \
        mpz_## NAME (z_, z_, z.z_);             \
    }                                           \
    Int Int::operator OP (const Int& z) const   \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%Zd%s", z_, "" DIV0TXT);    \
        Int w;                                  \
        mpz_## NAME (w.z_, z_, z.z_);           \
        return w;                               \
    }                                           \
    Int operator OP (long i, const Int& z)      \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%ld%s", i, "" DIV0TXT);     \
        return Int{i} OP z;                     \
    }
    MAKE_OP (+, add, "")
    MAKE_OP (-, sub, "")
    MAKE_OP (*, mul, "")
    MAKE_OP (/, tdiv_q, " / 0") // "t" = round to zero, rem has sign of divisor d.
    MAKE_OP (%, tdiv_r, " % 0") // n = q * d + r, 0 <= |r| < |d|.
    MAKE_OP (&, and, "")
    MAKE_OP (|, ior, "")
    MAKE_OP (^, xor, "")
#undef MAKE_OP

#define MAKE_OP(OP)                             \
    bool Int::operator OP (const Int& z) const  \
    {                                           \
        return mpz_cmp (z_, z.z_) OP 0;         \
    }                                           \
    bool Int::operator OP (long i) const        \
    {                                           \
        return mpz_cmp_si (z_, i) OP 0;         \
    }                                           \
    bool operator OP (long i, const Int& z)     \
    {                                           \
        return Int{i} OP z;                     \
    }
    MAKE_OP (==)
    MAKE_OP (!=)
    MAKE_OP (>)
    MAKE_OP (>=)
    MAKE_OP (<)
    MAKE_OP (<=)
#undef MAKE_OP

long Int::exact_log2 () const
{
    return (*this) > 0 && 1 == popcount()
        ? (long) mpz_sizeinbase (z_, 2) - 1
        : -1L;
}

long Int::ilog2 () const
{
    if ((*this) <= 0)
        return -1;
    return (long) mpz_sizeinbase (z_, 2) - 1;
}

long Int::bitsize () const
{
    return (*this) < 0
        ? -1L
        : (long) mpz_sizeinbase (z_, 2);
}

Int Int::setbit (int bitno, bool val) const
{
    Int z (*this);
    if (val)
        mpz_setbit (z.mpz(), bitno);
    else
        mpz_clrbit (z.mpz(), bitno);
    return z;
}

Int Int::factorial (unsigned long n, int step)
{
    assert (step >= 1);
    Int z;
    switch (step)
    {
        default:
            error ("todo: step = %d", step);
            break;
        case 1:
            mpz_fac_ui (z.z_, n);
            break;

        case 2:
            mpz_2fac_ui (z.z_, n);
            break;
    }
    return z;
}


Int Int::binom (long n, long k)
{
    if (k > n || k < 0 || n < 0)
        return 0;

    Int z;
    mpz_bin_uiui (z.z_, (long) n, (long) k);
    return z;
#if 0
    if (k > n - k)
        k = n - k;

    Int b{1};

    for (auto i = n - k + 1; i <= n; ++i)
        b *= i;

    return b / Int::factorial (k);
#endif // 0
}


#define MAKE_OP(NAME, OP)                       \
    Int Int::NAME () const                      \
    {                                           \
        Int w;                                  \
        mpz_## OP (w.z_, z_);                   \
        return w;                               \
    }                                           \
    Int NAME (const Int& z)                     \
    {                                           \
        return z.NAME ();                       \
    }
    MAKE_OP (abs, abs)
    MAKE_OP (sqrt, sqrt)
    MAKE_OP (next_prime, nextprime)
#undef MAKE_OP


#define MAKE_OP(NAME, OP, DIV0TXT)              \
    Int Int::NAME (const Int& z) const          \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%Zd%s", z_, "" DIV0TXT);    \
        Int w;                                  \
        mpz_## OP (w.z_, z_, z.z_);             \
        return w;                               \
    }                                           \
    Int NAME (const Int& z, const Int& w)       \
    {                                           \
        if (*DIV0TXT && w.is_zero())            \
            error ("%Zd%s", z.mpz(), "" DIV0TXT);\
        return z.NAME (w);                      \
    }
    MAKE_OP (trunc_div, tdiv_q, " trunc_div 0")
    MAKE_OP (floor_div, fdiv_q, " floor_div 0")
    MAKE_OP (ceil_div, cdiv_q, " ceil_div 0")
    MAKE_OP (gcd, gcd, "")
    MAKE_OP (lcm, lcm, "")
    MAKE_OP (mod, mod, "")
#undef MAKE_OP


Int Int::round_div (const Int& z) const
{
    if (z.is_zero())
        error ("%Zd round_div 0", z_);
    return ((*this) + z.abs() / 2).floor_div (z);
}

Int round_div (const Int& z, const Int& w)
{
    if (w.is_zero())
        error ("%Zd round_div 0", z.mpz());
    return z.round_div (w);
}

Int Int::divmod (const Int& a, Int& rem) const
{
    if (a.is_zero())
        error ("%Zd divmod 0", mpz());
    Int q;
    mpz_fdiv_qr (q.z_, rem.z_, z_, a.z_); // "f"loor div.
    return q;
}

Int Int::pow (unsigned long ex) const
{
    Int w;
    mpz_pow_ui (w.z_, z_, ex);
    return w;
}

Int Int::min (const Int& z) const { return (*this) < z ? (*this) : z; }
Int Int::max (const Int& z) const { return (*this) > z ? (*this) : z; }

Int min (const Int& z, const Int& w) { return z.min (w); }
Int max (const Int& z, const Int& w) { return z.max (w); }

Int Int::operator << (mp_bitcnt_t shift) const
{
    Int w;
    mpz_mul_2exp (w.z_, z_, shift);
    return w;
}

Int Int::operator >> (mp_bitcnt_t shift) const
{
    Int w;
    mpz_tdiv_q_2exp (w.z_, z_, shift);
    return w;
}

void Int::operator <<= (mp_bitcnt_t shift)
{
    *this = (*this) << shift;
}

void Int::operator >>= (mp_bitcnt_t shift)
{
    *this = (*this) >> shift;
}


Int Int::remove_factor (const Int& d, int *ord) const
{
    assert (d.abscmp (2) >= 0);
    assert (!is_zero());

    Int q;
    int o = (int) mpz_remove (q.z_, z_, d.z_);
    if (ord)
        *ord = o;

    return q;
}


int cmp (const Int& z, const Int& w) { return z.cmp (w); }

int sgn (const Int& z) { return z.sgn(); }
bool is_perfect_square (const Int& z) { return z.is_perfect_square(); }
int is_probab_prime (const Int& i, int tries) { return i.is_probab_prime (tries); }
Int pow (const Int& z, unsigned long ex) { return z.pow (ex); }

Int Int::random (const Int& n)
{
    if (n < 2)
        return Int (0L);

/*
    ??? mpz_urandomm and mpz_urandomb are crashing ???
    gmp_randstate_t rstat;
    static bool inited;
    if (!inited)
    {
        inited = true;
        gmp_randinit_default (rstat);
        gmp_randseed_ui (rstat, 42);
        out ("GMP randinit\n");
    }
    ??? mpz_urandomm and mpz_urandomb are crashing ???
*/
    Int x, b { RAND_MAX };
    do
        x = x * b + rand();
    while (x < n.abs());

    return x.mod (n);
}


std::ostream& operator << (std::ostream& out, const Int& z)
{
    z.print (stdout, 10);
    return out;
}


const Int Int::Zero (0L);
const Int Int::One (1L);
const Int Int::Two (2L);
