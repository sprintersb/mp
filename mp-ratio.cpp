#define DIAGNOSTIC_NO_FORMAT_CHECK
#include "diagnostic.h"

#include "mp-ratio.h"
#include "mp-int.h"
#include "mp-float.h"

#include <cstring>
#include <cassert>
#include <cstdlib>
#include <ostream>

Ratio::~Ratio ()
{
    clear();
}

void Ratio::clear()
{
    if (clear_p_)
    {
        mpq_clear (q_);
        clear_p_ = false;
    }
}

void Ratio::init()
{
    assert (!clear_p_);
    clear_p_ = true;
    mpq_init (q_);
}

Ratio::Ratio ()
{
    init();
    mpq_set_ui (q_, 0, 1);
}

Ratio::Ratio (const Int& z)
{
    init();
    mpq_set_z (q_, z.mpz());
}

Ratio::Ratio (const Int& z, const Int& n)
{
    init();
    if (n.is_zero())
        error ("%Zd / 0", z.mpz());
    mpq_set_num (q_, z.mpz());
    mpq_set_den (q_, n.mpz());
    mpq_canonicalize (q_);
}

Ratio& Ratio::copy (const mpq_t q)
{
    mpq_set (q_, q);
    return *this;
}

Ratio& Ratio::move (mpq_t q)
{
    clear();
    std::memcpy (&q_, q, sizeof (mpq_t));
    clear_p_ = true;
    return *this;
}

void Ratio::move_to (mpq_t q)
{
    assert (clear_p_);
    mpq_clear (q);
    init_move_to (q);
}

void Ratio::init_move_to (mpq_t q)
{
    assert (clear_p_);
    std::memcpy (q, &q_, sizeof (mpq_t));
    clear_p_ = false;
}

void Ratio::copy_to (mpq_t q) const
{
    mpq_set (q, q_);
}

void Ratio::init_copy_to (mpq_t q) const
{
    mpq_init (q);
    copy_to (q);
}


Ratio operator "" _Ratio (const char *str, size_t)
{
    return Ratio::from_str (str);
}

Ratio operator "" _Ratio (const char *str)
{
    return Ratio::from_str (str);
}

Ratio operator "" _Q (const char *str, size_t)
{
    return Ratio::from_str (str);
}

Ratio operator "" _Q (const char *str)
{
    return Ratio::from_str (str);
}

Ratio::Ratio (const Ratio& q)
{
    init();
    mpq_set (q_, q.q_);
}

Ratio::Ratio (Ratio&& q)
{
    std::memcpy (this, &q, sizeof (q));
    q.clear_p_ = false;
}

Ratio& Ratio::operator = (const Ratio& q)
{
    assert (this != &q);
    assert (clear_p_);
    mpq_set (q_, q.q_);
    return *this;
}

Ratio& Ratio::operator = (Ratio&& q)
{
    assert (this != &q);
    clear();
    std::memcpy (this, &q, sizeof (q));
    q.clear_p_ = false;
    return *this;
}

Int Ratio::num () const
{
    Int z;
    return z.copy (numref());
}

Int Ratio::den () const
{
    Int z;
    return z.copy (denref());
}

int Ratio::cmp (const Int& z) const { return mpq_cmp_z (q_, z.z_); }
bool Ratio::operator == (const Int& z) const { return cmp (z) == 0; }
bool Ratio::operator != (const Int& z) const { return cmp (z) != 0; }
bool Ratio::operator <  (const Int& z) const { return cmp (z) <  0; }
bool Ratio::operator <= (const Int& z) const { return cmp (z) <= 0; }
bool Ratio::operator >  (const Int& z) const { return cmp (z) >  0; }
bool Ratio::operator >= (const Int& z) const { return cmp (z) >= 0; }


Ratio::operator double() const
{
    return mpq_get_d (q_);
}

Ratio::operator Float() const
{
    return Float (*this);
}

char* Ratio::to_str (int base) const
{
    // By default GMP uses malloc, realloc and free.
    return mpq_get_str (nullptr /* alloc */, base, q_); 
}

Ratio Ratio::from_str (const char *str, int base)
{
    Ratio q; // = 0 / 1

    char upper = '0' + (base ? base : 10);
    if (str[0] >= '0' && ! str[1] && str[0] < upper)
        mpz_set_ui (q.numref(), (unsigned long) (str[0] - '0'));
    else if (str[0] == '-' && str[1] >= '0' && ! str[2] && str[1] < upper)
        mpz_set_si (q.numref(), (long) ('0' - str[1]));
    else if (! str[0])
        { /* 0 */ }
    else
        q.set_str (str, base);

    return q;
}

Ratio& Ratio::set_str (const char *str, int base)
{
    mpq_set_str (q_, str, base);
    mpq_canonicalize (q_);
    return *this;
}

void Ratio::print (FILE *stream, int base) const
{
    // By default GMP uses malloc, realloc and free.
    char *str = to_str (base); 
    fputs (str, stream);
    free (str);
}


Ratio& Ratio::operator ++ () // pre
{
    mpz_add (numref(), numref(), denref());
    if (0 == mpz_sgn (numref()))
        mpz_set_ui (denref(), 1ul);
    return *this;
}

Ratio& Ratio::operator -- () // pre
{
    mpz_sub (numref(), numref(), denref());
    if (0 == mpz_sgn (numref()))
        mpz_set_ui (denref(), 1ul);
    return *this;
}

Ratio Ratio::operator ++ (int) // post
{
    Ratio q { *this };
    ++(*this);
    return q;
}

Ratio Ratio::operator -- (int) // post
{
    Ratio q { *this };
    --(*this);
    return q;
}


Ratio Ratio::inv () const
{
    if (is_zero())
        error ("0.inv()");
    Ratio q;
    mpq_inv (q.q_, q_);
    return q;
}

Ratio Ratio::operator - () const
{
    Ratio w;
    mpq_neg (w.q_, q_);
    return w;
}

namespace
{
    // Just some fake mpq, mpz functions so we can use the macros below.

    void mpq__mymod (mpq_t m, const mpq_t a, const mpq_t b)
    {
        Ratio qa { Ratio().copy(a) };
        qa = qa.floor_mod (Ratio().copy(b));
        qa.move_to (m);
    }

    void mpz__my_rounddiv_q (mpz_t q, const mpz_t a, const mpz_t b)
    {
        Int z { Int().copy (a) };
        z = z.round_div (Int().copy (b));
        z.move_to (q);
    }
}; // ::anon

#define MAKE_OP(OP, NAME, DIV0TXT)              \
    void Ratio::operator OP ##= (const Ratio& q)\
    {                                           \
        if (*DIV0TXT && q.is_zero())            \
            error ("%Qd%s", q_, "" DIV0TXT);    \
        mpq_## NAME (q_, q_, q.q_);             \
    }                                           \
    Ratio Ratio::operator OP (const Ratio& q) const \
    {                                           \
        if (*DIV0TXT && q.is_zero())            \
            error ("%Qd%s", q_, "" DIV0TXT);    \
        Ratio w;                                \
        mpq_## NAME (w.q_, q_, q.q_);           \
        return w;                               \
    }                                           \
    void Ratio::operator OP ##= (const Int& z)  \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%Qd%s", q_, "" DIV0TXT);    \
        (*this) OP ##= Ratio (z);               \
    }                                           \
    Ratio Ratio::operator OP (const Int& z) const \
    {                                           \
        if (*DIV0TXT && z.is_zero())            \
            error ("%Qd%s", q_, "" DIV0TXT);    \
        return (*this) OP Ratio (z);            \
    }                                           \
    Ratio Int::operator OP (const Ratio& q) const \
    {                                           \
        if (*DIV0TXT && q.is_zero())            \
            error ("%Zd%s", z_, "" DIV0TXT);    \
        return Ratio (*this) OP q;              \
    }
    MAKE_OP (+, add, "")
    MAKE_OP (-, sub, "")
    MAKE_OP (*, mul, "")
    MAKE_OP (/, div, " / 0")
    MAKE_OP (%, _mymod, " % 0") // >= 0.
#undef MAKE_OP


#define MAKE_OP(NAME, OP)                               \
    Int Ratio::NAME () const                            \
    {                                                   \
        Int z;                                          \
        mpz_## OP (z.z_, mpq_numref (q_), mpq_denref (q_)); \
        return z;                                       \
    }                                                   \
    Int NAME (const Ratio& q)                           \
    {                                                   \
        return q.NAME ();                               \
    }                                                   \
    Ratio Ratio::NAME ##_mod (const Ratio& d) const     \
    {                                                   \
        if (d.is_zero())                                \
            error ("%Qd %s_mod 0", q_, #NAME);          \
        return (*this) - d * (*this / d).NAME();        \
    }
    MAKE_OP (trunc, tdiv_q)
    MAKE_OP (floor, fdiv_q)
    MAKE_OP (ceil, cdiv_q)
    MAKE_OP (round, _my_rounddiv_q)
#undef MAKE_OP

const mpz_t& Ratio::numref () const { return *(const mpz_t*) mpq_numref (q_); }
const mpz_t& Ratio::denref () const { return *(const mpz_t*) mpq_denref (q_); }
mpz_t& Ratio::numref () { return *(mpz_t*) mpq_numref (q_); }
mpz_t& Ratio::denref () { return *(mpz_t*) mpq_denref (q_); }

bool Ratio::is_perfect_square () const
{
    return (mpz_perfect_square_p (mpq_numref (q_))
            && mpz_perfect_square_p (mpq_denref (q_)));
}

Ratio Ratio::abs () const
{
    Ratio q (*this);
    mpz_abs (q.numref(), q.numref());
    return q;
}

Ratio Ratio::pow (long ex) const
{
    unsigned long uex = (unsigned long) ex;
    if (ex < 0)
    {
        if (is_zero())
            error ("0.pow(%ld)", ex);
        uex = -uex;
    }
    Int zz; mpz_pow_ui (zz.z_, numref(), uex);
    Int nn; mpz_pow_ui (nn.z_, denref(), uex);
    return ex < 0 ? Ratio (nn, zz) : Ratio (zz, nn);
}

Ratio Ratio::min (const Ratio& q) const { return (*this) < q ? (*this) : q; }
Ratio Ratio::max (const Ratio& q) const { return (*this) > q ? (*this) : q; }

Ratio min (const Ratio& q, const Ratio& w) { return q.min (w); }
Ratio max (const Ratio& q, const Ratio& w) { return q.max (w); }
Ratio gcd (const Ratio& q, const Ratio& w) { return q.gcd (w); }
Ratio lcm (const Ratio& q, const Ratio& w) { return q.lcm (w); }
Ratio abs (const Ratio& q) { return q.abs(); }

Ratio Ratio::operator << (mp_bitcnt_t shift) const
{
    Ratio w;
    mpq_mul_2exp (w.q_, q_, shift);
    return w;
    
}

Ratio Ratio::operator >> (mp_bitcnt_t shift) const
{
    Ratio w;
    mpq_div_2exp (w.q_, q_, shift);
    return w;
    
}

void Ratio::operator <<= (mp_bitcnt_t shift)
{
    *this = (*this) << shift;
}

void Ratio::operator >>= (mp_bitcnt_t shift)
{
    *this = (*this) >> shift;
}

Ratio Ratio::lcm (const Ratio& b) const
{
    Ratio ab { (*this * b).abs() };
    if (! ab.is_zero())
        ab /= gcd (b);
    return ab;
}

Ratio Ratio::gcd (const Ratio& b0) const
{
    Ratio a (*this);
    Ratio b (b0);

    while (!b.is_zero())
    {
        a = a.round_mod (b);
        mpq_swap (a.q_, b.q_);
    }
    return a.abs();
}

int cmp (const Ratio& q, const Ratio& w) { return q.cmp (w); }

int sgn (const Ratio& q) { return q.sgn(); }
bool is_perfect_square (const Ratio& q) { return q.is_perfect_square(); }
Ratio pow (const Ratio& q, long ex) { return q.pow (ex); }

Ratio Ratio::make (int i)
{
    Ratio r;
    mpq_set_si (r.q_, (long) i, 1L);
    return r;
}

Ratio Ratio::make (int i, int n)
{
    if (n == 0)
        error ("%d / 0", i);
    long zz = (long) (n < 0 ? -i : i);
    unsigned un = (unsigned) n;
    unsigned long nn = (unsigned long) (n < 0 ? -un : un);
    Ratio r;
    mpq_set_si (r.q_, zz, nn);
    return r;
}

std::ostream& operator << (std::ostream& out, const Ratio& q)
{
    q.print (stdout, 10);
    return out;
}


const Ratio Ratio::Zero;
const Ratio Ratio::One (Ratio::make(1));
const Ratio Ratio::Two (Ratio::make(2));
const Ratio Ratio::minusOne (Ratio::make(-1));
const Ratio Ratio::Half (Ratio::make(1,2));
