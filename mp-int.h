#ifndef MP_INT_H
#define MP_INT_H

/* From gmp.h:

   The prototypes for gmp_vprintf etc are provided only if va_list is defined,
   via an application having included <stdarg.h>.  Usually va_list is a typedef
   so can't be tested directly, but C99 specifies that va_start is a macro.

   A similar mechanism is used in <mpfr.h>.  Thus, include <cstdarg> now. */

#include <cstdarg>

#include <gmp.h>

#include <cstdio>
#include <iosfwd> // Just Forward decls for iostream.

class MInt;
class Ratio;
class GaussInt;
class Float;

template<class A> class Poly;
template<class X> class Factors;


class Int
{
    friend MInt;
    friend Ratio;
    mpz_t z_;
    bool clear_p_ = false;
    void clear();
public:
    ~Int () { clear(); }
    Int (long z = 0);
    const mpz_t& mpz() const { return z_; }
    mpz_t& mpz() { return z_; }
    // In order to initialize Int from mpz_t use
    // Int().copy(const mpz_t) or Int().move(mpz_t).
    Int (const Int&);
    Int (Int&&);
    Int& operator = (const Int&);
    Int& operator = (Int&&);
    Int& move (mpz_t);
    Int& copy (const mpz_t);
    // Use these if the target had been initialized.
    void move_to (mpz_t);
    void copy_to (mpz_t) const;
    // Use these if the target is not initialized yet.
    void init_move_to (mpz_t);
    void init_copy_to (mpz_t) const;
    static const Int Zero;
    static const Int One;
    static const Int Two;
    static Int make (int);
    explicit operator long() const;
    explicit operator double() const;
    explicit operator Float() const; // mp-float.cpp
    static Int from_Float (const Float&, int /*mpfr_rnd_t*/); // mp-float.cpp

    static Int from_str (const char*, int base = 0);
    Int& set_str (const char*, int base = 0);
    char* to_str (int base=10) const; // Uses malloc for string;
    void print (FILE*, int base=10) const;

    Int operator - () const;
    Int operator ~ () const;
    Int operator - (const Int&) const;
    Int operator + (const Int&) const;
    Int operator * (const Int&) const;
    Int operator / (const Int&) const;
    Int operator % (const Int&) const;
    Int operator & (const Int&) const;
    Int operator | (const Int&) const;
    Int operator ^ (const Int&) const;
    void operator += (const Int&);
    void operator -= (const Int&);
    void operator *= (const Int&);
    void operator /= (const Int&);
    void operator %= (const Int&);
    void operator &= (const Int&);
    void operator |= (const Int&);
    void operator ^= (const Int&);

    Ratio operator - (const Ratio&) const;
    Ratio operator + (const Ratio&) const;
    Ratio operator * (const Ratio&) const;
    Ratio operator / (const Ratio&) const;
    Ratio operator % (const Ratio&) const;

    bool operator == (const Int&) const;
    bool operator != (const Int&) const;
    bool operator <  (const Int&) const;
    bool operator <= (const Int&) const;
    bool operator >  (const Int&) const;
    bool operator >= (const Int&) const;

    bool operator == (long) const;
    bool operator != (long) const;
    bool operator <  (long) const;
    bool operator <= (long) const;
    bool operator >  (long) const;
    bool operator >= (long) const;
    int abscmp (const Int& z) const { return mpz_cmpabs (z_, z.z_); }
    int abscmp (long i) const
    {
        return mpz_cmpabs_ui (z_, i < 0 ? -(unsigned long) i : (unsigned long) i);
    }
    Int operator << (mp_bitcnt_t) const;
    Int operator >> (mp_bitcnt_t) const;
    void operator <<= (mp_bitcnt_t);
    void operator >>= (mp_bitcnt_t);
    Int& operator ++ ();
    Int& operator -- ();
    Int operator ++ (int);
    Int operator -- (int);

    int cmp (const Int& z) const { return mpz_cmp (z_, z.z_); }
    int sgn () const { return mpz_sgn (z_); }
    Int abs () const;
    Int sqrt () const;
    Int next_prime () const;
    // prime? 0=No, 1=Likely, 2=Yes.
    int is_probab_prime (int tries) const { return mpz_probab_prime_p (z_, tries); }
    bool is_perfect_square () const { return mpz_perfect_square_p (z_); }
    bool is_odd () const { return mpz_odd_p (z_); }
    bool is_even () const { return mpz_even_p (z_); }
    bool is_zero () const { return 0 == sgn(); }
    bool is_unit () const { return 0 == mpz_cmpabs_ui (z_, 1UL); }
    bool is_one () const  { return 0 == mpz_cmp_ui (z_, 1UL); }
    mp_bitcnt_t popcount () const { return mpz_popcount (z_); }
    long ilog2 () const;
    long exact_log2 () const;
    long bitsize () const; // 0 -> 1; negative -> -1.
    Int pow (unsigned long) const;
    Int divmod (const Int&, Int&) const; // with floor div.
    Int round_div (const Int&) const;
    Int trunc_div (const Int&) const;
    Int floor_div (const Int&) const;
    Int ceil_div (const Int&) const;
    Int phi () const; // mp-factor.cpp
    Int gcd (const Int&) const;
    Int lcm (const Int&) const;
    Int min (const Int&) const;
    Int max (const Int&) const;
    Int mod (const Int&) const; // Non-negative.
    int bit (mp_bitcnt_t bitno) const { return mpz_tstbit (z_, bitno); }
    Int setbit (int bitno, bool val = 1) const;
    long msbit () const { return ilog2(); }
    static Int factorial (unsigned long, int = 1);
    static Int binom (long, long);

    Int remove_factor (const Int&, int* = nullptr) const;
    void print() const { gmp_printf ("%Zd", z_); }

    static Int random (const Int&);
    GaussInt operator + (const GaussInt&) const;
    GaussInt operator - (const GaussInt&) const;
    GaussInt operator * (const GaussInt&) const;

    Poly<Int> operator + (const Poly<Int>&) const;
    Poly<Int> operator - (const Poly<Int>&) const;
    Poly<Int> operator * (const Poly<Int>&) const;
};

Int operator "" _Int (const char[], size_t);
Int operator "" _Int (const char*);
Int operator "" _Z (const char[], size_t);
Int operator "" _Z (const char*);

Int operator - (long, const Int&);
Int operator + (long, const Int&);
Int operator * (long, const Int&);
Int operator / (long, const Int&);
Int operator % (long, const Int&);
Int operator & (long, const Int&);
Int operator | (long, const Int&);
Int operator ^ (long, const Int&);

bool operator == (long, const Int&);
bool operator != (long, const Int&);
bool operator <  (long, const Int&);
bool operator <= (long, const Int&);
bool operator >  (long, const Int&);
bool operator >= (long, const Int&);

int cmp (const Int&, const Int&);
int sgn (const Int&);
Int abs (const Int&);
Int pow (const Int&, unsigned long);
Int sqrt (const Int&);
Int next_prime (const Int&);
bool is_perfect_square (const Int&);
Int phi (const Int&); // mp-factor.cpp
Int gcd (const Int&, const Int&);
Int lcm (const Int&, const Int&);
Int min (const Int&, const Int&);
Int max (const Int&, const Int&);
Int mod (const Int&, const Int&); // Always non-negative.
int is_probab_prime (const Int&, int); // 0: no, 1: pseudo-prime, 2: prime.

Int trunc_div (const Int&, const Int&);
Int round_div (const Int&, const Int&);
Int floor_div (const Int&, const Int&);
Int ceil_div (const Int&, const Int&);

std::ostream& operator << (std::ostream&, const Int&);

#endif // MP_INT_H
