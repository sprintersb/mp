#ifndef MP_RATIO_H
#define MP_RATIO_H

/* From gmp.h:

   The prototypes for gmp_vprintf etc are provided only if va_list is defined,
   via an application having included <stdarg.h>.  Usually va_list is a typedef
   so can't be tested directly, but C99 specifies that va_start is a macro.

   A similar mechanism is used in <mpfr.h>.  Thus, include <cstdarg> now. */

#include <cstdarg>

#include <gmp.h>

#include <cstdio>
#include <iosfwd> // Just Forward decls for iostream.

class Int;
class Float;
template<class A> class Poly;

class Ratio
{
    mpq_t q_;
    bool clear_p_ = false;
    void clear();
    void init();
public:
    ~Ratio ();
    const mpq_t& mpq () const { return q_; }
    mpq_t& mpq() { return q_; }
    // In order to initialize Ratio from mpq_t use
    // Ratio().copy(const mpq_t) or Ratio().move(mpq_t).
    Ratio ();
    //Ratio (int);
    Ratio (const Int&);
    Ratio (const Int&, const Int&);
    Ratio (const Ratio&);
    Ratio (Ratio&&);
    Ratio& operator = (const Ratio&);
    Ratio& operator = (Ratio&&);
    Ratio& move (mpq_t);
    Ratio& copy (const mpq_t);
    // Use these if the target had been initialized.
    void move_to (mpq_t);
    void copy_to (mpq_t) const;
    // Use these if the target is not initialized yet.
    void init_move_to (mpq_t);
    void init_copy_to (mpq_t) const;
    const mpz_t& numref () const;
    const mpz_t& denref () const;
    mpz_t& numref ();
    mpz_t& denref ();

    Int num () const; // Numerator
    Int den () const; // Denominator

    static const Ratio Two;
    static const Ratio One;
    static const Ratio Half;
    static const Ratio Zero;
    static const Ratio minusOne;

    static Ratio from_str (const char*, int base = 0);
    Ratio& set_str (const char*, int base = 0);
    char* to_str (int base=10) const; // Uses malloc for string;
    void print (FILE*, int base=10) const;

    Ratio inv () const;
    Ratio operator - () const;

#define MAKE_OP(OP)                         \
    void operator OP ##= (const Ratio&);    \
    void operator OP ##= (const Int&);      \
    Ratio operator OP (const Ratio&) const; \
    Ratio operator OP (const Int&) const;
    MAKE_OP (+)
    MAKE_OP (-)
    MAKE_OP (*)
    MAKE_OP (/)
    MAKE_OP (%)
#undef MAKE_OP

    int cmp (const Ratio& q) const { return mpq_cmp (q_, q.q_); }
    bool operator == (const Ratio& q) const { return mpq_equal (q_, q.q_); }
    bool operator != (const Ratio& q) const { return ! (q == *this); }
    bool operator <  (const Ratio& q) const { return cmp (q) <  0; }
    bool operator <= (const Ratio& q) const { return cmp (q) <= 0; }
    bool operator >  (const Ratio& q) const { return cmp (q) >  0; }
    bool operator >= (const Ratio& q) const { return cmp (q) >= 0; }
    int cmp (const Int&) const;
    bool operator == (const Int&) const;
    bool operator != (const Int&) const;
    bool operator <  (const Int&) const;
    bool operator <= (const Int&) const;
    bool operator >  (const Int&) const;
    bool operator >= (const Int&) const;

    Ratio operator << (mp_bitcnt_t) const;
    Ratio operator >> (mp_bitcnt_t) const;
    void operator <<= (mp_bitcnt_t);
    void operator >>= (mp_bitcnt_t);
    Ratio& operator ++ ();
    Ratio& operator -- ();
    Ratio operator ++ (int);
    Ratio operator -- (int);
    
    int sgn () const { return mpq_sgn (q_); }
    Ratio abs () const;
    //Ratio sqrt () const;
    bool is_integer() const { return 0 == mpz_cmp_ui (denref(), 1ul); }
    bool is_natural() const { return sgn() > 0 && is_integer(); }
    bool is_natural0() const { return sgn() >= 0 && is_integer(); }
    bool is_perfect_square () const;
    bool is_zero () const { return 0 == sgn(); }
    bool is_unit () const { return is_one() || 0 == mpq_cmp_si (q_, -1L, 0UL); }
    bool is_one () const  { return 0 == mpq_cmp_ui (q_, 1UL, 0UL); }
    Ratio pow (long) const;
    Int trunc () const; // Round to zero.
    Int round () const; // Round to nearest.
    Int floor () const; // Round down.
    Int ceil () const;  // Round up.
    Ratio trunc_mod (const Ratio&) const;
    Ratio round_mod (const Ratio&) const;
    Ratio floor_mod (const Ratio&) const;
    Ratio ceil_mod (const Ratio&) const;
    Ratio gcd (const Ratio&) const;
    Ratio lcm (const Ratio&) const;
    Ratio min (const Ratio&) const;
    Ratio max (const Ratio&) const;
    Ratio mod (const Ratio&) const; // Non-negative.
    explicit operator double() const;
    explicit operator Float() const;

    void print() const { gmp_printf ("%Qd", q_); }

    Poly<Ratio> operator + (const Poly<Ratio>&) const;
    Poly<Ratio> operator - (const Poly<Ratio>&) const;
    Poly<Ratio> operator * (const Poly<Ratio>&) const;

    static Ratio make (int);
    static Ratio make (int, int);
};


Ratio operator "" _Ratio (const char[], size_t);
Ratio operator "" _Ratio (const char*);
Ratio operator "" _Q (const char[], size_t);
Ratio operator "" _Q (const char*);

/*Ratio operator - (long, const Ratio&);
Ratio operator + (long, const Ratio&);
Ratio operator * (long, const Ratio&);
Ratio operator / (long, const Ratio&);
Ratio operator % (long, const Ratio&);

bool operator == (long, const Ratio&);
bool operator != (long, const Ratio&);
bool operator <  (long, const Ratio&);
bool operator <= (long, const Ratio&);
bool operator >  (long, const Ratio&);
bool operator >= (long, const Ratio&);*/

//int cmp (const Ratio&, const Ratio&);
int sgn (const Ratio&);
Ratio abs (const Ratio&);
Ratio pow (const Ratio&, long);
Ratio sqrt (const Ratio&);
bool perfect_square_p (const Ratio&);
Ratio gcd (const Ratio&, const Ratio&);
Ratio lcm (const Ratio&, const Ratio&);
Ratio min (const Ratio&, const Ratio&);
Ratio max (const Ratio&, const Ratio&);
Ratio mod (const Ratio&, const Ratio&); // Always non-negative.

Int trunc (const Ratio&);
Int round (const Ratio&);
Int floor (const Ratio&);
Int ceil (const Ratio&);

std::ostream& operator << (std::ostream&, const Ratio&);

#endif // MP_RATIO_H
