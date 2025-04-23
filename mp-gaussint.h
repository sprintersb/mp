#ifndef MP_GAUSSINT_H
#define MP_GAUSSINT_H

#include "mp-int.h"
#include <iostream>
#include <set>
#include <vector>
#include <cstdio>

template<class X> class Factors;
template<class X> class Power;

class GaussInt
{
    typedef GaussInt Self;
    Int x_, y_;

public:
    GaussInt () {}
    GaussInt (const Int& x) : x_(x) {}
    GaussInt (const Int& x, const Int& y) : x_(x), y_(y) {}
    Int& x() { return x_; }
    Int& y() { return y_; }
    const Int& x() const { return x_; }
    const Int& y() const { return y_; }
    Int& real() { return x_; }
    Int& imag() { return y_; }
    const Int& real() const { return x_; }
    const Int& imag() const { return y_; }
    bool operator == (const GaussInt& a) const
    {
        return x_ == a.x_ && y_ == a.y_;
    }
    bool operator != (const GaussInt& a) const { return ! ((*this) == a); }
    bool operator == (const Int& a) const { return x_ == a && y_.is_zero(); }
    bool operator != (const Int& a) const { return ! ((*this) == a); }
    bool is_real() const { return y_.is_zero(); }
    bool is_imag() const { return x_.is_zero(); }
    bool is_complex() const { return !is_real() && !is_imag(); }
    bool is_natural() const { return is_real() && x_.sgn() >= 0; }
    Self operator ~ () const { return Self (x_ , -y_); }
    Self operator - () const { return Self (-x_ , -y_); }
    Self operator - (const Self& a) const { return Self (x_-a.x_, y_-a.y_); }
    Self operator + (const Self& a) const { return Self (x_+a.x_, y_+a.y_); }
    Self operator - (const Int& a) const { return Self (x_ - a, y_); }
    Self operator + (const Int& a) const { return Self (x_ + a, y_); }
    Self operator - (long i) const { return Self (x_ - i, y_); }
    Self operator + (long i) const { return Self (x_ + i, y_); }
    Self operator * (const Int& a) const { return Self (x_ * a, y_ * a); }
    Self operator * (long i) const { return Self (x_ * i, y_ * i); }
    Self operator * (const Self& a) const
    {
        return Self (x_*a.x_ - y_*a.y_, x_*a.y_ + y_*a.x_);
    }
    void operator *= (const Self& a) { (*this) = (*this) * a; }
    Self operator % (const Self&) const;
    Self operator / (const Self&) const;
    Self operator / (const Int&) const;
    Self operator / (long) const;
    void operator /= (const Self&);
    void operator /= (const Int&);
    void operator /= (long);
    Self round_div (const Self&) const;
    Int sprod (const Self& a) const { return x_*a.x_ + y_*a.y_; }
    Self pow (int) const;
    enum Canonical { mod_void, mod_units, mod_negation, mod_units_conjugation,
                     mod_conjugation };
    Self canonical (Canonical) const;
    int cmp (const Self&, Canonical) const;
    Int norm () const { return sprod (*this); }
    bool is_one () const { return x_.is_one() && y_.is_zero(); }
    bool is_unit () const { return ((x_.is_zero() && y_.is_unit())
                                    || (y_.is_zero() && x_.is_unit())); }
    bool is_zero () const { return x_.is_zero() && y_.is_zero(); }
    int is_probab_prime (int tries) const
    {
        Int n = is_real() ? x_.abs() : is_imag() ? y_.abs() : norm();
        return n <= 1L ? 0 : n.is_probab_prime (tries);
    }
    Self next_prime() const;
    Self divmod (const Self&, Self&) const;
    Self gcd (const Self&) const;

    void print () const;
    char *to_str() const;
    
    static Self factor_integer_prime (const Int& p, int *p_report = nullptr);
    Factors<Self> factor (int verbose = 0) const;
    Self force_remove_factor (const Self&) const;

    struct ComparePythagorean { bool operator () (const Self&, const Self&) const; };
    struct CompareModuloUnits { bool operator () (const Self&, const Self&) const; };
    typedef std::set<Self, ComparePythagorean> PythagoreanTriples;

    static Factors<Self> factor (const Int&, int);
    static PythagoreanTriples solve_norm (const Factors<Int>&, bool primitive = false);
    static Self solve_norm_integral_part (const Factors<Int>&,
                                          std::vector<Power<Self>>&, bool primitive);
    static const GaussInt Zero;
    static const GaussInt One;
    static const GaussInt I;
    struct Stat;
};

/*GaussInt operator + (long, const GaussInt&);
GaussInt operator - (long, const GaussInt&);
GaussInt operator * (long, const GaussInt&);*/

std::ostream& operator << (std::ostream&, const GaussInt&);


struct GaussInt::Stat
{
    const char *percent = "\\%";
    int n_prime_a_plus_b = 0;
    int n_total = 0;
    void print (FILE*, const char *fmt) const;
    void operator += (const Stat&);
};

#endif // MP_GAUSSINT_H
