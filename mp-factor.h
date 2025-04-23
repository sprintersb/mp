#ifndef MP_FACTOR_H
#define MP_FACTOR_H

//#include "mp-int.h"

class Int;

#include <iostream>
#include <vector>
#include <set>

// Integer Factorization

// Holds 3 vector<Int>'s that represent:
// * prime factors
// * pseudo-prime factors
// * composite factors

template <class X>
class Power
{
public:
    typedef X base_type;
    class Sorter;
    X base_;
    mutable int expo_ = 0;
    Power () {}
    Power (const X& x) : Power(x,1) {}
    Power (const X& x, int e) : base_(x), expo_(e) {}
    void addto_expo (int n) const { expo_ += n; }
    X value() const { return base_.pow (expo_); }
    bool operator() (const Power&, const Power&) const; // Order for std::set.
};

template <class X>
std::ostream& operator << (std::ostream&, const Power<X>&);


template <class X>
class Factors
{
public:
    typedef std::vector<X> Xs;
    //typedef std::set<Int::Power, std::less<Int> > Powers;
    struct Sorter
    {
        bool operator() (const X& a, const X& b) const;
    };

    class Context;
    void set_powers();
    typedef std::set<Power<X>, Power<X>> Powers;
    enum Kind { Prime, PseudoPrime, Composite, Unit, Zero, Unknown };
    typedef int (*Callback) (const X&, Kind, Factors&);
    Xs primes_;
    Xs pseudo_primes_;
    Xs composites_;
    Powers powers_;
    Callback callback_ = nullptr;
    int stop_reason_ = 0; // != 0 : stopped by callback.
    X unit_, n_;
    int verbose_ = 0;
    int add (const X&, Kind = Unknown);
    int add_pseudo_prime (const X& d) { return add (d, PseudoPrime); }
    int add_prime (const X& d) { return add (d, Prime); }
    int add_composite (const X& d) { return add (d, Composite); }

    int is_fully_factored () const;
    X print (const char*, const Xs&) const;
    void print() const;
    void print_powers (bool with_unit = false, const char* end = "\n") const;
};

template <class X>
class Factors<X>::Context
{
public:
    Factors<X>::Callback callback_ = nullptr;
    const Factors<X> *known_ = nullptr;
    int verbose_ = 0;
};

#include "mp-int.h"
#include "mp-gaussint.h"

template <class X>
std::ostream& operator << (std::ostream&, const Factors<X>&);

extern template std::ostream& operator << (std::ostream&, const Factors<Int>&);
extern template std::ostream& operator << (std::ostream&, const Factors<GaussInt>&);

extern template class Factors<Int>;
extern template class Factors<GaussInt>;

Factors<Int> factor (const Int&, const Factors<Int>::Context* = nullptr);
Factors<GaussInt> factor (const GaussInt&, const Factors<GaussInt>::Context* = nullptr);

Factors<Int> factor_a_pow_n_minus_1 (const Int&, int, const Factors<Int>::Context* = nullptr);

Int remove_small_primes (const Int&, Factors<Int>&);
//Int find_factor_elliptic_curve (const Int&, int smooth);
Int phi (const Int&);
Int phi (const Factors<Int>&);
Int find_factor_Pollard_rho (const Int&, int n_loops, int verbose);
Int find_factor_Pollard_p_minus_1 (const Int&, int smooth, int verbose);
Int find_factor_ECM (const Int& n, int B1, int verbose, int *result);

int callback_stop_1_mod_4 (const Int&, Factors<Int>::Kind, Factors<Int>&);

#endif // MP_FACTOR_H
