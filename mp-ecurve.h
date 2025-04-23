#ifndef MP_ECURVE_H
#define MP_ECURVE_H

#include <iostream>

#include "mp-mint.h"

// Elliptic Curve mod Int and point on it.

class ECPoint;

// An Elliptic Curve E(a,b): y^2 = x^3 + ax + b  over  Z / n*Z

class ECurve
{
private:
    friend ECPoint;
    MInt a_, b_;
    mutable Int divisor_ { Int::One };
    mutable bool has_divisor_ = false;
    ECurve() {}
public:
    ECurve (const MInt&, const MInt&);
    MInt& a() { return a_; }
    MInt& b() { return b_; }
    const MInt& a() const { return a_; }
    const MInt& b() const { return b_; }
    Int divisor () const { return divisor_; }
    bool has_divisor () const { return has_divisor_; }
    void set_divisor (const Int&) const;
    MInt disc () const;
    MInt::Modulus modulus() const;
    void print (bool with_modulus = true) const;
    bool operator == (const ECurve&) const;
    bool operator != (const ECurve&) const;
};

std::ostream& operator << (std::ostream&, const ECurve&);

class ECPoint
{
public:
    typedef std::shared_ptr<const ECurve> Curve;
private:
    Curve E_ = nullptr;
    MInt x_, y_;
    bool zero_p_ = true;
    ECPoint () {};
public:
    ECPoint (Curve);
    ECPoint (Curve, const MInt&, const MInt&);

    const MInt& x() const { return x_; }
    const MInt& y() const { return y_; }
    Curve E() const { return E_; }
    bool is_zero() const { return zero_p_; }

    ECPoint operator - () const;
    ECPoint operator + (const ECPoint&) const;
    ECPoint operator - (const ECPoint& Q) const { return (*this) + (-Q); }
    ECPoint operator * (Int) const;
    void operator *= (const Int& k) { (*this) = (*this) * k; }
    void operator += (const ECPoint& Q) { (*this) = (*this) + Q; }
    void operator -= (const ECPoint& Q) { (*this) = (*this) + (-Q); }
    ECPoint add (const ECPoint&, Int&) const;
    void print (bool with_modulus = false) const;

    static ECPoint  random (MInt::Modulus m);
    static Curve make_curve (const MInt& a, const MInt& b);
    static Curve make_curve (const ECurve&);
};

ECPoint operator * (const Int&, const ECPoint&);

std::ostream& operator << (std::ostream&, const ECPoint&);

#endif // MP_ECURVE_H
