#ifndef MP_MINT_H
#define MP_MINT_H

#include "mp-int.h"

#include <memory> // shared_ptr
#include <vector>

#include <cstdio>
#include <iostream>

class ECPoint;
class ECurve;

class MInt
{
    friend ECurve;
    friend ECPoint;
public:
    typedef std::shared_ptr<const Int> Modulus;
private:
    Modulus m_;
    Int z_;
    static const Modulus dummy_modulus_;
public:
    MInt();
    MInt (Modulus, const Int&);
    static Modulus make_modulus (const Int&);
    static bool check_modulus;
    //MInt (const MInt&) = delete;
    //MInt (MInt&&);
    //MInt& operator = (const MInt&);
    //MInt& operator = (MInt&&);
    void print (FILE*, int base) const;

    Int& z() { return z_; }
    const Int& z() const { return z_; }
    Modulus modulus() const { return m_; }
    bool same_modulus (const MInt&) const;
    MInt operator - () const;
    void operator /= (const MInt&);
    void operator /= (const Int&);
    MInt operator / (const MInt&) const;
    MInt operator / (const Int&) const;
    bool is_zero (const Int&) const;
    bool is_zero () const;
    MInt inv (Int&) const;
    MInt div (const MInt&, Int&) const;
    MInt div (const Int&, Int&) const;

#define _MAKE_OP(OP)                                    \
    void operator OP##= (const MInt&);                  \
    void operator OP##= (const Int&);                   \
    void operator OP##= (long);                         \
    MInt operator OP (const MInt&) const;               \
    MInt operator OP (const Int&) const;                \
    MInt operator OP (long) const;                      \
    friend MInt operator OP (const Int&, const MInt&);  \
    friend MInt operator OP (long, const MInt&);
    _MAKE_OP (+)
    _MAKE_OP (-)
    _MAKE_OP (*)
#undef _MAKE_OP

    bool operator == (const MInt&) const;
    bool operator == (const Int&) const;
    friend bool operator == (const Int&, const MInt&);
    bool operator != (const MInt&) const;
    bool operator != (const Int&) const;
    friend bool operator != (const Int&, const MInt&);

    MInt pow (const Int&) const;

    enum SquareKind { NonSquare, Square, NotPrime };

    static MInt find_nonsquare (Modulus, SquareKind* = nullptr);
    MInt sqrt (SquareKind* = nullptr) const;
    int legendre (MInt::SquareKind* = nullptr) const;

    static MInt random (Modulus);
    static MInt make (int);
};

MInt pow (const MInt&, const Int&);

std::ostream& operator << (std::ostream&, const MInt&);

#endif // MP_MINT_H
