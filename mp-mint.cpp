#include "mp-mint.h"

#include <cstring>
#include <cassert>
#include <cstdlib>
#include <cstdio>

#define DIAGNOSTIC_NO_FORMAT_CHECK
#include "diagnostic.h"

const MInt::Modulus MInt::dummy_modulus_ { make_modulus (1) };

MInt::MInt() : m_(MInt::dummy_modulus_), z_(0L) {}

MInt::MInt (Modulus m, const Int& z) : m_(m), z_(z.mod(*m)) {}

void MInt::print (FILE *stream, int base) const
{
    // By default GMP uses malloc, realloc and free.
    z_.print (stream, base); 
}

auto MInt::make_modulus (const Int& m) -> Modulus
{
    if (m < 1)
        error ("muduli must be positive, got %Zd", m.mpz());
    return std::make_shared<Int> (m);
}

bool MInt::same_modulus (const MInt& a) const
{
    return m_ == a.m_ || *m_ == *a.m_;
}

#include "diagnostic.h"

#define M_SAME(M1, M2)                                  \
    do {                                                \
        if (MInt::check_modulus                         \
            && M1 != M2                                 \
            && (*M1) != (*M2))                          \
        {                                               \
            const char *_s = __PRETTY_FUNCTION__;       \
            error_at (__FILE__, __LINE__, "%s: modulus %Zd != %Zd",\
                      _s, (*M1).mpz(), (*M2).mpz());    \
        }                                               \
    } while(0)

MInt MInt::operator - () const
{
    return MInt (m_, z_.is_zero() ? z_ : *m_ - z_);
}

#define MAKE_OP(OP)                                 \
    void MInt::operator OP ##= (const MInt& b)      \
    {                                               \
        M_SAME (m_, b.m_);                          \
        z_ = (z_ OP b.z_).mod (*m_);                \
    }                                               \
    void MInt::operator OP ##= (const Int& i)       \
    {                                               \
        z_ = (z_ OP i.mod (*m_)).mod (*m_);         \
    }                                               \
    void MInt::operator OP ##= (long i)             \
    {                                               \
        z_ = (z_ OP i).mod (*m_);                   \
    }                                               \
    MInt MInt::operator OP (const MInt& b) const    \
    {                                               \
        M_SAME (m_, b.m_);                          \
        return MInt { m_, z_ OP b.z_ };             \
    }                                               \
    MInt MInt::operator OP (const Int& i) const     \
    {                                               \
        return MInt { m_, z_ OP i.mod (*m_) };      \
    }                                               \
    MInt MInt::operator OP (long i) const           \
    {                                               \
        return MInt { m_, z_ OP i };                \
    }                                               \
    MInt operator OP (long i, const MInt& z)        \
    {                                               \
        return MInt { z.m_, i } OP z;               \
    }
    MAKE_OP (+)
    MAKE_OP (-)
    MAKE_OP (*)
#undef MAKE_OP

// Invert mod m_, return divisor if not invertible.
MInt MInt::inv (Int& d) const
{
    Int s;
    mpz_gcdext (d.z_, s.z_, nullptr, z_.z_, m_->z_);
    return MInt { m_, s };
}

bool MInt::is_zero (const Int& b) const
{
    static Int xzero (0);
    return mpz_congruent_p (xzero.z_, b.z_, m_->z_);
}

bool MInt::is_zero () const
{
    return is_zero (z_);
}

MInt MInt::div (const Int& b, Int& d) const
{
    if (is_zero (b))
        error ("division by zero mod %Zd", m_->mpz());

    return div (MInt{m_,b}, d);
}

MInt MInt::div (const MInt& b, Int& d) const
{
    M_SAME (m_, b.m_);
    if (b.is_zero ())
        error ("division by zero mod %Zd", m_->mpz());

    return (*this) * b.inv (d);
}

MInt MInt::operator / (const MInt& b) const
{
    M_SAME (m_, b.m_);
    if (b.is_zero ())
        error ("division by zero mod %Zd", m_->mpz());

    Int d;
    MInt q { div (b, d) };
    if (d != 1)
        error ("%Zd has divisor %Zd", m_->mpz(), d.mpz());
    return q;
}

MInt MInt::operator / (const Int& b) const
{
    if (is_zero (b))
        error ("division by zero mod %Zd", m_->mpz());

    Int d;
    MInt q { div (b, d) };
    if (d != 1)
        error ("%Zd has divisor %Zd", m_->mpz(), d.mpz());
    return q;
}

void MInt::operator /= (const MInt& b)
{
    M_SAME (m_, b.m_);
    if (b.is_zero())
        error ("division by zero mod %Zd", m_->mpz());

    (*this) = (*this) / b;
}

void MInt::operator /= (const Int& b)
{
    if (is_zero (b))
        error ("division by zero mod %Zd", m_->mpz());

    (*this) = (*this) / b;
}


#define MAKE_OP(OP)                                 \
    bool MInt::operator OP (const MInt& b) const    \
    {                                               \
        M_SAME (m_, b.m_);                          \
        return z_ OP b.z_;                          \
    }                                               \
    bool MInt::operator OP (const Int& i) const     \
    {                                               \
        return z_ OP i.mod (*m_);                   \
    }                                               \
    bool operator OP (const Int& i, const MInt& b)  \
    {                                               \
        return i.mod (*b.m_) OP b.z_;               \
    }
    MAKE_OP (==)
    MAKE_OP (!=)
#undef MAKE_OP

MInt MInt::pow (const Int& ex) const
{
    MInt b { m_, 0 };
    mpz_powm (b.z_.z_, z_.z_, ex.z_, m_->z_);
    return b;
}

MInt pow (const MInt& b, const Int& ex)
{
    return b.pow (ex);
}


MInt MInt::random (Modulus m)
{
    return MInt { m, Int::random (*m) };
}


std::ostream& operator << (std::ostream& ost, const MInt& z)
{
    z.print (stdout, 10);
    return ost;
}


int MInt::legendre (MInt::SquareKind *p_report) const
{
    MInt::SquareKind dummy;
    MInt::SquareKind& kind = p_report ? *p_report : dummy;

    if (m_->is_even())
    {
        kind = *m_ == 2 ? MInt::Square : MInt::NotPrime;
        return 1 & (long) z_;
    }

    int leg = 0;

    if (*this == 0)
    {
        kind = MInt::Square;
    }
    else
    {
        MInt q = this->pow (*m_ >> 1);

        if (q != 1 && q != -1)
            kind = MInt::NotPrime;
        else
        {
            kind = q == 1 ? MInt::Square : MInt::NonSquare;
            leg = kind == MInt::Square ? 1 : -1;
        }
    }

    return leg;
}

static MInt sequence_lucas (const Int& n, const MInt& P, const MInt& Q);

MInt MInt::sqrt (MInt::SquareKind *p_report) const
{
    MInt::SquareKind dummy;
    MInt::SquareKind& kind = p_report ? *p_report : dummy;
    const MInt& self = *this;

    int leg = legendre (&kind);

    if (kind == MInt::NotPrime)
        goto _not_prime;

    if (kind == MInt::NonSquare)
        return MInt();

    if (*m_ == 2)
        return self;

    assert (m_->is_odd());

    if (m_->bit(1) == 1)
    {
        // 3 mod 4.

        MInt q = self.pow ((*m_ + 1) / 4);
        if (q * q == self)
        {
            kind = MInt::Square;
            return q;
        }

        goto _not_prime;
    }

    if (m_->bit(1) == 0)
    {
        // 1 mod 4.

        MInt r (m_, 0);
        do {
            r += 1;
            leg = (r*r - 4*self).legendre (&kind);
        } while (kind == MInt::Square && !r.is_zero());

        if (kind == MInt::NotPrime || r.is_zero())
            goto _not_prime;

        assert (leg == -1);

        Int d;
        MInt q = sequence_lucas (*m_ >> 2, r, self);
        q *= self.div (2 * r, d);
        if (d > 1)
            goto _not_prime;
        
        assert (q*q == self);
        return q;
    }

_not_prime:;
    kind = MInt::NotPrime;
    return MInt();
}


MInt MInt::find_nonsquare (MInt::Modulus m, MInt::SquareKind *p_report)
{
    MInt::SquareKind dummy;
    MInt::SquareKind& kind = p_report ? *p_report : dummy;

    if (m->is_even() || *m == 1)
    {
        kind = *m == 2 ? MInt::Square : MInt::NotPrime;
        return MInt (m, 0);
    }

    for (MInt i (m, 2); i.z() < *m; i += 1)
    {
        (void) i.legendre (&kind);
        if (kind != MInt::Square)
            return i;
    }

    kind = MInt::NotPrime;
    return MInt();
}



MInt sequence_lucas (const Int& n, const MInt& P, const MInt& Q)
{
    assert (n >= 1);

    MInt::Modulus m = P.modulus();
    MInt two (m, 2);
    MInt V0 = two;
    MInt V1 = P;
    MInt V2 = P * V1 - Q * V0;

    Int d; // Potential factor of M.
    // Wn = V_{2n} / Q^n
    // W1 = V2 / Q
    MInt W1 = V2.div (Q, d);
    // W2 = W1^2 - 2  (wg. W_{2n} = W_n^2 - 2)
    MInt W2 = W1 * W1 - two;

    // (Aj, Bj)  with  (A0, B0) = (W1, W2)
    MInt Aj = W1, Bj = W2;
    MInt Aj1 = V0, Bj1 = V0; // Just to satisfy some constructor.

    for (int bitno = n.msbit() - 1; bitno >= 0; --bitno)
    {
        if (n.bit (bitno) == 0)
        {
            Aj1 = Aj * Aj - two;
            Bj1 = Aj * Bj - W1;
        }
        else
        {
            Aj1 = Aj * Bj - W1;
            Bj1 = Bj * Bj - two;
        }

        Aj = Aj1;
        Bj = Bj1;
        //cout << bitno << ": (" << Aj << ", " << Bj << ")" << endl;
    }
    // Assert (At, Bt) = (Wn, W{n+1})  mit t = msb 
    return Aj + Bj;
}


bool MInt::check_modulus = true;
