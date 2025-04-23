#include "mp-ecurve.h"

#include "diagnostic.h"

#include <cstdlib>

#define M_SAME(M1, M2)                                  \
    do {                                                \
        if (M1 != M2 && (*M1) != (*M2))                 \
        {                                               \
            const char *_s = __PRETTY_FUNCTION__;       \
            const char *_sm1 = (*M1).to_str();          \
            const char *_sm2 = (*M2).to_str();          \
            error_at (__FILE__, __LINE__, "%s: modulus %s != %s", _s, _sm1, _sm2); \
        }                                               \
    } while(0)

ECurve::ECurve (const MInt& a, const MInt& b) : a_(a), b_(b)
{
    M_SAME (a.m_, b.m_);
}

MInt::Modulus ECurve::modulus() const
{
    return a_.m_;
}

MInt ECurve::disc () const
{
    return (-4) * a_.pow(3) - 27 * b_ * b_;
}

void ECurve::print (bool with_modulus) const
{
    char *a = a_.z_.to_str();
    char *b = b_.z_.to_str();
    char *m = with_modulus ? modulus()->to_str() : nullptr;
    printf ("E(a,b): a = %s; b = %s", a, b);
    if (m)
        printf (" (mod %s)", m);
    printf ("\n");
    free (a);
    free (b);
    free (m);
}

bool ECurve::operator == (const ECurve& E) const
{
    return (a_.same_modulus (E.a_)
            && b_.same_modulus (E.b_)
            && a_ == E.a_
            && b_ == E.b_);
}

bool ECurve::operator != (const ECurve& E) const
{
    return ! ((*this) == E);
}

void ECurve::set_divisor (const Int& d) const
{
    if (!has_divisor_
        && d > Int::One
        && d.mod (*modulus()) != Int::Zero)
    {
        if (*modulus() % d != Int::Zero)
            error ("%s does not divide %s", d.to_str(), modulus()->to_str());

        has_divisor_ = true;
        divisor_ = d;
        out ("E: divisor = %s | %s\n", d.to_str(), modulus()->to_str());
    }
}

std::ostream& operator << (std::ostream &ost, const ECurve& E)
{
    E.print();
    return ost;
}

ECPoint::ECPoint (Curve E)
    : E_(E), x_(MInt{E->modulus(),0}), y_(MInt{E->modulus(),0}), zero_p_(true)
{}


ECPoint::ECPoint (Curve E, const MInt& x, const MInt& y)
    : E_(E), x_(x), y_(y), zero_p_(false)
{
    M_SAME (x.m_, y.m_);
    M_SAME (E->modulus(), x.m_);
}


ECPoint ECPoint::random (MInt::Modulus m)
{
    out_context = true;
    MInt x, y, disc;
    ECurve E;
    ECPoint P;

    do {
        E.a_ = MInt::random (m);
        //out ("E.a = %s\n", E.a_.z_.to_str());
        x = MInt::random (m);
        y = MInt::random (m);
        E.b_ = y*y - x * (x*x + E.a_);
        P = ECPoint (ECPoint::make_curve (E), x, y);
    } while (P.is_zero() || E.disc() == 0);

    //std::cout << "x = " << x << " (mod " << *x.m_ << ")" << std::endl;
    return P;
}

void ECPoint::print (bool with_modulus) const
{
    char *x = x_.z_.to_str();
    char *y = y_.z_.to_str();
    char *m = with_modulus ? E_->modulus()->to_str() : nullptr;
    if (zero_p_)
        printf ("(x,y) = (0)");
    else
        printf ("(x,y) = (%s, %s)", x, y);
    if (m)
        printf ("; m = %s", m);
    free (x);
    free (y);
    free (m);
}


auto ECPoint::make_curve (const MInt& a, const MInt& b) -> Curve
{
    M_SAME (a.m_, b.m_);
    return std::make_shared<ECurve> (ECurve { a, b });
}

auto ECPoint::make_curve (const ECurve& E) -> Curve
{
    return std::make_shared<ECurve> (E);
}

ECPoint ECPoint::operator - () const
{
    return zero_p_ ? *this : ECPoint { E_, x_, -y_ };
}

ECPoint ECPoint::operator + (const ECPoint& Q) const
{
    if (E_->has_divisor())
        return ECPoint (E_);

    Int d;
    ECPoint R { add (Q, d) };

    if (d > Int::One)
    {
        E_->set_divisor (d);
        info (0, "E(a,b): P + Q -> %s divides %s", d.to_str(), E_->modulus()->to_str());
    }
    return R;
}

ECPoint ECPoint::add (const ECPoint& Q, Int& d) const
{
    d = 1;
    const ECPoint& P = *this;

    if (P.is_zero())    return Q;
    if (Q.is_zero())    return P;

    MInt s;

    if (P.x_ == Q.x_)
    {
        if (P.y_ == -Q.y_)
            // P + (-P) = 0.
            return ECPoint { E_ };

        // 2*P.
        s = (3 * P.x_*P.x_ + E_->a_).div (2 * P.y_, d);
    }
    else
    {
        s = (P.y_ - Q.y_).div (P.x_ - Q.x_, d);
    }

    MInt x = s*s - P.x_ - Q.x_;
    MInt y = s * (P.x_ - x) - P.y_;

    return ECPoint (E_, x, y);
}


ECPoint operator * (const Int& k, const ECPoint& P)
{
    return P * k;
}

ECPoint ECPoint::operator * (Int k) const
{
    ECPoint PP = (*this), Q { E_ };

    if (k < 0)
    {
        PP = -PP;
        k = -k;
    }
    
    //std::cout << "* PP = " << PP << std::endl;
    //std::cout << "* k  = " << k << std::endl;

    while (! E_->has_divisor())
    {
        //out ("k = %ld\n", (long) k);
        if (k.is_odd())
            Q += PP;

        if (k <= Int::One)
            break;
        k >>= 1;
        
        PP += PP;
    }

    return Q;
}

std::ostream& operator << (std::ostream& ost, const ECPoint& P)
{
    P.print();
    return ost;
}
