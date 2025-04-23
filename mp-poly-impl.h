#include <cassert>

#include "mp-int.h"


namespace
{
    int _ilog2 (int i)
    {       return i <= 0
            ? -1
            : 8 * (int) sizeof (i) - 1 - __builtin_clz ((unsigned) i);
    }

    int _next_pow2 (int i)
    {
        if (i <= 1)
            return 1;
        int p2 = 1 << _ilog2 (i);
        return p2 == i ? p2 : (p2 << 1);
    }
} // ::anon


template<class A>
Poly<A>::Poly (const As& a)
{
    a_ = a;
    set_deg();
}

template<class A>
Poly<A>::Poly (As&& a)
{
    std::swap (a_, a);
    set_deg();
}

template<class A>
Poly<A>::Poly (std::initializer_list<A> as)
{
    assert (deg_ == 0);
    for (int i = (int) as.size() - 1; i >= 0; --i)
    {
        const A& x = as.begin()[i];
        if (deg_ == 0 && (i == 0 || ! _is_zero<A>(x)))
            a_.reserve (1 + (deg_ = i));
        if (i <= deg_)
            a_.emplace (a_.begin(), x);
    }
}

template<class A>
int Poly<A>::set_deg ()
{
    deg_ = std::max<int> (0, (int) a_.size() - 1);
    while (deg_ > 0 && _is_zero<A> (a_[deg_]))
        deg_--;
    a_.resize (1 + deg_);
    return deg_;
}

template<class A>
A& Poly<A>::at (int i)
{
    if (i < 0 || i > deg_)
        fatal ("%d not in [0, %d], trace: %s", i, deg_, trace_str (4));
    return a_[i];
}

template<class A>
const A& Poly<A>::at (int i) const
{
    if (i < 0 || i > deg_)
        fatal ("%d not in [0, %d], trace: %s", i, deg_, trace_str (4));
    return a_[i];
}

template<class A>
void Poly<A>::set (int i, const A& ai)
{
    if (i > deg_ && ! _is_zero<A> (ai))
    {
        a_.resize (1 + i, Scalar::zero());
        deg_ = i;
    }
    if (i <= deg_)
        at(i) = ai;
    if (i == deg_)
        set_deg();
}

template<class A>
void Poly<A>::clear (int i)
{
    if (i <= deg_)
        at(i) = Scalar::zero();

    if (i == deg_)
        set_deg();
}


template<class A>
auto Poly<A>::operator - () const -> Poly
{
    return Poly() - (*this);
}

template<class A>
auto Poly<A>::operator - (const Poly& q) const -> Poly
{
    As r;
    int d = 0;
    for (int i = std::max (deg_, q.deg_); i >= 0; --i)
    {
        A a = i > q.deg_ ? at(i) : i > deg_ ? -q[i] : at(i) - q[i];
        if (d == 0 && (i == 0 || ! _is_zero<A> (a)))
            r.resize (d = 1 + i);
        if (d)
            r[i] = std::move (a);
    }

    return Poly (std::move (r));
}

template<class A>
auto Poly<A>::operator + (const Poly& q) const -> Poly
{
    As r;
    int d = 0;
    for (int i = std::max (deg_, q.deg_); i >= 0; --i)
    {
        A a = i > q.deg_ ? at(i) : i > deg_ ? q[i] : at(i) + q[i];
        if (d == 0 && (i == 0 || ! _is_zero<A> (a)))
            r.resize (d = 1 + i);
        if (d)
            r[i] = std::move (a);
    }

    return Poly (std::move (r));
}

template<class A>
auto Poly<A>::operator * (const Poly& q) const -> Poly
{
    As r;
    int d = 0;
    for (int i = deg_ + q.deg_; i >= 0; --i)
    {
        A a = Scalar::zero();
        for (int k = 0; k <= deg_; ++k)
            if (i - k >= 0 && i - k <= q.deg_)
                a += at(k) * q[i-k];
        if (d == 0 && (i == 0 || ! _is_zero<A> (a)))
            r.resize (d = 1 + i);
        if (d)
            r[i] = std::move (a);
    }

    return Poly (std::move (r));
}


template<class A>
auto Poly<A>::operator *= (const A& a) -> Poly&
{
    for (int i = 0; i <= deg_; i++)
        at(i) *= a;
    set_deg();
    return *this;
}

template<class A>
auto Poly<A>::operator /= (const A& a) -> Poly&
{
    if (_is_zero<A> (a))
        fatal ("Division by zero");

    for (int i = 0; i <= deg_; i++)
        at(i) /= a;
    set_deg();
    return *this;
}

template<class A>
bool Poly<A>::operator == (const Poly& q) const
{
    if (deg_ != q.deg_)
        return false;

    for (int i = 0; i <= deg_; ++i)
        if (at(i) == q[i])
            continue;
        else if (_isnan<A>(at(i)) && _isnan<A>(q[i]))
            continue;
        else
            return false;

    return true;
}

template<class A>
bool Poly<A>::polynomial_p() const
{
    for (int i = 0; i <= deg_; ++i)
        if (! _isnumber<A> (at(i)))
            return false;

    return true;
}

template<class A>
bool Poly<A>::zero_p() const
{
    return deg_ == 0 && _is_zero<A> (at(0));
}

// To which order  x  divides the polynomial.  ord_x of 0-polynomial is = 0.
template<class A>
int Poly<A>::ord_x() const
{
    for (int i = 0; i <= deg_; ++i)
        if (!_is_zero<A> (at(i)))
            return std::max (0, i - 1);

    assert (deg_ == 0);
    return 0;
}


template<class A>
template <class X>
auto Poly<A>::operator () (const X& x) const -> result_type<X>
{
    using Y = result_type<X>;

    Y y ( _make<Y,A>(at(deg_)) );

    for (int i = deg_ - 1; i >= 0; --i)
        y = _fma (x, y, _make<Y,A>(at(i)) );
    return y;
}

template<class A>
A Poly<A>::operator () (const A& x) const
{
    A y { at(deg_) };

    for (int i = deg_ - 1; i >= 0; --i)
        y = _fma (y, x, at(i));

    return y;
}

template<class A>
Poly<A> Poly<A>::operator () (const Poly& p) const
{
    Poly q { at(deg_) };

    for (int i = deg_ - 1; i >= 0; --i)
        q = q * p + at(i);

    return q;
}

template<class A>
Poly<A> Poly<A>::operator << (int s) const
{
    if (s < 0 || s > 10000)
        fatal ("bad shift offset %d", s);
    Poly q;
    q.a_.resize (s + 1 + deg_, Scalar::zero());
    for (int i = s + deg_; i >= s; --i)
        q.a_[i] = at (i - s);
    q.set_deg();
    return q;
}

template<class A>
Poly<A> Poly<A>::operator >> (int s) const
{
    if (s < 0 || s > 10000)
        fatal ("bad shift offset %d", s);
    Poly q;
    int d = std::max (0, deg_ - s);
    q.a_.resize (1 + d);
    if (s <= deg_)
    {
        for (int i = deg_; i >= s; --i)
            q.a_[i - s] = at (i);
    }
    q.set_deg();
    return q;
}

template<class A>
Poly<A>& Poly<A>::operator <<= (int s)
{
    if (s < 0 || s > 10000)
        fatal ("bad shift offset %d", s);
    a_.resize (s + 1 + deg_, Scalar::zero());
    for (int i = s + deg_; i >= 0; --i)
        if (i >= s)
            a_[i] = at (i - s);
        else
            a_[i] = Scalar::zero();
    deg_ += s;
    set_deg();
    return *this;
}

template<class A>
Poly<A> Poly<A>::D () const
{
    Poly q;
    if (deg_ == 0)
        return q;
    q.a_.resize (deg_);
    for (int i = 0; i < deg_; ++i)
        q.a_[i] = at(i+1) * Scalar::get(i + 1);
    q.set_deg();
    return q;
}

template<class A>
Poly<A> Poly<A>::crop (int lo, int hi, int stride) const
{
    int d = (hi - lo) / stride;

    Poly f;
    f.a_.resize (1 + d);

    for (int i = 0; i <= d; i++, lo += stride)
    {
        f.a_[i] = at (lo);
    }
    f.set_deg();

    return f;
}

template<class A>
Poly<A> Poly<A>::pow (const Int& ex) const
{
    if (ex < 0)
        error ("todo Poly<A>::pow (x < 0), x = %ld", (long) ex);

    Poly q { Scalar::one() };
    if (!ex.is_zero())
        for (int b = ex.msbit(); ; --b)
        {
            if (ex.bit (b))
                q *= *this;
            if (b == 0)
                break;
            q *= q;
        }
    return q;
}

template<class A>
Poly<A> Poly<A>::swap () const
{
    Poly p;
    p.a_.resize (1 + deg_);

    for (int i = 0; i <= deg_; ++i)
        p.a_[i] = at (deg_ - i);

    p.set_deg();

    return p;
}

template<class A>
auto Poly<A>::height() const -> scalar_real_type
{
    using R = scalar_real_type;
    R h { _make<R>(0) };
    for (const auto& a : a_)
        h = _max<R> (h, _abs<R,A>(a));
    return h;
}


// Return Self - p * c * x^r with the assertion
// that the highest order term will vanish.
template<class A>
Poly<A> Poly<A>::msub (const Poly& p, const A& c, int r) const
{
    if (deg() - r < p.deg())
        fatal ("deg=%d, r=%d, p.deg=%d, need deg - r >= p.deg", deg(), r, p.deg());

    // 0 ... r - 1
    Poly q = crop (0, r - 1);

    // deg() - 1 ... r
    for (int k = deg() - 1; k >= r; --k)
        q.set (k, at(k) - c * p[k - r]);

    return q;
}

template<class A>
Poly<A> Poly<A>::divmod (const Poly& p, Poly& rem) const
{
    if (p.zero_p())
        fatal ("divmod by 0, trace: %s", trace_str (4));

    A q_k;
    Poly q{}, r = *this;
    int deg_q = deg() - p.deg();

    if (deg_q >= 0)
    {
        for (int k = deg_q; k >= 0; --k)
        {
            // Eliminate r[k + p.deg] with p's leading term.

            if (r.deg() >= k + p.deg())
            {
                q_k = r[k + p.deg()] / p[p.deg()];
                q.set (k, q_k);

                if (! _is_zero<A>(q_k))
                    r = r.msub (p, q_k, k);
            }
            assert (r.deg() < k + p.deg());
        }
    }

    rem = r;
    return q;
}

template<class A>
Poly<A> Poly<A>::operator / (const Poly& p) const
{
    if (p.zero_p())
        fatal ("div by 0, trace: %s", trace_str (4));
    Poly rem;
    return divmod (p, rem);
}

template<class A>
Poly<A> Poly<A>::operator % (const Poly& p) const
{
    if (p.zero_p())
        fatal ("mod by 0, trace: %s", trace_str (4));
    Poly rem;
    (void) divmod (p, rem);
    return rem;
}


template<class A>
std::ostream& operator << (std::ostream& ost, const Poly<A>& q)
{
    return q.print (ost);
}

/////////////////////////////////////////////////////////////////////
// RationalFunction = RationalFunction <op> RationalFunction

template<class A>
auto RationalFunction<A>::operator - () const -> Self
{
    return Self { -p_, q_ };
}

template<class A>
auto RationalFunction<A>::operator - (const Self& b) const -> Self
{
    return Self { p_ * b.q_ - q_ * b.p_, q_ * b.q_ };
}

template<class A>
auto RationalFunction<A>::operator + (const Self& b) const -> Self
{
    return Self { p_ * b.q_ + q_ * b.p_, q_ * b.q_ };
}

template<class A>
auto RationalFunction<A>::operator * (const Self& b) const -> Self
{
    return Self { p_ * b.p_, q_ * b.q_ };
}

template<class A>
auto RationalFunction<A>::operator / (const Self& b) const -> Self
{
    return Self { p_ * b.q_, q_ * b.p_ };
}


/////////////////////////////////////////////////////////////////////
// RationalFunction = RationalFunction <op> Poly

template<class A>
auto RationalFunction<A>::operator - (const Poly<A>& b) const -> Self
{
    return Self { p_ - q_ * b, q_ };
}

template<class A>
auto RationalFunction<A>::operator + (const Poly<A>& b) const -> Self
{
    return Self { p_ + q_ * b, q_ };
}

template<class A>
auto RationalFunction<A>::operator * (const Poly<A>& b) const -> Self
{
    return Self { p_ * b, q_ };
}

template<class A>
auto RationalFunction<A>::operator / (const Poly<A>& b) const -> Self
{
    return Self { p_, q_ * b };
}

/////////////////////////////////////////////////////////////////////
// Poly<A>::Scalar

template<class A>
A Poly<A>::Scalar::get (int i)
{
    return _make<A> (i);
}


/////////////////////////////////////////////////////////////////////
// RationalFunction = Poly <op> RationalFunction

template<class A>
RationalFunction<A> Poly<A>::operator - (const RationalFunction<A>& b) const
{
    return RationalFunction<A> { (*this) * b.q_ - b.p_, b.q_ };
}

template<class A>
RationalFunction<A> Poly<A>::operator + (const RationalFunction<A>& b) const
{
    return RationalFunction<A> { (*this) * b.q_ + b.p_, b.q_ };
}

template<class A>
RationalFunction<A> Poly<A>::operator * (const RationalFunction<A>& b) const
{
    return RationalFunction<A> { (*this) * b.p_, b.q_ };
}

template<class A>
RationalFunction<A> Poly<A>::operator / (const RationalFunction<A>& b) const
{
    return RationalFunction<A> { (*this) * b.q_, b.p_ };
}


/////////////////////////////////////////////////////////////////////
// Eval at

template<class A>
auto RationalFunction<A>::at_1_over_x() const -> Self
{
    auto p = p_.swap();
    auto q = q_.swap();
    auto m = std::max (p_.deg(), q_.deg());
    return Self { p << (m - p_.deg()), q << (m - q_.deg()) };
}

template<class A>
auto RationalFunction<A>::operator () (const Poly<A>& b) const -> Self
{
    return Self { p_(b), q_(b) };
}

template<class A>
auto RationalFunction<A>::operator () (const Self& b) const -> Self
{
    return p_(b) / q_(b);
}

template<class A>
bool RationalFunction<A>::is_polynomial () const
{
    return q_.deg() == 0 && q_[0] == Poly<A>::Scalar::one();
}

template<class A>
auto RationalFunction<A>::normalize_p (int i) const -> Self
{
    assert (i >= 0 && i <= p_.deg());
    Poly<A> p = p_ / p_[i];
    p.set (i, Poly<A>::Scalar::one());
    return Self { p, q_ / p_[i] };
}

template<class A>
auto RationalFunction<A>::normalize_q (int i) const -> Self
{
    assert (i >= 0 && i <= q_.deg());
    Poly<A> q = q_ / q_[i];
    q.set (i, Poly<A>::Scalar::one());
    return Self { p_ / q_[i], q };
}

template<class A>
void RationalFunction<A>::print_VHDL_Tables (
    std::ostream& ost, const char* varP, const char *varQ, const char *typ,
    int prec2) const
{
    int max_deg = std::max (p_.deg(), q_.deg());
    int n_coeff = 1 + max_deg;
    int tabsize = _next_pow2 (n_coeff) - 1;

    ost << "TYPE " << typ << " IS ARRAY(0 TO " << tabsize << ")"
        << " OF std_logic_vector(31 DOWNTO 0);"
        << " -- Fuer Tabellen " << varP << " und " << varQ << "." << std::endl;
    ost << "CONSTANT MaxDeg : INTEGER := " << max_deg << ";"
        << " -- MAX (Grad P, Grad Q) = " << max_deg << std::endl;
    p_.print_VHDL_Table (ost, varP, typ, prec2, n_coeff);
    q_.print_VHDL_Table (ost, varQ, typ, prec2, n_coeff);
}
