#include "mp-factor.h"
#include "mp-mint.h"

#define DIAGNOSTIC_NO_FORMAT_CHECK
#include "diagnostic.h"

#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <set>

template <>
bool Power<Int>::operator() (const Power& a, const Power& b) const
{
    return a.base_ < b.base_;
}

template <>
bool Power<GaussInt>::operator() (const Power& pa, const Power& pb) const
{
    const GaussInt& a = pa.base_;
    const GaussInt& b = pb.base_;
    Int na = a.norm();
    Int nb = b.norm();

    if (na != nb)
        return na < nb;

    if (a.x() != b.x())
        return a.x() < b.x();

    return a.y() < b.y();
}

template <class X>
int Factors<X>::is_fully_factored () const
{
    if (stop_reason_ || ! composites_.empty())
        return 0;
    return 1 + pseudo_primes_.empty();
}


template <class X>
void Factors<X>::print_powers (bool with_unit, const char* end) const
{
    //cout << "size: " << ps.size() << endl;
    bool first = true;
    if (with_unit && ! unit_.is_one())
    {
        std::cout << unit_;
        first = false;
    }

    for (const auto& p : powers_)
    {
        if (!first)
            std::cout << " * ";
        first = false;
        std::cout << p;
    }
    std::cout << end;
}


template <class X>
int Factors<X>::add (const X& d, Kind kind)
{
    assert (! d.is_zero());
    if (d.is_unit())
        return stop_reason_;

    switch (kind)
    {
        default:
            assert (0);
        case Unknown:
        {
            switch (d.is_probab_prime (30))
            {
                default: assert (0);
                case 2: return add (d, Prime);
                case 1: return add (d, PseudoPrime);
                case 0: return add (d, Composite);
            }
        }

        case Prime:
            primes_.push_back (d);
            break;

        case PseudoPrime:
            pseudo_primes_.push_back (d);
            break;

        case Composite:
            composites_.push_back (d);
            break;
    }

    if (!stop_reason_ && callback_)
        stop_reason_ = callback_ (d, kind, *this);

    return stop_reason_;
}


static void factor_no_small_primes (const Int&, Factors<Int>&);

template<>
std::ostream& operator << (std::ostream& ost, const Power<Int>& ip)
{
    bool paren = ip.base_ >= 0;
    const char *l = "(" + paren;
    const char *r = ")" + paren;
    
    ost << l << ip.base_ << r;
 
    if (ip.base_ == 0 || ip.expo_ != 1)
        ost << "^" << ip.expo_;
    return ost;
}

template<>
std::ostream& operator << (std::ostream& ost, const Power<GaussInt>& ip)
{
    ost << ip.base_;

    if (ip.base_ == 0 || ip.expo_ != 1)
        ost << "^" << ip.expo_;
    return ost;
}

template<class X>
X Factors<X>::print (const char *text, const Xs& qs) const
{
    X f { X::One };
    bool first = true;
    printf ("%s[%d]%s: ", text, (int) qs.size(), "               " + strlen (text));
    if (qs.empty())
        printf ("(none)");
    else
        for (const auto& q : qs)
        {
            if (!first)
                printf (" * ");
            q.print();
            f *= q;
            first = false;
        }
    out ("\n");
    return f;
}


template<class X>
void Factors<X>::print() const
{
    X f1 = print ("Primes", primes_);
    X f2 = print ("Pseudo-Primes", pseudo_primes_);
    X f3 = print ("Composite", composites_);
    std::cout << "Product: " << (f1 * f2 * f3) << std::endl;
    if (stop_reason_)
        printf ("Stop Reason: %d\n", stop_reason_);
    out ("");
}

template <class X>
std::ostream& operator << (std::ostream& ost, const Factors<X>& f)
{
    f.print();
    return ost;
}

int callback_stop_1_mod_4 (const Int& d, Factors<Int>::Kind, Factors<Int>&)
{
    switch ((long) d.mod(4))
    {
        case 0: return 2;
        case 1: return 0;
        case 2: return 2;
        case 3: return 3;
    }
    assert (0);
    return 0;
}


template <class X>
void Factors<X>::set_powers ()
{
    assert (powers_.empty());

    unit_ = n_;
    for (const auto& ds : { primes_, pseudo_primes_, composites_ })
        for (const X& d : ds)
        {
            auto f = powers_.find (d);
            if (f == powers_.end())
                powers_.insert (Power<X> (d));
            else
                f->addto_expo (1);
            unit_ = unit_ / d;
        }
}


template std::ostream& operator << (std::ostream&, const Factors<Int>&);
template std::ostream& operator << (std::ostream&, const Factors<GaussInt>&);

template class Factors<Int>;
template class Factors<GaussInt>;

Int Int::phi() const
{
    return ::phi (*this);
}

Int phi (const Int& n)
{
    return phi (factor (n));
}

Int phi (const Factors<Int>& ns)
{
    if (!ns.is_fully_factored())
        error ("need all prime factors of %Zd to compute phi(n)", ns.n_.mpz());

    Int phi_n = 1;
    for (const auto& be : ns.powers_)
        phi_n *= be.value() - be.base_.pow(be.expo_ - 1);

    return phi_n;
}


static void remove_known (Int& n, Factors<Int>& fs, const Factors<Int> *know)
{
    int ord;
    for (const auto& p : know->primes_)
    {
        n = n.remove_factor (p, &ord);
        while (ord--)
            if (fs.add_prime (p))
                return;
    }
    for (const auto& p : know->pseudo_primes_)
    {
        n = n.remove_factor (p, &ord);
        while (ord--)
            if (fs.add_pseudo_prime (p))
                return;
    }
}

Factors<Int> factor (const Int& n, const Factors<Int>::Context *ctx)
{
    assert (n >= 1);
    Factors<Int> f;
    f.n_ = n;
    Int r = n;
    if (ctx)
    {
        f.callback_ = ctx->callback_;
        f.verbose_ = ctx->verbose_;
        if (ctx->known_)
            remove_known (r, f, ctx->known_);
    }

    if (r != 1L && !f.stop_reason_)
        r = remove_small_primes (r, f);

    if (r != 1L && !f.stop_reason_)
        factor_no_small_primes (r, f);

    f.set_powers();

    return f;
}

Factors<GaussInt>
factor (const GaussInt& g, const Factors<GaussInt>::Context *ctx)
{
    return g.factor (ctx ? ctx->verbose_ : 0);
}

Factors<Int> factor_a_pow_n_minus_1 (const Int& a, int n,
                                     const Factors<Int>::Context* ctx0)
{
    Int a_n_1 = a.pow(n) - 1;
    Factors<Int> ns = factor (n);//if (n > 1 && !is_probab_prime (n, 25))

    Factors<Int> known;
    Factors<Int>::Context ctx;
    if (ctx0)
        ctx = *ctx0;
        
    //Factors<Int>::Context* ctx;
    if (!ctx.known_
        && ns.powers_.size() > 0)
    {
        int p = (long) ns.powers_.begin()->base_;
        if (p != n)
        {
            Int a_np_1 = a.pow (n / p) - 1;
            assert (a_n_1 % a_np_1 == 0);
            //out ("n = %d, p = %d, n_p = %d\n", n, p, n / p);
            known = factor_a_pow_n_minus_1 (a, n / p);
            ctx.known_ = &known;
            out ("known = ");
            known.print_powers();
        }
    }

    return factor (a_n_1, ctx.known_ ? &ctx : nullptr);
}


#define ARRAY_SIZE(X) (sizeof (X) / sizeof (*X))

static Int find_special_factor (const Int& n, int verbose)
{
    static const char* const cheat[] =
    {
        // For around 10^50.
        "100000000000000000000000000000000000000000000000027", "2587066943291159687641",
        "1028226990156474555135022141325892534444833", "5271211709451252581209",
        "127481664041362394758866743333714482261873", "39155968826707001",
        "1922500949186781137243516985555423693351563", "56041143257705677",
        "1966413654776418767451921186140716561135800527", "137513949200517854220769",
        "20000000000000000000000000000000000000000000000009", "189198096947243195029",
        "2045910225459306845615614386840705429845738369", "12260450030546096531623",
        "236342924129194496045982879318576081150706429", "511887232216464623",
        "250436615450294139302958204090270865721", "272555294827355857",
        "2555192150449713818479149632052330335241210139", "13407792350136185977247",
        "28816517166575606472420287750213818557375991", "36357872341129309374913",
        "3052596233096248359229524710766506914130467", "4151251613742235519351",
        "319489172085439979705025423782179081201729", "25883975674529228197",
        "33333333333333333333333333333333333333333333333323", "89220954141915940523",
        "44081494843861927921757710319531203753", "81339056487339199",
        "50000000000000000000000000000000000000000000000021", "34963514888753421623",
        "6952539324605300866272494767171277335820303", "198564406462335103",
        "77942322681215900233826968043647700701480904131", "2643874458618067",
        "78554595443833464257659073055773762765121759623", "18733296487763313457", 
        "96805421103581800580832526621490803484995159729", "317115696148294961",
        "99999999999999999999999999999999999999999999999941", "198901294156453439",
    };

    static long min_ilog2 = LONG_MAX;

    if (min_ilog2 == LONG_MAX)
        for (int i = 0; i < (int) ARRAY_SIZE (cheat); i += 2)
            min_ilog2 = std::min (min_ilog2, Int::from_str(cheat[i]).ilog2());

    if (n.ilog2() < min_ilog2)
        return Int::One;

    for (int i = 0; i < (int) ARRAY_SIZE (cheat); i += 2)
    {
        Int d = Int::from_str (cheat[i]);
        if (n == d)
        {
            d.set_str (cheat[i+1]);
            if (verbose & 1)
                out ("Cheat(%Zd -> %Zd\n", n.mpz(), d.mpz());
            assert (n % d == 0);
            return d;
        }
    }

    return Int::One;
}

Int remove_small_primes (const Int& n, Factors<Int>& f)
{
    if (n.is_one())
        return n;

    Int r { n };
    for (Int d = 2; d <= 1000; d = d.next_prime())
    {
        while ((r % d).is_zero())
        {
            f.add_prime (d);
            r /= d;
            if (f.stop_reason_ || r.is_one())
            {
                f.add (r);
                return Int::One;
            }
        }

        if (d * d > n)
        {
            assert (r > 1 && r <= n && n % r == 0);
            f.add_prime (r);
            return Int::One;
        }
    }
    return r;
}


void factor_no_small_primes (const Int& n, Factors<Int>& f)
{
    if (f.stop_reason_)
        return;

    switch (n.is_probab_prime (25))
    {
        case 0:
            break;
        case 1:
            f.add_pseudo_prime (n);
            return;
        case 2:
            f.add_prime (n);
            return;
        default: assert (0);
    }

    // Composite.

    assert (f.n_ != 0);
    if (f.callback_
        && !f.stop_reason_
        && n > 1
        && ((f.stop_reason_ = f.callback_ (n, Factors<Int>::Composite, f))))
    {
        return;
    }

    Int q;

    q = find_factor_Pollard_p_minus_1 (n, 20000 /* smooth */, f.verbose_);
    if (f.stop_reason_ || (q > Int::One && q < n))
    {
done:;
        if (!f.stop_reason_)
            factor_no_small_primes (q, f);
        if (!f.stop_reason_)
            factor_no_small_primes (n / q, f);
        return;
    }

    if (n.is_perfect_square())
    {
        q = n.sqrt();
        goto done;
    }

    q = find_special_factor (n, f.verbose_);
    if (q > Int::One && q < n)
        goto done;
/*
    Int n_2root4 = 4 * n.sqrt().sqrt();
    int n_rho = (long) n_2root4.min (100000);

    for (int i = 0; i < 10; i++)
    {
        q = find_factor_Pollard_rho (n, n_rho, 1+f.verbose);
        if (f.stop_reason_ || (q > Int::One && q < n))
            goto done;
    }
*/
/*
    int smooth = 100000;

    for (int i = 0; i < 0; i++)
    {
        q = find_factor_elliptic_curve (n, smooth);
        if (q != 1)
            goto found;
    }
*/

    bool up = true;
    int b1 = 10000;
    for (int i = 0; i < 50; ++i)
    {
        int result;
        q = find_factor_ECM (n, b1, f.verbose_, &result);

        if (f.stop_reason_ || (q > Int::One && q < n))
            goto done;
        if (up)
        {
            up = !(result && q == n);
            if (up)
                b1 = std::min (2500000.0, b1 * 1.2);
        }
        if (!up)
        {
            b1 *= 0.9;
        }
    }

    f.add_composite (n);
}


Int find_factor_Pollard_p_minus_1 (const Int& n, int smooth, int verbose)
{
    MInt::Modulus m = MInt::make_modulus (n);

    if (verbose & 1) out ("Pollard.p-1(%Zd -> %d", n.mpz(), smooth);

    MInt a0{m, 2}, a_ex;
    Int k, ex, phi_n = 1, d = 1;

    do // Cold loop if we find a factor of phi(n).
    {
        MInt a = a0;
        for (k = 2; k < smooth; k = k.next_prime())
        {
            int i = std::log (smooth) / std::log ((long) k);
            ex = k.pow (i);
            a_ex = a.pow (ex);
            if (a_ex == 1)
            {
                if (verbose & 1) out (" -> %ld", (long) k);
                break;
            }
            a = a_ex;
        }

        d = a == 1 ? Int::One : n.gcd (a.z() - 1);

        if (d > 1 || a_ex != 1)
            break;

        phi_n *= k;
        a0 = a0.pow (phi_n);
        assert (phi_n > 1 && k > 1);
    } while (phi_n < n);
    
    if (verbose & 1)
    {
        out (")");
        if (d > 1)
            out (" [%Zd]", d.mpz());
        out ("\n");
    }
 
    return d;
}


Int find_factor_Pollard_rho (const Int& n, int n_loops, int verbose)
{
    Int x = 1;
    Int y = 1;
    Int c = 1 + Int::random (n - 1);

    if (verbose & 1) out ("Pollard.Rho(%Zd) # %d: %Zd\n", n.mpz(), n_loops, c.mpz());

    while (n_loops--)
    {
        x = (x * x + c).mod (n);
        y = (y * y + c).mod (n);
        y = (y * y + c).mod (n);
        Int g = n.gcd (x - y);
        if (g != 1 && g != n)
            return g;
    }

    return 1;
}

#include "ecm.h"

Int find_factor_ECM (const Int& n0, int b1, int verbose, int *result)
{
    Int f, n = n0;
    ecm_params ep;
    ecm_init (ep);

    if (verbose & 1) out ("ECM-GMP(%Zd B1=%d)", n.mpz(), b1);

    ep->verbose = (verbose & 2) != 0;
    *result = ecm_factor (f.mpz(), n.mpz(), b1, ep);
    if (verbose & 1)
    {
        if (*result || (f > 1 && f < n))
            out (" %d: %Zd", *result, f.mpz());
        out ("\n");
    }

    ecm_clear (ep);
    return f;
}

#if 0

Int find_factor_elliptic_curve (const Int& n, int smooth)
{
    MInt::Modulus m = MInt::make_modulus (n);

    ECPoint P = ECPoint::random (m);
    std::cout << *P.E();
    fflush(stdout);
    //printf ("Order = %ld\n", (long) order (P));
    for (Int k = 2; k < smooth && !P.is_zero(); k = k.next_prime())
    {
        int i = std::log (smooth) / std::log ((long) k);
        P *= k.pow (i);
        //k *= i;
        //cout << i << ": " << P << " = * " << k << endl;
        if (P.E()->has_divisor())
            return P.E()->divisor();
    }

    return 1;
}
#endif
