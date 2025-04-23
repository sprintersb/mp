#include "diagnostic.h"

#include "mp-approx.h"
#include "mp-matrix.h"

/*
    Approx::Target

        .n      Degree of the numerator polynomial.

        .m      Degree of the denominator polynomial.  If .m = 0, then the
                approximation will be a polynomial of degree .n, otherwise
                it will be a rational function of degree [.n/.m].

        .f      The function to be approximated as Approx::Func, i.e.
                Float (*)(const Float&).

        .a, .b  Interval over which .f is to be approxed.
*/

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm> // std::sort


static const double _pi = 3.1415926535897932384626433832795;

using std::cout;
using std::endl;

Approx::Verbose Approx::verbose;

#define VERB Approx::verbose

std::ostream& operator << (std::ostream& ost, const Approx::DValue& v)
{
    return ost << "(" << v.x << ", " << v.y << ")";
}

namespace
{
    template <class X> std::string _tostr (const X& x)
    {
        std::stringstream buffer;
        buffer << x;
        return buffer.str();
    }
} // anon


struct Approx::Quality
{
    struct SignGroup;

    const Approx *approx;
    // Holds x, delta f(x) Value's with delta != 0.
    Tree tree;
    // Holds x, delta f(x) Value's with delta = 0.
    Tree tree_delta0;
    int n_groups = 0;
    // Allow that more SignGroups than needed, and ditch them after
    // make_groups().
    int allow_more_groups = 3;
    //mutable
    std::vector<SignGroup> groups;

    Quality (const Approx *a) : approx(a) {}
    Approx::Failure init (int flags);
    void try_flip (Tree::Node*);
    const Tree::Node* add (const Value&);
    const Tree::Node* add (const F&);
    void make_groups();
    void ditch_excess_groups();
    void find_worse (int flags);
    void find_zeros();
    FVector<F> adjust_worst_place0() const;
    static bool sign_flip (const F& a, const F& b)
    {
        return a.sgn() * b.sgn() < 0;
    }
    static bool sign_flip (const Tree::Node*, const Tree::Node*);
    double xacos (const F&) const;
    FMatrix<F> place0_matrix (const Target&) const;
    FVector<F> place0_vector (const Target&) const;
};

struct Approx::Quality::SignGroup
{
    int id;
    int sign = 0;
    Value worst;
    SignGroup *prev = nullptr;
    SignGroup *next = nullptr;
    const Approx *approx = nullptr;
    Quality *quality; // Owned by Approx::Result
    const Tree::Node *worst_node = nullptr;
    const Tree::Node *lo = nullptr;
    const Tree::Node *hi = nullptr;
    // Node with y ~= 0 at the left of this group.  This sign might not be the
    // same like .sign.
    const Tree::Node *node_zero = nullptr;
    Value zero;
    bool find_zero();
    Failure find_worse();
    bool find_worse (int how, const Values&);
    void set_worst (const Tree::Node*);
    bool maybe_set_worst (const Tree::Node*);
    Failure worst_fix_neighbours();
    F width() const;
    // x's that are used for width computation and whose delta(x) might have
    // other sign than .sign.
    const F& width_x_left() const;
    const F& width_x_right() const;
    bool is_worse (const F& y) const
    {
        return (sign < 0 && y < worst.y) || (sign > 0 && y > worst.y);
    };
};

struct Approx::Remez
{
    static FVector<F> vector (const Values&, const Target&);
    static FMatrix<F> matrix (const Values&, const Target&, const F&);
    static Poly<F> poly (int, const F*);
    static Poly<F> poly_p (const F*, const Target::PolynomialDescriptor&);
    static Poly<F> poly_q (const F*, const Target::PolynomialDescriptor&);
};

struct Approx::Place0
{
    static FVector<F> vector (const Values&);
    static FMatrix<F> matrix (const Quality&, const Target&);
    static Poly<F> poly (int, const F*);
};


// Return a polynomial that hits all v = (x,y) in v's.
auto Approx::Value::polynomial (const Values& vs) -> Poly<F>
{
    int n = (int) vs.size();
    FMatrix<F> m { Alloc (n, n) };
    FVector<F> b { Alloc (n) };

    for (int i = 0; i < n; ++i)
    {
        b[i] = vs[i].y;
        m[i][0] = 1;
        for (int j = 1; j < n; ++j)
            m[i][j] = vs[i].x * m[i][j-1];
    }

    FMatrix<F>::push_verbose (std::max (0, Approx::verbose.level - 3));
    FVector<F> x = m.solve (b, Algo::QR_vectors);
    FMatrix<F>::pop_verbose();
    Poly<F> p;
    for (int i = 0; i < n; ++i)
        p.set (i, x[i]);
    return p;
}

// Right-hand side for linear system with matrix from below.
// For explanation see there.
auto Approx::Remez::vector (const Values& vs, const Target& t) -> FVector<F>
{
    int dim = t.symmetry
        ? t.p.n_non0 + t.q.n_non0
        : vs.size();

    FVector<F> v { Alloc (dim) };
    for (int i = 0; i < dim; ++i)
        v[i] = vs[i].y;
    return v;
}

auto Approx::Remez::matrix (const Values& vs, const Target& t, const F& E_old)
    -> FMatrix<F>
{
    const int dim = t.p.n_non0 + t.q.n_non0;

    if (! t.symmetry
        && dim != (int) vs.size())
    {
        for (int i = 0; i < (int) vs.size(); ++i)
            printf ("v[%d] = (%.3f, %.ee)\n",
                    i, (double) vs[i].x, (double) vs[i].y);
        fatal ("dim=%d != vs.size=%d", dim, (int) vs.size());
    }

    if (t.symmetry
        && dim > (int) vs.size())
    {
        fatal ("dim=%d > vs.size=%d", dim, (int) vs.size());
    }

    if (t.symmetry)
        info (VERB.symm, "using v[%d] = %.5e ... v[%d] = %.5e",
              0, (double) vs[0].x, dim - 1, (double) vs[dim - 1].x);

#if 1
    FMatrix<F> m { Alloc (dim, dim) };
    F one(1);

    // Holds pre-computed powers of x.
    std::vector<F> xk; // x^k
    xk.resize (1 + std::max (t.p.deg, t.q.deg));
    xk[0] = one;

    // Number of columns is p.n_non0 + q.n_non0:
    //
    // Values  | #            | Column number(s)
    // --------+--------------+-------------------
    // E       | 1            | 0
    // p0...pn | p.n_non0     | 1 ... p.n_non0
    // q1...qm | q.n_non0 - 1 | p.n_non0 + 1 ... p.n_non0 + q.n_non0 - 1

    // Set i-th row of M.
    for (int i = 0; i < dim; ++i)
    {
        const F& x = vs[i].x;

        F Gi = t.error_absolute() ? one : vs[i].y;
        if (i & 1) Gi = -Gi;
        m[i][0] = -Gi;

        // Pre-compute k-th powers of x in xk[].
        for (int k = 1; k <= std::max (t.p.deg, t.q.deg); ++k)
            xk[k] = k == 1 ? x : x * xk[k - 1];

        /* We want to solve this equation for x = xi's and E simultaneously:

                          p(x)
                        --------   -  f(x)  =  Gi
                        1 + q(x)

           Gi = (-1)^i * E          ;  absolute error
           Gi = (-1)^i * E * f(xi)  ;  relative error

           Multiplying through by the denominator:

                   p(x) + q(x) * (-Gi - f(x)) - Gi  =  f(x)

                                 qqqqqqqqqqqq   EE

           where we have to guess the "-Gi" part in the factor for q.
           One guess is based on E_old, the E from the previous iteration.  */

        int kp = t.p.symmetry == SymmetryOdd;
        int sp = 1 + !! t.p.symmetry;
        for (int j = 1; j <= t.p.n_non0; ++j, kp += sp)
            m[i][j] = xk[kp];

        int kq = 0;
        int sq = 1 + !! t.q.symmetry;
        // Skip the lowest power of q() and set implicit to 1.
        kq += sq;
        for (int j = 1; j <= t.q.n_non0 - 1; ++j, kq += sq)
            m[i][j + t.p.n_non0] = xk[kq] * (E_old * (-Gi) - vs[i].y);
    } // loop i
#else
    int nm = std::max (t.n, t.m);
    FMatrix<F> m { Alloc (dim, dim) };
    F one(1), xj;

    // Number of columns is n + m + 2:
    // Values  | #     | Column number(s)
    // --------+-------+--------------------
    // E       | 1     | 0
    // p0...pn | n + 1 | 1 ... n + 1
    // q1...qm | m     | n + 2 ... n + m + 1
    for (int i = 0; i < dim; ++i)
    {
        const F& fi = t.error_absolute() ? one : vs[i].y;
        m[i][0] = (i & 1) ? fi : -fi;
        m[i][1] = one; // x^0
        const F& x = vs[i].x;
        xj = x; // x^j
        for (int j = 1; j <= nm; ++j, xj *= x)
        {
            if (j <= t.n)
                m[i][j + 1] = xj;
            if (j <= t.m)
                m[i][j + 1 + t.n] = xj * (E_old * m[i][0] - vs[i].y);
        }
    }
#endif // 1
    return m;
}

// Turn v[0] ... v[deg] into a polynomial.
auto Approx::Remez::poly (int deg, const F *v) -> Poly<F>
{
    Poly<F> p;
    for (int d = 0; d <= deg; ++d)
        p.set (d, *v++);
    if (! p.polynomial_p())
        error ("This is not a polynomial: %s", _tostr(p).c_str());

    return p;
}

auto Approx::Remez::poly_p (const F *v, const Target::PolynomialDescriptor &d)
    -> Poly<F>
{
    Poly<F> p;

    int kp = d.symmetry == SymmetryOdd;
    int sp = 1 + !! d.symmetry;
    for (int j = 0; j < d.n_non0; ++j, kp += sp)
        p.set (kp, *v++);

    if (! p.polynomial_p())
        error ("This is not a polynomial: %s", _tostr(p).c_str());

    if (p.deg() != d.deg)
        warning ("p.deg = %d has not expected degree of %d", p.deg(), d.deg);

    return p;
}

auto Approx::Remez::poly_q (const F *v, const Target::PolynomialDescriptor &d)
    -> Poly<F>
{
    Poly<F> q;

    assert (d.symmetry != SymmetryOdd);
    int sq = 1 + !! d.symmetry;
    q.set (0, F{1});
    // Skip the lowest power of q() because we set it to 1 above.
    int kq = sq;
    for (int j = 1; j < d.n_non0; ++j, kq += sq)
        q.set (kq, *v++);

    if (! q.polynomial_p())
        error ("This is not a polynomial: %s", _tostr(q).c_str());

    if (q.deg() != d.deg)
        warning ("q.deg = %d has not expected degree of %d", q.deg(), d.deg);

    return q;
}


void Approx::Target::set_poly_descriptors() const
{
    p.symmetry = q.symmetry = SymmetryNone;

    if (symmetry)
    {
        if (m & 1)
            error ("denominator cannot have symmetry and odd degree in [%d/%d]",
                   n, m);
        if (m
            && (n + m + symmetry) % 2 != 0)
        {
            const char *sym = symmetry == SymmetryOdd ? "odd" : "even";
            error ("requesting %s symmetry, but degree [%d/%d] is not %s",
                   sym, n, m, sym);
        }
        if (n == 0 && m == 0
            && symmetry == SymmetryOdd)
        {
            // FIXME: The 0-polynomial would actually do.
            error ("requesting odd symmetry, but max degree is [0/0]");
        }

        p.symmetry = symmetry;
        q.symmetry = SymmetryEven;
    }

    p.deg = Target::deg_poly (n, p.symmetry);
    q.deg = Target::deg_poly (m, q.symmetry);
    p.n_non0 = Target::n_nonzero_poly_coeff (p.deg, p.symmetry);
    q.n_non0 = Target::n_nonzero_poly_coeff (q.deg, q.symmetry);
}

auto Approx::Target::xs_start() const -> Fs
{
    assert (n >= 0 && m >= 0);
    set_poly_descriptors();

    if (0 & symmetry)
    {
        auto xs = xs_start (n + m + 2, SymmetryNone);
        out ("========= xs vanilla ==========\n");
        for (auto &x : xs)
            out ("x: %f\n", (double) x);

        out ("========= xs 1 + vanilla ==========\n");
        xs = xs_start (n + m + 2 + 1, SymmetryNone);
        for (auto &x : xs)
            out ("x: %f\n", (double) x);

        out ("========= xs %d ==========\n", symmetry);
        xs = xs_start (p.n_non0 + q.n_non0, symmetry);
        for (auto &x : xs)
            out ("x: %f\n", (double) x);
    }

    //auto xs = xs_start (n + m + 2, symmetry);
    auto xs = xs_start (n + m + 2 + !!symmetry, SymmetryNone);
    return xs;
}

auto Approx::Target::xs_start (int nn, Symmetry sym) const -> Fs
{
    const F pi = F::pi();
    Fs xs;

    if (! sym)
    {
        F s = (b - a) >> 1;
        F t = pi / (nn - 1);
        xs.push_back (a);
        for (int i = nn - 2; i > 0; --i)
            xs.push_back (a + s * (cos (t * i) + 1));
    }
    else if (sym == SymmetryEven)
    {
        F s = b;
        F t = (pi >> 1) / (nn - 1);
        xs.push_back (F{0});
        for (int i = nn - 2; i > 0; --i)
        {
            xs.push_back (s * cos (t * i));
            //printf ("cos pi * %.4f\n", (double) (t*i/pi));
        }
    }
    else if (sym == SymmetryOdd && nn > 1)
    {
        F s = b;
        F t = (pi >> 1) / (2 * nn - 1);
        for (int i = 2 * nn - 2; i > 0; i -= 2)
        {
            xs.push_back (s * cos (t * (i)));
            //printf ("%d: cos pi * %.4f\n", i, (double) (t*i/pi));
        }
    }

    xs.push_back (b);
/*
    info (1, "xs[%d=%d], sym = %d", (int) xs.size(), nn, sym);
    for (int i = 0; i < (int) xs.size(); ++i)
        printf ("x[%d] = %.4f\n", i, (double) xs[i]);

    printf ("[");
    for (int i = 0; i < (int) xs.size(); ++i)
        printf ("%s%.4f", i ? ", " : "", (double) xs[i]);
        printf ("]\n");
*/

    if ((int) xs.size() != nn)
        error ("size=%d != nn=%d", (int) xs.size(), nn);

    return xs;
}

// Return a std::vector of  (x, f(x))  pairs.
auto Approx::get_values (const Fs& xs) const -> Values
{
    Values vs;
    for (const auto& x : xs)
    {
        F y = target.f(x);
        if (! y.number_p())
            errorX ("x = %e maps to Nan (x = %RNa)", (double) x, x.mpfr());

        vs.push_back (Value { x, y });
    }

    return vs;
}

bool Approx::Target::check_symmetry() const
{
    if (! symmetry)
        return true;

    if (abscmp (a + b, 1e-6)  > 0)
        error ("expect a = -b for symmetric target, is a=%f, b=%f",
               (double) a, (double) b);

    F eps { 1e-10 };
    for (int i = 0; i < 30; ++i)
    {
        F x = random();
        F ym = f (-x);
        F yp = f (x);
        F diff = symmetry == SymmetryOdd ? ym + yp : ym - yp;
        if (abs (diff) > eps)
            error ("broken symmetry %d at x=%.6f: f(-x)=%.6f, f(x)=%.6f",
                   symmetry, (double) x, (double) ym, (double) yp);
    }

    return true;
}

void Approx::remez()
{
    const bool rational = result.is_rational;
    assert (rational == !! target.m);

    if (Approx::verbose.level > 2)
        for (const auto& x : xs_remez)
            cout << x << endl;

    target.set_poly_descriptors();

    info (VERB.symm,
          "sym %d -> [%d/%d]", target.symmetry, target.p.deg, target.q.deg);
    info (VERB.non0, "non-0 = %d, %d", target.p.n_non0, target.q.n_non0);
    target.check_symmetry();

    Values vs = get_values (xs_remez);

    const F& E = rational ? result.E : F::Zero;
    if (rational)
        info (VERB.E, "using E = %.5e", (double) E);

    FMatrix<F> m = Remez::matrix (vs, target, E);
    if (! m.matrix_p())
        error ("%s", "no matrix");

    if (Approx::verbose.level > 2)
        cout << m << endl;

    FMatrix<F>::push_verbose (std::max (0, Approx::verbose.level - 1));
    FVector<F> v = Remez::vector (vs, target);
    if (! v.vector_p())
        error ("not a vector: v[] = %s", _tostr(v).c_str());
    FVector<F> w = m.solve (v, Algo::QR_vectors);
    if (! w.vector_p())
    {
        //Format<Float>::format = " %RNe ";
        info (1, "m[][] =\n%s", _tostr(m).c_str());

        error ("not a vector: w[] = %s", _tostr(w).c_str());
    }
    FMatrix<F>::pop_verbose();

    if (VERB.E)
    {
        cout << w << endl;
        if (rational)
            cout << "E_old = " << (double) result.E << endl;
        cout << "E     = " << (double) w[0] << endl;
    }

    assert (w.dim() == target.p.n_non0 + target.q.n_non0);

    // Because &w[i] is not the same like &w[0] + i.  The former is only
    // valid if w[i] is valid, which might not be the case when there are
    // no non-trivial q's.
    const F *w0 = &w[0];
    result.E = w[0];
    result.poly.p = Remez::poly_p (w0 + 1, target.p);
    result.poly.q = Remez::poly_q (w0 + 1 + target.p.n_non0, target.q);

    if (rational)
        result.check_denominator();

    result.have |= HavePoly;

    double bits_E = -std::log2 ((double) result.E.abs());
    int bits_F = F::GetPrecision (2);
    if (! std::isinf (bits_E) && bits_E > bits_F - 2)
        warning ("E = %.2f bits, Float = %d bits: consider increasing"
                 " Float precision", bits_E, bits_F);

    if (VERB.pq)
    {
        cout << "p(x) = " << result.poly.p << endl;
        if (rational)
            cout << "q(x) = " << result.poly.q << endl;
    }
}


std::ostream& operator << (std::ostream& ost, const Approx::Value& v)
{
    return ost << "(" << v.x << ";" << v.y << ")";
}


bool Approx::Quality::sign_flip (const Tree::Node *a, const Tree::Node *b)
{
    assert (a && b);
    const F& ya = a->t_.y;
    const F& yb = b->t_.y;
    assert (!ya.is_zero() && !yb.is_zero());
    return sign_flip (ya, yb);
}

auto Approx::delta (const F& x, int *pdenom_sign) const -> F
{
    F fx = target.f (x);
    F rx = result (x, pdenom_sign);

    if (! fx.number_p())
        error ("f(%e) = %f is not a number", (double) x, (double) fx);

    F rx_fx = rx - fx;

    if (target.error_relative() && fx.is_zero())
        errorX ("at x = %RNf: relative error requires target function has no"
                " zeroes", x.mpfr());
    return target.error_relative() ? rx_fx / fx : rx_fx;
}

Approx::Result::~Result ()
{
    delete quality;
    quality = nullptr;
}

Approx::Result::Result (const Target& t) : target(t), is_rational(t.m)
{
    if (t.n < 0)
        error ("Target: n = %d: degree of numerator must be >= 0", t.n);

    if (t.m < 0)
        error ("Target: m = %d: degree of denominator must be >= 0", t.m);

    if (t.a >= t.b)
        error ("Target: [a, b] must be an interval, but is [%e, %e]",
               (double) t.a, (double) t.b);
}

auto Approx::Result::operator () (const F& x0, int *psign) const -> F
{
    F x = x0.saturate (target.a, target.b);
    //if (poly.p.deg() != target.n)
    //    error ("deg(p) = %d != %d\n", poly.p.deg(), target.n);
    assert (poly.p.deg() <= target.n);
    assert (poly.q.deg() <= target.m);

    if (! is_rational)
    {
        if (psign) *psign = SignUnknown;
        return poly.p(x);
    }
    else
    {
        F q = poly.q(x);
        if (psign) *psign = q.sgn();
        return poly.p(x) / q;
    }
}

#define NEED(f, hs...)                                                      \
    do {                                                                    \
        for (Approx::Have _h : { hs })                                      \
            if (! ((f) & _h))                                               \
                error ("need have 0x%x for: %s", _h, __PRETTY_FUNCTION__);  \
    } while (0)

Float Approx::Result::delta_height (const Result& r, bool relative) const
{
    NEED (have, HavePoly);

    F dh = (poly.p - r.poly.p).height();
    if (is_rational)
        dh = dh.max ((poly.q - r.poly.q).height());

    if (relative)
    {
        F h = r.poly.p.height();
        if (is_rational)
            h = h.max (r.poly.q.height());
        dh /= h;
    }

    return dh;
}

// Easy search for zero of the denominator.  In general, everything is fine
// and the denominator won't have zeros.
void Approx::Result::check_denominator() const
{
    double step = 1.0 / 32;
    int qa_sgn = poly.q (target.a) . sgn();
    const auto &style = Poly<Float>::Style::PythonHorner;
    for (double t = 0; t <= 1.0; t += step)
    {
        F x = F{t}.linear (target.a, target.b);
        F y = poly.q (x);

        if (y.sgn() != qa_sgn)
        {
            cout << "p = "; poly.p.print (cout, style, "x", 0); cout << endl;
            cout << "q = "; poly.q.print (cout, style, "x", 0); cout << endl;
            error ("denominator has a zero near %.5f", (double) x);
        }
    }
}


// Add x, y = delta f(x) to the tree, keeping track of sign groups so far.
// Except in corner cases, this should only be used via  add (const F&) below.
auto Approx::Quality::add (const Value& v) -> const Tree::Node*
{
    if (v.y.is_zero())
    {
        tree_delta0.add (v, true);
        return nullptr;
    }

    bool changed;
    const Tree::Node *n = tree.add (v, true, &changed);
    if (!changed)
        return nullptr;

    const Tree::Node *l = n->prev();
    const Tree::Node *r = n->next();
    assert (!l || l->t_.x < n->t_.x);
    assert (!r || r->t_.x > n->t_.x);

    n_groups += 0?0
        : !l && !r ? 1
        : l && !r  ? sign_flip (l, n)
        : r && !l  ? sign_flip (n, r)
        : sign_flip (l, n) && sign_flip (n, r) ? 2
        : 0;

    return n;
}

// Add x, delta f(x) to the tree and keep track of the sign groups so far.
auto Approx::Quality::add (const F& x) -> const Tree::Node*
{
    int sgn_den_y;
    Float y = approx->delta (x, & sgn_den_y);

    return add (Value { x, y, sgn_den_y });

    //cout << "groups = " << n_groups << " Tree #" << tree.size() << ": ";
    //tree.dump (cout);
    //cout << endl;
}

// If we have not enough sign groups, try flipping the sign of a value
// at the interval endpoint, so that "phantom value" pops a new sign group.
void Approx::Quality::try_flip (Tree::Node *n)
{
    assert (! n->next() || ! n->prev());

    const Tree::Node *neigh = n->next() ? n->next() : n->prev();
    F &y = n->value().y;
    const F &w = neigh->value().y;

    if (y.sgn()
        && y.sgn() == w.sgn()
        // y is "very" small compared to w.
        && abscmp (y << 16, w) < 0)
    {
        // y is "almost" zero, which might be the reason for why we did not
        // find enough groups.  Flipping the sign of y pops one more group,
        // which might be enough to get Remez started.  This might occur
        // in symmetric cases.
        y = -y;
        n_groups += 1;
        const F& x = n->value().x;
        const F& a = approx->target.a;
        const F& b = approx->target.b;
        infoX (1, "flipped y at x%s = %.20f to %.20RNe",
               x == a ? " = a" : x == b ? " = b" : "", (double) x, y.mpfr());
    }
}

// Random in (a,b).
auto Approx::Target::random() const -> F
{
    assert (a < b);
    F f;
    do
    {
        double x = random_range<double> ((double) a, (double) b);
        f = (F) x;
    } while (f <= a || f >= b);

    return f;
}

// Random in (a,b) that's cos-distributed (more likely at interval ends).
auto Approx::Target::random_cos() const -> F
{
    assert (a < b);
    F f;
    do
    {
        double x = std::cos (random_range<double> (0, _pi));
        f = ((b - a) * F{x} + b + a) >> 1;
    } while (f <= a || f >= b);

    return f;
}

// Init .tree with enough values so we can begin grouping them by their signs
// and are getting as much such groups as we need.
auto Approx::Quality::init (int flags) -> Failure
{
    const Target& target = approx->target;

    int groups_needed = 2 + target.n + target.m;
    if (target.symmetry)
    {
        info (VERB.symm, "allow +1 group due to symmetry %d", target.symmetry);
        ++ groups_needed;
    }

    double d = 0.0000001;
    bool have_a = add (target.a);
    bool have_b = add (target.b);
    add (F{0+d}.linear (target.a, target.b));
    add (F{1-d}.linear (target.a, target.b));

    for (const F& x : approx->xs_remez)
        add (x);

    for (const F& x : target.xs_start (2 * groups_needed - 1, target.symmetry))
        add (x);

    for (int i = 0; n_groups < groups_needed && i < 10 * groups_needed; ++i)
        add (target.random_cos());

    if (tree.empty())
    {
        info (1, "found exact formula at %d points", tree_delta0.size());
        return Failure::NoDeltas;
    }

    // Last resort. Sometimes, in (symmetric?) cases, we have too few groups
    // because delta at the end of the target interval is zero.
    // Try adding phantom values.

    if (n_groups < groups_needed
        && ! have_a)
    {
        const F &ya = tree.first()->value().y;
        add (Value (target.a, -ya, SignUnknown));
    }

    if (n_groups < groups_needed
        && ! have_b)
    {
        const F &yb = tree.last()->value().y;
        add (Value (target.b, -yb, SignUnknown));
    }

    if (n_groups < groups_needed && tree.size() > 3)
        try_flip (tree.first());

    if (n_groups < groups_needed && tree.size() > 3)
        try_flip (tree.last());

    // Nothing worked, we are about to fail.

    if (n_groups < groups_needed)
    {
        const Value& first = tree.first()->t_;
        const Value& last  = tree.last()->t_;
        const char *is_a = first.x == target.a ? " = a" : " != a";
        const char *is_b = last.x  == target.b ? " = b" : " != b";
        cout << "Tree.first.x = " << first.x << is_a << endl;
        cout << "Tree.last.x  = " << last.x  << is_b << endl;
    }

    if (n_groups > groups_needed
        && n_groups <= groups_needed + allow_more_groups)
    {
        info (1, "must ditch %d groups", n_groups - groups_needed);
    }
    else if (n_groups != groups_needed)
    {
        //for (auto& g : groups)
        //    cout << "worst[" << g.id << "] " << g.worst.y << endl;

        cout << "Tree #" << tree.size() << endl;
        for (const Tree::Node *n = tree.first(); n; n = n->next())
            cout << n->t_ << endl;

        if (flags & Flag::OnBadGroupNumber_write_deltas)
            approx->write_deltas ("deltas-fail.data", 777);

        if (flags & Flag::OnBadGroupNumber_error)
            error ("need %d groups for degree %d/%d, found %d",
                groups_needed, target.n, target.m, n_groups);
        return n_groups < groups_needed ? TooFewMaxima : TooManyMaxima;
    }

    assert (n_groups >= groups_needed);
    assert (n_groups <= groups_needed + allow_more_groups);

    return Success;
}

bool Approx::Target::operator == (const Target& t) const
{
    return a == t.a && b == t.b && m == t.m && n == t.n && f == t.f;
}

auto Approx::rate() -> const Result::Rating*
{
    return result.rate (this);
}

auto Approx::Result::find_worst (int flags) -> Failure
{
    Approx *approx = (Approx*) ((char*) this - offsetof (Approx, result));
    assert (target == approx->target);
    Approx::Quality *q = new Approx::Quality (approx);
    quality = q;

    have |= HaveQuality;

    rating.failure = failure = q->init (flags);
    if (failure)
        return failure;

    q->make_groups();
    q->ditch_excess_groups();

    have |= HaveQualityGroups;

    if (Approx::verbose.level > 2)
        for (auto& g : q->groups)
        {
            out ("lo[%d] = %f, %e\n", g.id, (double) g.lo->t_.x, (double) g.lo->t_.y);
            out ("hi[%d] = %f, %e\n", g.id, (double) g.hi->t_.x, (double) g.hi->t_.y);
            //out ("lo[%d] = %.66RNf, %.66RNe\n", g.id, g.lo->t_.x, g.lo->t_.y);
            //out ("hi[%d] = %.66RNf, %.66RNe\n", g.id, g.hi->t_.x, g.hi->t_.y);
        }

    //for (auto& g : q->groups)
    //    cout << g.worst.x << ": " << (double) g.worst.y <<  endl;

    if (flags & Flag::FindZeros)
    {
        q->find_zeros();
        have |= HaveQualityZeros;
    }

    for (int i = 0; i < q->n_groups; ++i)
    {
        for (int j = 0; j < 10; ++j)
        {
            Failure fail = q->groups[i].find_worse();
            if (fail == Failure::Void)
                // Nothing found.
                break;
            else if (fail != Failure::Success)
                return fail;
        }
    }

    if (Approx::verbose.level)
    {
        Float::Format::push ("% .15RNe");
        for (auto& g : q->groups)
            cout << "worst[" << g.id << "] " << g.worst.y
                 << " (x=" << g.worst.x << ")" << endl;
        Float::Format::pop ();
    }

    assert (xs.empty());
    for (auto& g : q->groups)
    {
        values.push_back (g.worst);
        xs.push_back (g.worst.x);
    }

    have |= HaveWorst;

    return Failure::Success;
}

auto Approx::Result::rate (const Approx *approx) -> const Rating*
{
    assert (approx == (Approx*) ((char*) this - offsetof (Approx, result)));
    assert (target == approx->target);
    assert (quality);
    NEED (have, HaveQuality);

    rating.set (quality, *this);

    have |= HaveRating;

    return & rating;
}


void Approx::Quality::SignGroup::set_worst (const Tree::Node *node)
{
    const Value& w = node->t_;

    if (0)
    {
        Value dv { w.x - worst.x, w.y - worst.y };
        cout << "worst["<< id <<"] " << worst.y << " -> "<< w.y << endl;
        cout << "delta["<< id <<"] " << dv << endl;
    }
    worst = w;
    worst_node = node;
}


bool Approx::Quality::SignGroup::maybe_set_worst (const Tree::Node *node)
{
    if (is_worse (node->t_.y))
    {
        set_worst (node);
        return true;
    }
    return false;
}

// Make sure WORST is surrounded by nodes of same sign.
// This might not be possible, e.g. when we has to add phantom values.
auto Approx::Quality::SignGroup::worst_fix_neighbours() -> Failure
{
    int n_loops = 0;
    bool changed_worst;

    do
    {
        changed_worst = false;

        for (Tree::Side side : { Tree::Left, Tree::Right })
        {
            assert (worst_node);

            const Tree::Node *neighbour = worst_node->neighbour(side);
            if (!neighbour)
                continue;

            bool changed = false;

            Value v = neighbour->t_;
            while (Quality::sign_flip (worst.y, v.y) || v.y.is_zero())
            {
                //cout << "worst.y = " << worst.y << endl;
                //cout << "v.y     = " << v.y << endl;
                v.x = (worst.x + v.x) >> 1;
                if (v.x == worst.x)
                    break;
                v.y = approx->delta (v.x, & v.sign_denom_y);

                if (v.has_denom_zero (worst))
                {
                    double mi = (double) v.x.min (worst.x);
                    double ma = (double) v.x.max (worst.x);
                    error ("approximation denominator has a zero in "
                           "[%f, %f], id = %d", mi, ma, id);
                    return Failure::DenomZero;
                }
                changed = true;
                assert (! worst.y.is_zero());
                assert (worst.x != v.x);
            }

            if (changed)
            {
                bool added;
                Tree::Node *neu = quality->tree.add (v, true, &added);
                if (0 && !added)
                {
                    warningX ("%.20f, %RNe", (double) worst.x, worst.y.mpfr());
                    warningX ("%.20f, %RNe", (double) v.x, v.y.mpfr());
                }
                assert (neu);
                if (!added)
                    continue;

                /*cout << "v.x:     " << v.x << endl;
                  cout << "worst.x: " << worst.x << ", " << worst_node->t_.x << endl;
                  cout << "neigh.x: " << neighbour->t_.x << endl;
                  cout << "neu.x  : " << neu->t_.x << endl;*/
                assert (neu->neighbour((Tree::Side) !side) == worst_node);
                assert (neu == worst_node->neighbour(side));

                //cout << "add: " << v << endl;

                changed_worst = maybe_set_worst (neu);
                if (changed_worst)
                    break;
            }

            // Avoid skew, i.e. that worst.x leans too much to some side.

            const Tree::Node *l = worst_node->prev();
            const Tree::Node *r = worst_node->next();
            if (!changed_worst
                && l && l->t_.y.sgn() == sign
                && r && r->t_.y.sgn() == sign)
            {
                double skew = 0.05;
                Value v;
                double t = (double) worst.x.where_in (l->t_.x, r->t_.x);
                v.x = 0?0
                    : t < 2*skew   ? Float{skew}.linear (l->t_.x, r->t_.x)
                    : t > 1-2*skew ? Float{1-skew}.linear (l->t_.x, r->t_.x)
                    : v.x;

                if (v.x.number_p())
                {
                    bool added;
                    v.y = approx->delta (v.x, & v.sign_denom_y);
                    Tree::Node *neu = quality->tree.add (v, true, &added);
                    if (added)
                        changed_worst = maybe_set_worst (neu);
                }
            }
        } // for Tree::Side
    } while (changed_worst && ++n_loops < 100);

    return Failure::Success;
}

auto Approx::Quality::SignGroup::find_worse() -> Failure
{
    //cout << "-- " << id << endl;
    Failure fail = worst_fix_neighbours();
    if (fail)
        return fail;

    const Tree::Node *left  = worst_node->prev();
    const Tree::Node *right = worst_node->next();
    if (!left || !right)
        return Failure::Void;
    const Value& v0 = left->t_;
    const Value& v1 = worst;
    const Value& v2 = right->t_;

    if (v0.has_denom_zero (v1))
    {
        error ("approximation denominator has a zero in [%f, %f], "
               "id = %d", (double) v0.x, (double) v1.x, id);
        return Failure::DenomZero;
    }
    if (v1.has_denom_zero (v2))
    {
        error ("approximation denominator has a zero in [%f, %f], "
               "id = %d", (double) v1.x, (double) v2.x, id);
        return Failure::DenomZero;
    }

    if (! (v0.y.sgn() == v1.y.sgn() && v1.y.sgn() == v2.y.sgn()))
    {
        Float::format.str = "% .10RNe ";
        cout << "v0 = " << v0 << endl;
        cout << "v1 = " << v1 << endl;
        cout << "v2 = " << v2 << endl;
        // Run into assert() below.
    }
    assert (v0.y.sgn() == v1.y.sgn() && v1.y.sgn() == v2.y.sgn());
    assert (v0.x < v1.x && v1.x < v2.x);

    /*cout << v0 << endl;
    cout << v1 << endl;
    cout << v2 << endl;*/

    for (int how = 0; how < 3; ++how)
    {
        if (find_worse (how, Values { v0, v1, v2 }))
            return Failure::Success;
    }

    return Failure::Void;
}

bool Approx::Quality::SignGroup::find_worse (int how, const Values& v)
{
    assert (v.size() == 3);

    Value w;

    switch (how)
    {
        default:
            assert (0);

        case 0:
        {
            // P is a parabola passing through V[0...2].  It's derivative
            // is P'(x) = P_1 + 2 * P_2 * x  with zero  x' = -P_1 / (2 * P_2).
            Poly<F> p = Value::polynomial (v);
            if (p.deg() < 2)
                return false;

            w.x = -p[1] / (p[2] << 1);
            break;
        }

        case 1:
        case 2:
        {
            // Determine crossing point of tangents at adjacent points.
            int shift = 13;
            F h { std::ldexp (1.0, -1 - shift) };
            const F& x0 = v[how - 1].x, &x1 = v[how].x;
            const F& y0 = v[how - 1].y, &y1 = v[how].y;
            F m0 = (approx->delta (x0 + h) - approx->delta (x0 - h)) >> shift;
            F m1 = (approx->delta (x1 + h) - approx->delta (x1 - h)) >> shift;
            if (m0 * m1 > -1e-16)
                return false;
            w.x = ((m0*x0 - y0) - (m1*x1 - y1)) / (m0 - m1);
            info (1, "m = %g, %g\n", (double) m0 , (double) m1);
            assert (0);
            break;
        }
    }

    w.y = approx->delta (w.x, & w.sign_denom_y);

    if (v[0].x < w.x && w.x < v[2].x)
    {
        bool changed;
        Tree::Node *neu = quality->tree.add (w, true, &changed);
        if (!changed)
            return false;
        maybe_set_worst (neu);
    }
    else
    {
        /*info (how, "outside %g < %g < %g",
              (double) v[0].x, (double) w.x, (double) v[2].x);*/
        assert (isfinite (w.x));
        assert (isfinite (w.y));
        return false;
        error ("todo");
    }

    return true;
}


void Approx::Quality::make_groups()
{
    groups.resize (n_groups);
    for (int i = 0; i < n_groups; ++i)
    {
        SignGroup &g = groups[i];
        g.id = i;
        g.prev = i == 0 ? nullptr : &groups[i-1];
        g.next = i == n_groups - 1 ? nullptr : &groups[i+1];
        g.quality = this;
        g.approx = approx;
    }

    int prev_sign = 0;
    SignGroup *g = nullptr;

    for (const Tree::Node *n = tree.first(); n; n = n->next())
    {
        const F& y = n->t_.y;
        bool is_worse = true;
        if (y.sgn() != prev_sign)
        {
            assert (y.sgn());
            g = prev_sign ? g->next : &groups[0];
            g->sign = prev_sign = y.sgn();
            g->lo = n;
            if (g->prev)
                g->prev->hi = g->lo->prev();
        }
        else
            is_worse = g->is_worse (y);

        if (is_worse)
        {
            g->worst_node = n;
            g->worst = n->t_;
        }
        //cout << n->t_.x << " in " << g->id <<  endl;
    }
    assert (! g->next);
    g->hi = tree.last();

    for (auto& g : groups)
        assert (g.worst_node);
}

void Approx::Quality::ditch_excess_groups ()
{
    const Target& target = approx->target;
    int groups_needed = 2 + target.n + target.m + !!target.symmetry;

    if (n_groups == groups_needed)
        return;

    assert (n_groups > groups_needed);
    assert (n_groups <= groups_needed + allow_more_groups);
    info (1, "have %d groups, need %d", n_groups, groups_needed);

    while (n_groups > groups_needed)
    {
        info (1, "ditching group");
        assert (n_groups == (int) groups.size());
        const SignGroup& g0 = groups[0];
        const SignGroup& g1 = groups[n_groups - 1];
        info (1, "g0.worst = %f", (double) g0.worst.y);
        info (1, "g1.worst = %f", (double) g1.worst.y);
        if (g0.worst.y.abscmp (g1.worst.y) < 0)
            // std::vector does not have pop_front...
            groups.erase (groups.begin());
        else
            groups.pop_back();
        --n_groups;
    }
}

void Approx::Quality::find_zeros()
{
    for (SignGroup& g : groups)
        if (g.prev)
        {
            int more = 8 * (g.id == 1 || ! g.next);
            for (int i = 0; i < 8 + more; ++i)
                if (g.find_zero())
                    break;
        }
}


bool Approx::Quality::SignGroup::find_zero()
{
    const Value& v0 = prev->hi->t_;
    const Value& v1 = lo->t_;

    if (Approx::verbose.level > 2)
    {
        cout << id << ": " << v0.x << " -> " << (double) v0.y << endl;
        cout << id << ": " << v1.x << " -> " << (double) v1.y << endl;
    }

    assert (sign_flip (v0.y, v1.y));
    int sign_denom_y;

    F x = (v0.x * v1.y - v1.x * v0.y) / (v1.y - v0.y);
    F y = approx->delta (x, &sign_denom_y);
    if (Approx::verbose.level > 2)
        cout << id << "= " << x << " -> " << (double) y << endl;
    assert (v0.x < x && x < v1.x);

    const Tree::Node *n
        = quality->tree.add (zero = Value { x, y, sign_denom_y });
    assert (n);

    if (y.sgn() == prev->sign)
        prev->hi = node_zero = n;
    else
        lo = node_zero = n;

    return y.is_zero();
}


auto Approx::Quality::SignGroup::width_x_left() const -> const F&
{
    assert (id == 0 || node_zero);
    return prev ? zero.x : lo->t_.x;
}

auto Approx::Quality::SignGroup::width_x_right() const -> const F&
{
    assert (id == 0 || node_zero);
    return next ? next->zero.x : hi->t_.x;
}

auto Approx::Quality::SignGroup::width() const -> F
{
    assert (id == 0 || node_zero);
    return width_x_right() - width_x_left();
}

void Approx::Result::adjust_worst_place0()
{
    NEED (have, HaveQuality, HaveQualityGroups, HaveQualityZeros);

    FVector<F> a = quality->adjust_worst_place0();
    F new_x0 { quality->groups[0].width_x_left() };

    place0.E = a[quality->n_groups];
    place0.E = place0.E.setsign (quality->groups[0].sign);
    for (auto& g : quality->groups)
    {
        const F& x0 = g.width_x_left();
        const F& x1 = g.width_x_right();
        // Relative position in [0,1] where WORST.x is located.
        F t = g.worst.x.where_in (x0, x1);
        assert (t >= 0 && t <= 1);
        F new_width = a[g.id] * (x1 - x0);
        F new_x1 = new_x0 + new_width;
        F new_worst_x = new_x0 + t * new_width;
        place0.xs.push_back (new_worst_x);
        cout << "new_worst[" << g.id << "] = " << (double) new_worst_x << endl;
        new_x0 = new_x1;
    }
}

auto Approx::Quality::place0_matrix (const Target&) const -> FMatrix<F>
{
    int dim = 1 + n_groups;
    FMatrix<F> m = FMatrix<F>::make (dim, dim, F(0));
    for (int i = 0; i < n_groups; ++i)
    {
        const F& y = groups[i].worst.y;
        m[i][i] = y.abs();
        m[i][dim - 1] = F(-1);
        m[dim - 1][i] = groups[i].width();
    }
    cout << m << endl;
    return m;
}

auto Approx::Quality::place0_vector (const Target& t) const -> FVector<F>
{
    return FVector<F>::e_k (1 + n_groups, n_groups, t.b - t.a);
}

auto Approx::Quality::adjust_worst_place0() const -> FVector<F>
{
    assert (groups[1].node_zero);

    for (auto& g : groups)
    {
        if (g.prev)
            cout << "zero[" << g.id << "]  = " << DValue (g.node_zero->t_)<< endl;
        cout << "worst[" << g.id << "] = " << DValue (g.worst) << endl;
    }
    FMatrix<F> m = place0_matrix (approx->target);
    FVector<F> v = place0_vector (approx->target);
    FVector<F> a = m.solve (v, Algo::QR_vectors);

    cout << a << endl;

    return a;
}

double Approx::xacos (const F& x) const
{
    assert (result.quality);
    return result.quality->xacos (x);
}

double Approx::Quality::xacos (const F& x0) const
{
    assert (x0.in_range (approx->target.a, approx->target.b));

    // Find the two groups between whose extrama X is located.

    const SignGroup *lo = & groups[0];
    const SignGroup *hi = & groups[n_groups - 1];

    assert (lo != hi && lo->next != hi);

    // There are cases when the maximum of the outermost group is not
    // at the interval border; rectify that now.
    const F x = x0.saturate (lo->worst.x, hi->worst.x);

    while (lo->next != hi)
    {
        if (0 && Approx::verbose.level >= 2)
            out ("[%d,%d] ", lo->id, hi->id);
        const SignGroup *g = & groups[(lo->id + hi->id) / 2];
        assert (g != lo && g != hi);
        if (x > g->worst.x)
            lo = g;
        else
            hi = g;
    }
    if (0 && Approx::verbose.level >= 2)
        out ("[%d,%d]\n", lo->id, hi->id);

    if (! (x >= lo->worst.x && x <= hi->worst.x))
    {
        for (auto& g : groups)
            out ("[%d] = %f\n", g.id, (double) g.worst.x);

        warning ("need %g <= %g <= %g", (double) lo->worst.x,
                 (double) x, (double) hi->worst.x);
    }
    assert (x >= lo->worst.x && x <= hi->worst.x);

    F d = approx->delta (x);
    bool in_hi = d.sgn() == hi->sign;
    // Normalize D so it is in the range [-1,1] appropriate for acos().
    double y = (double) (d / abs (in_hi ? hi->worst.y : lo->worst.y));
    double acos_y = std::acos (y < -1 ? -1 : y > 1 ? 1 : y);
    // Extend to 2*pi, the full period of cos().  This is needed if we are
    // in [pi,2pi] where cos rises again.
    if (lo->sign < 0) acos_y = 2 * _pi - acos_y;
    //out ("xacos(%f).0 = %f\n", y, acos_y);
    //if (groups[0].sign < 0) acos_y -= pi;

    bool neg0 = groups[0].sign < 0;
    // Determine the appropriate branch of acos.  We have N_GROUPS groups,
    // thus the possible range of values is  [0, pi * (N_GROUPS - 1)].
    // Then formalize to [0,1] and finally to [-1,1].

    int n_2pi = (lo->id + neg0) / 2;
    int n_pi = 2 * n_2pi - neg0;
    acos_y += _pi * n_pi;
    acos_y /= _pi * (n_groups - 1);

    return -std::cos (_pi * acos_y);

    acos_y = 2 * acos_y - 1;
    //out ("xacos(%f).2 = %f\n", y, acos_y);
    acos_y = std::sin (acos_y * _pi / 2);
    return acos_y;
}


double Approx::Result::bits_for_print() const
{
    NEED (have, HaveRating);

    double prec2 = rating.prec2;

    // Evaluating a polynomial of degree N with Horner takes
    // 2 * N - 1 operations.  Assume each operation loses 0.5 bits
    // of precision.

    return prec2 + target.n + target.m;
}


void Approx::Result::Rating::set (const Quality* q, Result& result)
{
    int i = 0;
    for (auto& g : q->groups)
    {
        if (!i || abscmp (g.worst.y, delta.worst.y) > 0) delta.worst = g.worst;
        if (!i || abscmp (g.worst.y, delta.best.y)  < 0) delta.best = g.worst;
        i++;
    }
    F worst = delta.worst.y.abs();
    F best  = delta.best.y.abs();
    quality = worst / best;
    double dworst = (double) worst;
    quality_p10 = -log10 ((double) (quality - 1_Float));
    prec2 = -std::log2 (dworst);
    prec10 = -std::log10 (dworst);

    height_p = result.poly.p.height();
    height_q = result.poly.q.height();

    height = result.is_rational
        ? std::max (height_p, height_q)
        : height_p;
}

void Approx::Result::Rating::print() const
{
    cout << "|delta| (worst): " << (double) delta.worst.y.abs() << endl;
    cout << "|delta| (best) : " << (double) delta.best.y.abs() << endl;
    cout << "prec bits    : " << prec2 << endl;
    cout << "prec digits  : " << prec10 << endl;
    cout << "quality      : " << quality << endl;
    cout << "qualit.p10   : " << quality_p10 << endl;
}

void Approx::write_deltas (FILE *fout, int n, const char *label) const
{
    char lab[100];
    assert (n > 1);

    if (!label)
    {
        const char *s_arel = target.error_relative() ? "relative" : "absolute";
        if (result.target.m == 0)
            sprintf (lab, "MiniMax degree %d polynomial %s error",
                     result.target.n, s_arel);
        else
            sprintf (lab, "MiniMax rational [%d/%d] %s error",
                     result.target.n, result.target.m, s_arel);
        label = lab;
    }

    int deg = result.target.n + result.target.m;
    int prec2 = Float::GetPrecision (2);
    int prec2_for_print = std::isnan (result.rating.prec2)
        ? prec2
        : (int) result.rating.prec2 + 2*10 + deg;
    push_precision<Float> (std::min (prec2, prec2_for_print), 2);

    fprintf (fout, "\"x\" \"%s\" \"x\" \"predict\" \"x\" \"max\"\n", label);

    int m = n + (int) result.values.size();
    DValue* vs = new DValue[m];

    const F& a0 = target.a;
    const F& b0 = target.b;

    F y, x, dx = (b0 - a0) / (n - 1);
    for (int i = 0; i < n; ++i)
    {
        x = (a0 + dx * F(i)).saturate (a0, b0) ;//-1e-10;
        vs[i] = DValue (x, delta (x));
    }
    // Make sure we are drawing the maxima.
    if (flags & Flag::PlotMaxima)
        for (int i = n; i < m; ++i)
        {
            x = result.values[i - n].x;
            vs[i] = DValue (x, result.values[i - n].y);
        }

    std::sort (vs, vs + m);
    //const Quality *q = result.quality;
    for (int i = 0; i < m; ++i)
    {
        fprintf (fout, "% .5e % .5e", vs[i].x, vs[i].y);
        // x's used as input.
        if (i < (int) xs_remez.size())
        {
            x = xs_remez[i];
            DValue v (x, delta (x/*s_remez[i]*/));
            fprintf (fout, " % .5e % .5e", v.x, v.y);
        }
        // Determined to be the local extrema of delta().
        if (i < (int) result.values.size())
        {
            x = result.values[i].x;
            DValue v { x, result.values[i].y };
            fprintf (fout, " % .5e % .5e", v.x, v.y);
        }
        /*if (q
            && i < (int) q->groups.size() - 1
            && q->groups[i + 1].node_zero)
        {
            DValue v (q->groups[i + 1].node_zero->t_);
            fprintf (fout, " % .5e % .5e", v.x, v.y);
        }*/
        fprintf (fout, "\n");
    }

    delete[] vs;
    pop_precision<Float>();
}
void Approx::write_deltas (const char *fname, int n, const char *label) const
{
    if (!strstr (fname, ".data"))
        error ("suggest file name \"%s\" with \".data\" suffix", fname);

    FILE *fout = fopen (fname, "w");
    if (fout)
    {
        write_deltas (fout, n, label);
        fclose (fout);
    }
}

void Approx::write_deltas (FILE *fout, Func f, const AsinContext &ctx,
                           AsinContext::Func g, int n, const char* head) const
{
    auto delta = [] (const Float& x, Func f, const AsinContext &ctx,
                     AsinContext::Func g) -> Float
    {
        push_precision<Float> (200, 2);
        Float fx = f(x);
        pop_precision<Float> ();
        Float y = (g(x, ctx) - fx) / fx;
        return ctx.ulp ? y << ctx.ulp : y;
    };

    int deg = result.target.n + result.target.m;
    int prec2 = Float::GetPrecision (2);
    int prec2_for_print = std::isnan (result.rating.prec2)
        ? prec2
        : (int) result.rating.prec2 + 2*10 + deg;
        push_precision<Float> (std::min (prec2, prec2_for_print), 2);

    fprintf (fout, "\"x\" \"%s\" \"y\" \"bababa\"\n", head);

    const F a0 { -1.0 };
    const F b0 { +1.0 };
    Float eps = Float { 0.00001 };
    const Float extra[6] =
        {
            a0+eps, F{-0.5}-eps, F{-0.5}+eps, F{0.5}-eps, F{0.5}+eps, b0-eps
        };
    int m = n + (int) (sizeof (extra) / sizeof (*extra));

    DValue* vs = new DValue[m];

    Float x, dx = (b0 - a0) / (n - 1);
    for (int i = 0; i < m; ++i)
    {
        x = i < n
            ? (a0 + dx * F(i)).saturate (a0, b0)
            : extra[i - n];
        vs[i] = DValue (x, delta (x,f,ctx,g));
    }

    std::sort (vs, vs + m);

    for (int i = 0; i < m; ++i)
    {
        fprintf (fout, "% .5e % .5e\n", vs[i].x, vs[i].y);
    }

    delete[] vs;
    pop_precision<Float>();
}
void Approx::write_deltas_asin (FILE *fout, AsinContext& ctx, int n) const
{
    char head[200];
    sprintf (head, "%s relative error (degree [%d/%d] MiniMax%s)",
             ctx.asin_p ? "arcsin" : "arccos", target.n, target.m,
             ctx.ulp ? ", ULPs" : "");
    if (ctx.asin_p)
        write_deltas (fout, ::asin, ctx, ctx.as, n, head);
    else
        write_deltas (fout, ::acos, ctx, ctx.ac, n, head);
}
bool Approx::write_deltas_asin (const char *fname, int n, bool asin_p,
                                int ulp) const
{
    if (target.a > F{0} || target.b < F{0.5})
    {
        warning ("%s: target interval does not contain [0, 0.5]", __FUNCTION__);
        return false;
    }

    if (!strstr (fname, ".data"))
        error ("suggest file name \"%s\" with \".data\" suffix", fname);

    using AC = const AsinContext&;
    AsinContext ctx
        {
            .asin_p = asin_p,
            .ulp = ulp,
            .p = result.poly.p,
            .q = result.poly.q,
            .pi = Float::pi(),
            .pi2 = Float::pi() >> 1,
            // a(x) = arcsin(q) / q; q = sqrt(x/2), x in [0, 0.5]
            .a = [](const Float& x, AC _) -> Float
            {
                return _.p(x) / _.q(x);
            },
            .c = [](const Float& x, AC _) -> Float // arccos in [0.5, 1].
            {
                Float x1 = Float::One - x;
                return _.a(x1, _) * (x1 << 1).sqrt();
            },
            .as = [] (const Float& x, AC _) -> Float // arcsin in [-1, 1]
            {
                return x.abscmp (0.5) <= 0
                    ? x * _.a ((x * x) << 1, _)
                    : x.sgn() > 0 ? _.pi2 - _.c(x, _) : _.c(-x, _) - _.pi2;
            },
            .ac = [] (const Float& x, AC _) -> Float // arccos in [-1, 1]
            {
                return x.abscmp (0.5) <= 0
                    ? _.pi2 - _.as (x, _)
                    : x.sgn() > 0 ? _.c(x, _) : _.pi - _.c(-x, _);
            }
        };

    FILE *fout = fopen (fname, "w");
    if (fout)
    {
        write_deltas_asin (fout, ctx, n);
        fclose (fout);
    }
    return true;
}


void Approx::write_xacos (FILE *fout, int n, bool sub_x) const
{
    assert (n > 1);
    Poly<double> h = phase_poly (6, 17, sub_x);

    F x, dx = (target.b - target.a) / (n - 1);
    for (int i = 0; i < n; ++i)
    {
        x = target.a + dx * F(i);
        double px = 2.0 * i / (n - 1) - 1;
        double py = xacos (x) - px * sub_x;
        fprintf (fout, "%.16e %.16e %.16e\n", px, py, h (px));
    }
}
void Approx::write_xacos (const char *fname, int n, bool sub_x) const
{
    if (!strstr (fname, ".data"))
        error ("suggest file name \"%s\" with \".data\" suffix", fname);

    FILE *fout = fopen (fname, "w");
    if (fout)
    {
        write_xacos (fout, n, sub_x);
        fclose (fout);
    }
}


Poly<double> Approx::phase_poly (int n, int n_points, bool sub_x) const
{
    F x, dx = (target.b - target.a) / (n_points - 1);
    DValues vs;
    for (int i = 0; i < n_points; ++i)
    {
        x = target.a + dx * F(i);
        double px = 2.0 * i / (n_points - 1) - 1;
        double py = xacos (x) - px * sub_x;

        vs.push_back (DValue { px, py });
    }

    return best_polynomial (n, vs);
}

template<class A>
void _best_polynomial_row (int r, FVector<A>& row, FVector<A>& b,
                           const A& x, const A& y)
{
    b[r] = y;
    row[0] = 1;
    A xc = x;
    for (int c = 1; c < row.dim(); ++c, xc *= x)
        row[c] = xc;
}

namespace
{

template<class A, class V>
Poly<A> _best_polynomial (int n, const std::vector<V>& vs,
                          const std::vector<V>* fixed)
{
    int n_vs = (int) vs.size();
    int n_fixed = fixed ? (int) fixed->size() : 0;

    if (n_fixed)
        error ("todo: fixed");

    if (n_vs <= n)
        error ("need at least %d points for best_polynomial of degree %d, "
               "have %d", 1 + n, n, n_vs);

    FMatrix<A> m { Alloc (n_vs, n + 1) };
    FVector<A> b { Alloc (n_vs) };

    for (int i = 0; i < n_vs; ++i)
        _best_polynomial_row<A> (i, m[i], b, vs[i].x, vs[i].y);

    FMatrix<A> mT = m.T();
    FVector<A> a = (mT*m).solve (mT*b, Algo::QR_vectors);

    Poly<A> p;
    for (int i = 0; i < a.dim(); ++i)
        p.set (i, a[i]);

    return p;
}
}; // ::anon

Poly<double>
Approx::best_polynomial (int n, const DValues& vs, const DValues* fixed)
{
    return _best_polynomial<double, DValue> (n, vs, fixed);
}

Poly<Float>
Approx::best_polynomial (int n, const Values& vs, const Values* fixed)
{
    return _best_polynomial<Float, Value> (n, vs, fixed);
}
