#ifndef MP_APPROX_H
#define MP_APPROX_H

#include "mp-float.h"
#include "mp-poly.h"
#include "rb-tree.h"

#include <vector>
#include <cstdio>
#include <cmath>

#include "diagnostic.h"

class Approx
{
public:
    typedef Float F;
    typedef std::vector<F> Fs;
    struct Result;
    struct Target;
    typedef F (*Func)(const F&);
    struct Value;
    typedef std::vector<Value> Values;
    static constexpr int SignUnknown = -2;
    struct Value
    {
        F x, y;
        static Poly<F> polynomial (const Values&);
        // Order when we put Value's in a Tree.
        bool operator < (const Value& v) const { return x < v.x; };
        // May hold the sign of the denominator of y when the approximation
        // target is a rational function.
        int sign_denom_y = SignUnknown;
        Value (const F& x0, const F& y0) : x(x0), y(y0) {}
        Value (const F& x0, const F& y0, int s)
            : x(x0), y(y0), sign_denom_y(s) {}
        Value () {}
        bool has_denom_zero (const Value& v) const
        {
            return (sign_denom_y != SignUnknown && v.sign_denom_y != SignUnknown
                    && (sign_denom_y != v.sign_denom_y || sign_denom_y == 0));
        }
    };
    struct DValue
    {
        double x, y;
        DValue () {};
        DValue (const F& x, const F& y) : x((double) x), y((double) y) {};
        DValue (const Value& v) : DValue (v.x, v.y) {};
        bool operator < (const DValue& v) const { return x < v.x; };
    };
    typedef std::vector<DValue> DValues;
    enum Norm
    {
        NormQ0, NormQm, NormP0, NormPn
    };
    enum Flag
    {
        OnBadGroupNumber_error = 1 << 0,
        OnBadGroupNumber_write_deltas = 1 << 1,
        FindZeros = 1 << 2,
        AdjustZeros = 1 << 3,
        ErrorRelative = 1 << 4,
        PlotMaxima = 1 << 5,
        Null = 0
    };
    int flags = PlotMaxima;
    typedef RBTree<Value> Tree;
    struct Remez;
    struct Place0;
    struct Quality;
    enum Symmetry
    {
        SymmetryNone = 0, SymmetryOdd = 1, SymmetryEven = 2
    };
    struct Target
    {
        void set_degree (int n, int m, Symmetry sym)
        {
            this->n = n;
            this->m = m;
            this->symmetry = sym;
            set_poly_descriptors();
        }
        friend Approx;
        // Degrees of numerator / denominator polynomial degrees.
        // If m = 0 then it's approximation by polynomials, otherwise by
        // rational functions.
        int n = 0;
        int m = 0;
        int flags = Flag::PlotMaxima;
        Norm norm = NormQ0;
        Approx::Symmetry symmetry = SymmetryNone;
        // Approximate f(x) over [a,b].
        F a, b;
        Func f = nullptr;
        const char *f_str = "<undefined>";

        // Values used internally.
        struct PolynomialDescriptor
        {
            int deg;
            Symmetry symmetry = SymmetryNone;
            int n_non0;
        };
        mutable PolynomialDescriptor p;
        mutable PolynomialDescriptor q;
        void set_poly_descriptors() const;
        bool check_symmetry() const;

        F random() const;
        F random_cos() const;
        Fs xs_start () const;
        Fs xs_start (int, Symmetry) const;
        bool operator == (const Target& t) const;
        bool error_relative() const { return flags & Flag::ErrorRelative; }
        bool error_absolute() const { return ! error_relative(); }
        // Number of non-zero coefficients in polynomial of degree <= deg
        // and symmetry sym.
        static int n_nonzero_poly_coeff (int deg, Symmetry sym)
        {
            switch (sym)
            {
                case SymmetryNone: return deg + 1;
                case SymmetryOdd:  return (deg + 1) / 2;
                case SymmetryEven: return 1 + deg / 2;
                default: fatal ("unreachable");
            }
            return 0;
        }
        // Maximum possible degree for a polynomial of symmetry sym and
        // degree <= deg;
        static int deg_poly (int deg, Symmetry sym)
        {
            switch (sym)
            {
                case SymmetryNone: return deg;
                case SymmetryOdd:  return deg == 0 ? 0 : (deg - 1) | 1;
                case SymmetryEven: return deg & ~1;
                default: fatal ("unreachable");
            }
            return 0;
        }

    } target;
    enum Failure
    {
        Success = 0, TooFewMaxima, TooManyMaxima, NoDeltas, DenomZero,
        Void
    };
    enum Have
    {
        HavePoly = 1 << 0,
        HaveQuality = 1 << 1,
        HaveQualityGroups = 1 << 2,
        HaveQualityZeros = 1 << 3,
        HaveWorst = 1 << 4,
        HaveRating = 1 << 5
    };
    Fs xs_remez;
    struct Result
    {
        Result() = delete;
        Target target;
        struct { Poly<F> p, q; } poly;
        //__attribute__((deprecated))
        F E;
        bool is_rational;
        Values values;
        Fs xs;
        struct
        {
            Fs xs;
            F E;
        } place0;
        const Quality *quality = nullptr;
        Failure failure = Failure::Void;
        Result (const Target&);
        ~Result();
        int have = 0;
        double bits_for_print() const; // Precision.
        struct Rating
        {
            Failure failure = Failure::Void;
            struct { Value worst, best; } delta;
            F quality = 1000;
            F height = 0, height_p = 0, height_q = 1;
            double quality_p10 = -1000;
            double prec2 = std::nan ("");
            double prec10 = std::nan ("");
            void set (const Quality*, Result&);
            void print() const;
        } rating;
        F operator () (const F&, int* = nullptr) const;
        F f (const F&, const void*) const;
        Failure find_worst (int flags = 0);
        const Rating* rate (const Approx*);
        void adjust_worst_place0();
        F delta_height (const Result&, bool relative = false) const;
        void check_denominator() const;
    } result;

    static struct Verbose
    {
        int level;
        int E;
        int symm;
        int non0;
        int pq;
    } verbose;

    Approx (const Approx&) = delete;
    Approx (Approx&&) = default;
    Approx& operator = (const Approx&) = delete;
    Approx& operator = (Approx&&) = default;

    Approx (const Result& r)
        : flags(r.target.flags), target(r.target), xs_remez(r.xs),
          result(r.target)
    {
        result.E = r.E;
    }

    Approx (const Target& t, const Fs& xs, const F& e_old)
        : flags(t.flags), target(t), xs_remez(xs), result(t)
    {
        result.poly.q = Poly<F> { 1 };
        result.E = e_old;
    }

    void remez();
    const Result::Rating* rate();

    Values get_values (const Fs&) const;
    F delta (const F&, int* pdenom_sign = nullptr) const;
    double xacos (const F&) const;
    void write_deltas (FILE*, int, const char* = nullptr) const;
    void write_deltas (const char*, int, const char* = nullptr) const;
    void write_xacos (FILE*, int, bool) const;
    void write_xacos (const char*, int, bool) const;
    // Assume the target function is something like asinq2 over [0, 0.5],
    // and write deltas for asin and acos using this approximation.
    struct AsinContext
    {
        typedef Float (*Func)(const Float&, const AsinContext&);

        bool asin_p;
        int ulp;
        const Poly<Float> &p, &q;
        const Float &pi, &pi2;
        Func a, c, as, ac;
    };
    void write_deltas (FILE*, Func, const AsinContext&, AsinContext::Func, int, const char*) const;
    bool write_deltas_asin (const char*, int, bool, int ulp) const;
    void write_deltas_asin (FILE*, AsinContext&, int) const;

    static Poly<Float> best_polynomial (int, const Values&, const Values *fixed = nullptr);
    static Poly<double> best_polynomial (int, const DValues&, const DValues *fixed = nullptr);

    Poly<double> phase_poly (int, int n_points, bool minus_x) const;
};

#endif // MP_APPROX_H
