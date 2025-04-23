#include "mp-cheby.h"

#include "diagnostic.h"
#include "rb-tree.h"

#include <cassert>
#include <utility>

namespace
{
    struct Comparator
    {
        bool operator () (const Poly<Ratio>& a, const Poly<Ratio>& b) const
        {
            return a.deg() < b.deg();
        }

        int operator () (const Poly<Ratio>& a, const int& b) const
        {
            return b - a.deg();
        }
    };

    using TnTree = RBTree<Poly<Ratio>, Comparator>;

    TnTree _Tn_tree;

    const TnTree::Node* _Tn_node (int n)
    {
        assert (n >= 0);
        const TnTree::Node *node = _Tn_tree.search<Comparator,int> (n);
        if (!node)
        {
            Poly<Ratio> p;

            if (n == 0) p = Poly<Ratio> { 1_Ratio };
            if (n == 1) p = Poly<Ratio> { 0_Ratio, 1_Ratio };
            if (n > 1)
            {
                const TnTree::Node *n2 = _Tn_node (n - 2);
                const TnTree::Node *n1 = _Tn_node (n - 1);
                assert (n1 && n2);
                // T_n(x) = 2 * x * T_{n-1}(x) - T_{n-2}(x).
                p = 2_Ratio * (n1->t_ << 1) - n2->t_;
            }

            node = _Tn_tree.add (std::move(p), true);
            assert (node);
        }
        return node;
    }

    FVector<Float> _to_vec (const Poly<Float>& p, int dim)
    {
        FVector<Float> v { Alloc (dim) };
        for (int i = 0; i < dim; ++i)
            v[i] = i > p.deg() ? 0_Float : p[i];
        return v;
    }

    Poly<Float> _to_poly (const FVector<Float>& v)
    {
        Poly<Float> p;
        for (int i = v.dim() - 1; i >= 0; --i)
            p.set (i, v[i]);
        return p;
    }

    template<class F> F _cos (const F& x) { return std::cos (x); }
    template<class F> F _pi() { return 3.1415926535897932384626433832795; }

    template<> Float _cos (const Float& x) { return x.cos(); }
    template<> Float _pi() { return Float::pi(); }
}; // ::anon


namespace Cheby
{
    const Poly<Ratio>& Tn (int n)
    {
        assert (n >= 0);
        return _Tn_node (n)->t_;
    }

    void Tn_clear()
    {
        _Tn_tree.remove_all();
    }

    FMatrix<Float> Tn_matrix (int n)
    {
        int dim  = 1 + n;
        FMatrix<Float> m { FMatrix<Float>::make (dim, dim, 0_Float) };

        for (int c = 0; c < dim; ++c)
        {
            const Poly<Ratio>& tn = Cheby::Tn (c);

            for (int i = 0; i <= c; ++i)
                if (! tn[i].is_zero())
                    m[i][c] = (Float) tn[i];
        }
        return m;
    }

    FVector<Float> basis_Tn_to_xn (const FVector<Float>& v)
    {
        auto m = Cheby::Tn_matrix (v.dim() - 1);
        auto w = m * v;
        return w;
    }

    FVector<Float> basis_xn_to_Tn (const FVector<Float>& v)
    {
        auto m = Cheby::Tn_matrix (v.dim() - 1);
        auto w = m.solve_triangle (v);
        return w;
    }

    Poly<Float> basis_Tn_to_xn (const Poly<Float>& p, int dim)
    {
        if (dim != -1 && dim <= p.deg())
            error ("dim %d must be at least %d for a polynomial of degree %d",
                   dim, 1 + p.deg(), p.deg());

        auto v = _to_vec (p, dim == -1 ? 1 + p.deg() : dim);
        auto w = Cheby::basis_Tn_to_xn (v);
        return _to_poly (w);
    }

    Poly<Float> basis_xn_to_Tn (const Poly<Float>& p, int dim)
    {
        if (dim != -1 && dim <= p.deg())
            error ("dim %d must be at least %d for a polynomial of degree %d",
                   dim, 1 + p.deg(), p.deg());

        auto v = _to_vec (p, dim == -1 ? 1 + p.deg() : dim);
        auto w = Cheby::basis_xn_to_Tn (v);
        return _to_poly (w);
    }

    template<class F>
    F Tn_zero (int n, int k)
    {
        assert (k >= 0 && k < n);
        if (n == 2 * k + 1)
            return F { 0.0 };
        F kk { (double) (2 * (n - k) - 1) };
        F nn { (double) (2 * n) };
        return _cos<F> (kk / nn * _pi<F>());
    }

    template<class F>
    F Tn_maximum (int n, int k)
    {
        assert (k >= 0 && k <= n);
        if (k == 0)
            return F { -1.0 };
        if (k == n)
            return F { 1.0 };
        if (2 * k == n)
            return F { 0.0 };
        F kk { (double) (n - k) };
        F nn { (double) n };
        return _cos<F> (kk / nn * _pi<F>());
    }

    template double Tn_zero (int, int);
    template Float  Tn_zero (int, int);

    template double Tn_maximum (int, int);
    template Float  Tn_maximum (int, int);

}; // ::Cheby
