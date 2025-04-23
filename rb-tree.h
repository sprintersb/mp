#ifndef RB_TREE_H
#define RB_TREE_H

#include <cassert>
#include <cstring>
#include <iostream>
#include <utility>      // std::move
#include <functional>   // std::less

/*  Red-Black Tree implementation.

    In the remainder, "Node" is short for RBTree<T>::Node.


    Complexity / Execution Times
    ============================

    If N is the current size (in number of nodes) in the tree,
    then the following execution times are consumed:

    Time O(1):
    *   int size() const;
    *   bool empty() const;

    Time O(log n):

    In RBTree<T>:
    *   Node* add (const T&, bool unique = false, bool *pchanged = nullptr);
    *   Node* add (T&&, bool unique = false, bool *pchanged = nullptr);
    *   Node* add (Node*);
    *   Node* remove (Node*, bool delet = true);
    *   [const] Node* find (const T&) [const];
    *   template<class S, class V>
        [const] Node* search (const V& val) [const]
    *   template<class S, class V>
        Node* search_and_remove (const V& val, bool delet = true);
    o   [const] Node* extremum (RBTree<T>::Side) [const];
    o   [const] Node* first()  [const];
    o   [const] Node* last() [const];

    In RBTree<T>::Node:
    o   [const] T& value() [const];
    o   [const] Node* next() [const];
    o   [const] Node* prev() [const];
    o   [const] Node* neighbour (RBTree<T>::Side) [const];

    Time O(n):
    *   remove_all();


    Tree / Node Operations and their Effects on Nodes' Payloads
    ===========================================================

    Methods listed with  o  are pure tree operations, they do not need the
    payload value of any node, e.g. for comparision.

    The payload values t_ held in the existing nodes are NEITHER COPIED NOR
    MOVED NOR CHANGED in any way.  This means you can have some Node* NODE,
    and after whatever insertion or deletion you are performing on the tree,
    NODE's t_ field will be unaltered.  The only exception is when NODE itself
    is removed from the tree and deleted as requested by  remove (NODE, true).


    Finding a Node
    ==============

    In order to find nodes / values in the tree, there are several ways:

    * If a _Compare type was specified as 2nd template argument, then it's
         bool operator () (const T& a, const T& b) const;
      operator is invoked.  The operator shall return true iff  a < b.

    * If no _Compare type was specified as 2nd argument, then T's
         bool operator < (const T&) const;
      operator is used to determine how two T's compare.

    * A node can be searched with
         template<class S, class V> [const] Node* search (const V& val) [const];
      S is a type providing
         int operator () (const T&, const V&) const; // Return < 0, 0, > 0.
      that will be used to <=> compare nodes' payloads against VAL.  search()
      returns the first node encountered for which 0 (equal) is returned,
      or nullptr if no such node is found.  In order to work as expected,
      this method must be implemented in such a way that it is compatible with
      the usual comparison as used by  add()  and  find().  Suppose for example
      a tree that stores polynomials of different degrees.  In order to find
      such a polynomial by it's degree,  search()  can be used.

    * The complete tree could be traversed, which is quite slow of course.
      Example:

        using T = RBTree<X>;
        T::traverse_[const_]func f =
            [] ([const] T::Node& n, void *v) -> int
            {
                if (n.t_ == ...)
                {
                    * ([const] T::Node**) (v) = &n;
                    return 1; // Stop traversal.
                }
                return 0; // Don't stop yet.
            };
        [const] T::Node *node = nullptr;
        int stopped = t.traverse (T::TopDown, f, &node);
*/

template<class T, typename _Compare = std::less<T>>
class RBTree
{
public:
    class Node;
    enum Color { Red, Black };
    enum Traverse { TopDown, DownTop, LeftRight, RightLeft };
    enum Side { Left, Right };
    typedef T value_type;
    typedef Node node_type;
    typedef _Compare value_compare;
    typedef int (*traverse_func) (Node&, void*);
    typedef int (*traverse_const_func) (const Node&, void*);

    class Node
    {
        friend RBTree;
    private:
        Color color_ = Red;
    public:
        T t_;
        Node *dad = nullptr;
        Node *son[2] = { nullptr, nullptr };

        Node (const T& t) : t_(t) {}
        Node (T&& t) : t_(std::move(t)) {}

        T& value () { return t_; }
        const T& value () const { return t_; }

        bool is (Side s) const { return this == dad->son[s]; }
        Side side() const { return (Side) (this == dad->son[Right]); }
        int n_sons() const { return (!!son[0]) + (!!son[1]); }

        Color color() const { return color_; }

private:
        static Color color (const Node *n)
        {
            return n ? n->color_ : Black;
        }

        const Node* opa() const
        {
            assert (dad);
            assert (dad->dad);
            return dad->dad;
        }

        const Node* bruder() const
        {
            assert (dad);
            return dad->son[!side()];
        }

        const Node* onkel() const
        {
            assert (dad);
            assert (dad->dad);
            return dad->bruder();
        }

#define NON_CONST(T, METH) \
        T* METH() { return const_cast<T*> (((const T*) this)->METH()); }

#define NON_CONST_s(T, METH) \
        T* METH(Side s) { return const_cast<T*> (((const T*) this)->METH(s)); }

        NON_CONST (Node, opa)
        NON_CONST (Node, bruder)
        NON_CONST (Node, onkel)

        int traverse (const int *id, traverse_const_func f, void *data) const
        {
            int stop = 0;
            for (int i = 0; i < 3 && !stop; ++i)
                if (id[i] < 0)
                    stop = f (*this, data);
                else if (son[id[i]])
                    stop = son[id[i]]->traverse (id, f, data);
            return stop;
        }

        int traverse (Traverse how, traverse_const_func f, void *data) const
        {
            static const int ids[][3] =
            {
                { -1, 0, 1 }, // TopDown
                { 0, 1, -1 }, // DownTop
                { 0, -1, 1 }, // LeftRight
                { 1, -1, 0 }, // RightLeft
            };
            return traverse (ids[how], f, data);
        }

        void insert (RBTree *t)
        {
            if (!dad)
            {
                color_ = Black;
                return;
            }
            else if (dad->color_ == Black)
                return;

            if (color (onkel()) == Red)
            {
                dad->color_ = Black;
                onkel()->color_ = Black;
                opa()->color_ = Red;
                opa()->insert (t);
                return;
            }

            Node *next = this;
            Side s = (Side) !side();
            if (dad->is(s))
            {
                dad->rotate (t, s);
                next = son[s];
            }
            next->insert_finish (t);
        }

        void insert_finish (RBTree *t)
        {
            dad->color_ = Black;
            opa()->color_ = Red;
            bool ll = is(Left) && dad->is(Left);
            assert (ll || (is(Right) && dad->is(Right)));
            opa()->rotate (t, (Side) ll);
        }

        void rotate (RBTree *t, Side s)
        {
            Node *n = son[!s];

            replace (t, n);

            son[!s] = n->son[s];
            if (son[!s])
                son[!s]->dad = this;
            n->son[s] = this;
            dad = n;
        }

        void replace (RBTree *t, Node *neu)
        {
            if (! dad)
                t->root_ = neu;
            else
                dad->son[side()] = neu;
            if (neu)
                neu->dad = dad;
        }

        std::ostream& dump (std::ostream& ost, const char *sep) const
        {
            auto& l = son[0];
            auto& r = son[1];

            if (dad && !l && !r)
                return ost;

            ost << sep << t_ << "{";
            if (l) ost << l->t_;
            ost << ";";
            if (r) ost << r->t_;
            ost << "}";

            if (l) l->dump (ost, sep);
            if (r) r->dump (ost, sep);
            return ost;
        }

        const Node* extremum (Side s) const
        {
            for (const Node *n = this; ; n = n->son[s])
                if (! n->son[s])
                    return n;
            assert (0);
        }
        NON_CONST_s (Node, extremum)

        void remove_start (RBTree *t)
        {
            assert (n_sons() <= 1);
            Node *kid = son[!! son[1]];
            if (color_ == Black)
            {
                color_ = color (kid);
                remove (t);
            }
            replace (t, kid);
            if (kid && ! dad)
                // Root is Black.
                kid->color_ = Black;
        }

        void remove (RBTree *t)
        {
            if (! dad)
                return;

            if (color (bruder()) == Red)
            {
                dad->color_ = Red;
                bruder()->color_ = Black;
                dad->rotate (t, side());
            }

            Node* bro = bruder();
            if (color (bro) == Black
                && color (bro->son[0]) == Black
                && color (bro->son[1]) == Black)
            {
                bro->color_ = Red;
                if (color (dad) == Black)
                    dad->remove (t);
                else
                    dad->color_ = Black;
                return;
            }

            Side s = side();
            if (color (bro) == Black
                && color (bro->son[s]) == Red
                && color (bro->son[!s]) == Black)
            {
                bro->color_ = Red;
                bro->son[s]->color_ = Black;
                bro->rotate (t, (Side) !s);
            }

            bro = bruder();
            bro->color_ = dad->color_;
            dad->color_ = Black;

            s = side();
            assert (color (bro->son[!s]) == Red);
            bro->son[!s]->color_ = Black;
            dad->rotate (t, s);
        }
public:
        const Node* neighbour (Side s) const
        {
            if (son[s])
                return son[s]->extremum ((Side) !s);

            for (const Node *d = this; d && d->dad; d = d->dad)
                if (d->is((Side) !s))
                    return d->dad;

            return nullptr;
        }
        NON_CONST_s (Node, neighbour)

        const Node* next() const { return neighbour (Right); }
        const Node* prev() const { return neighbour (Left); }
        NON_CONST (Node, prev)
        NON_CONST (Node, next)
    }; // Node

    Node *root_ = nullptr;
    int size_ = 0;

    RBTree() {}
    ~RBTree() { kill(); }

    RBTree (RBTree&& t)
    {
        std::memcpy (this, &t, sizeof (*this));
        t.root_ = nullptr;
        t.size_ = 0;
    }

    RBTree& operator = (RBTree&& t)
    {
        kill();
        std::memcpy (this, &t, sizeof (*this));
        t.root_ = nullptr;
        t.size_ = 0;
        return *this;
    }

    // Todo.
    RBTree (const RBTree&) = delete;
    RBTree& operator = (const RBTree&) = delete;

    bool empty() const { return !root_; }
    int size() const { return size_; }

    int traverse (Traverse how, traverse_const_func f, void *data = nullptr) const
    {
        return root_ ? root_->traverse (how, f, data) : 0;
    }
    int traverse (Traverse how, traverse_func f, void *data = nullptr)
    {
#if 0
        auto cf = (traverse_const_func) f;
#else
        traverse_const_func cf;
        __asm ("" : "=r" (cf) : "0" (f));
#endif
        return root_ ? ((const Node*) root_)->traverse (how, cf, data) : 0;
    }

private:
    void kill ()
    {
        traverse (DownTop, [] (Node& n, void*) -> int
        {
            delete &n;
            return 0;
        });
    }

    Node** find_insert (const T& t, bool unique, Node **dad, bool& exists)
    {
        const value_compare& less = value_compare();

        *dad = nullptr;
        Node **son = &root_;
        exists = false;

        for (Node *node = *son; node; node = *son)
        {
            bool is_less = less (t, node->t_);
            if ( (exists = unique && ! is_less && ! less (node->t_, t)) )
                break;
            son = & node->son[is_less ? Left : Right];
            if (! *son)
            {
                *dad = node;
                break;
            }
        }

        return son;
    }

    Node* do_insert (Node *neu, Node *dad)
    {
        neu->dad = dad;
        neu->insert (this);
        size_++;
        return neu;
    }

public:

    Node* add (const T& t, bool unique = false, bool *pchanged = nullptr)
    {
        bool exists;
        Node *dad, **son = find_insert (t, unique, &dad, exists);

        if (pchanged)
            *pchanged = !exists || !unique;

        return exists
            ? *son
            : do_insert ((*son = new Node (t)), dad);
    }

    Node* add (T&& t, bool unique = false, bool *pchanged = nullptr)
    {
        bool exists;
        Node *dad, **son = find_insert (t, unique, &dad, exists);

        if (pchanged)
            *pchanged = !exists || !unique;

        return exists
            ? *son
            : do_insert ((*son = new Node (std::move (t))), dad);
    }

    Node* add (Node *neu)
    {
        Node *dad, **son = find_insert (neu->t_, false, &dad);
        assert (son);

        neu->son[0] = neu->son[1] = nullptr;
        neu->color_ = Red;
        return do_insert ((*son = neu), dad);
    }

    const Node* find (const T& t) const
    {
        const value_compare& compa = value_compare();
        for (const Node *n = root_; n; )
        {
            bool less = compa (t, n->t_);
            if (!less && !compa (n->t_, t))
                return n;
            n = n->son[!less];
        }
        return nullptr;
    }

    Node* find (const T& t)
    {
        return const_cast<Node*> (((const RBTree*) this)->find (t));
    }

    template<class S, class V>
    const Node* search (V const& val) const
    {
        const S& s = S();
        for (Node *n = root_; n; )
        {
            int star = s (n->t_, val);
            if (star == 0)
                return n;
            n = n->son[star < 0];
        }

        return nullptr;
    }
    template<class S, class V>
    Node* search (const V& val)
    {
        return const_cast<Node*> (((const RBTree*) this)->search<S,V> (val));
    }

    const Node* extremum (Side s) const
    {
        return root_ ? root_->extremum (s) : nullptr;
    }
    Node* extremum (Side s)
    {
        return root_ ? root_->extremum (s) : nullptr;
    }
    const Node* first()  const { return extremum (Left); }
    const Node* last() const { return extremum (Right); }
    Node* first()  { return extremum (Left); }
    Node* last() { return extremum (Right); }

    void remove_all()
    {
        kill();
        root_ = nullptr;
        size_ = 0;
    }

    Node* remove (Node *n, bool delet = true)
    {
        if (!n)
            return nullptr;

        if (n->n_sons() == 2)
        {
            // The left neighbour has 1 son at most.
            Node *left = n->prev();
#if 0
            // Move from left neighbour and then delete that instead.
            // FIXME:  We could re-wire the nodes to avoid moving the payload.
            n->t_ = std::move (left->t_);
            n = left;
#else
            // Rewire such that n and it's left neighbour switch places.
            // This avoids moving or copying the payload, however exchanging
            // is more contrived than "simply" copying the payload like above.
            exchange (n, left);
#endif
        }

        n->remove_start (this);
        if (--size_ == 0)
            root_ = nullptr;

        if (delet)
        {
            delete n;
            n = nullptr;
        }
        else
            n->son[0] = n->son[1] = n->dad = nullptr;

        return n;
    }

    template<class S, class V>
    Node* search_and_remove (const V& val, bool delet = true)
    {
        Node *n = search<S,V> (val);
        return remove (n, delet);
    }


private:
    // Rewire A and B in such a way that A and B are effectively swapped,
    // but without touching their contents t_ in any way.
    void exchange (Node *a, Node *b)
    {
        // This follows from the caller's purpose.
        assert (a != b);
        assert (a->n_sons() == 2);
        assert (b->n_sons() < 2 && b->dad);

        Node* *dad_a_son = a->dad ? & a->dad->son[a->side()] : & root_;
        Node* *dad_b_son = & b->dad->son[b->side()];
        assert (dad_a_son && a == *dad_a_son);
        assert (dad_b_son && b == *dad_b_son);

        // Fix nodes in the vincinity of A and B, including root_ as needed.

        bool b_sonof_a = b->dad == a;

        if (b_sonof_a)
            assert (b->side() == Left);
        else
        {
            a->son[0]->dad = b;
            *dad_b_son = a;
        }

        *dad_a_son = b;
        a->son[1]->dad = b;
        if (b->son[0]) b->son[0]->dad = a;
        if (b->son[1]) b->son[1]->dad = a;

        // Fix nodes A and B themselves.

        std::swap (a->color_, b->color_);
        std::swap (a->son, b->son);
        std::swap (a->dad, b->dad);

        if (b_sonof_a)
        {
            assert (a == a->dad);
            assert (b == b->son[0]);
            b->son[0] = a;
            a->dad = b;
        }
    }

public:
    std::ostream& dump (std::ostream& ost, const char *sep = ", ") const
    {
        if (root_)
            root_->dump (ost, sep);
        else
            ost << "(null)";
        return ost;
    }
}; // RBTree

namespace
{
    template<class T>
    std::ostream& operator << (std::ostream& ost, const RBTree<T>& tree)
    {
        return tree.dump (ost);
    }

    template<class T, typename C>
    std::ostream& operator << (std::ostream& ost, const RBTree<T,C>& tree)
    {
        return tree.dump (ost);
    }
};

#endif // RB_TREE_H
