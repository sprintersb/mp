#ifndef RB_STRINGMAP_H
#define RB_STRINGMAP_H

#include "rb-tree.h"

#include <cassert>
#include <cstring>
#include <utility>      // std::move

/*  template<class T, bool copy_id>
    class RBStringMap
    is a mapping from ID strings to T's, aka RBStringMap::value_type.

    "nullptr" is a valid ID like any other string (const char*).
    The IDs impose a linear order as implied by  std::strcmp (id1, id2)
    if neither ID is nullptr, and otherwise by  id1 < id2.

    In the remainder, "element" stands for RBStringMap::element_type which
    describes one  ID -> T (=value_type) mapping.  RBStringMap is an RBTree
    of such mappings.

    RBStringMap inherits the invariance property from RBTree, namely
    neither element_type nor value_type are changed after an element_type
    node has been created:  They will neither be copied nor moved and will
    thus remain at the same memory location no matter which set operations
    will be performed.  The only exception is when an element is deleted,
    for example when  remove()  is called with a request to delete the found
    element.

    copy_id determines whether the IDs will be copies of the strings provided
    by, say, add (const char* id, T const&).  If copy_id is "true", then
    a copy of the ID will be made and managed.  If copy_id is "false", then
    no such copy is made, and the user takes responsibility that the lifetime
    of the IDs extend over the entire lifetime of the set's elements.

    Element Actions
    ===============

    const char* id() const;
        return the ID associated to the element.

    [const] T& value () [const];
        Return a (const) reference to the element's payload of value_type T.

    [const] Element* next() [const];
    [const] Element* prev() [const];
        Return the address of the next / previous element in the set, or
        nullptr if it is the last resp. first element.

    RBStringMap Actions
    ===================

    bool empty () const;
        Return true iff the set is empty.

    int size () const;
        return the number of elements in the set.

    Element* add (const char *id, T const&, bool *added_p = nullptr);
    Element* add (const char *id, T&&, bool *added_p = nullptr);
        Insert T into the set provided the ID does not already exist.
        Return the address of the newly created Element.  Otherwise, if ID
        already exists, return the respective Element's address.
        If  added_p  is non-null, then set it to true iff a new element has been
        added.

    [const] Element* find (const char *id) [const];
        return the address of the element with ID id, or nullptr if no such
        element exists.

    Element* remove (const char *id, bool delet = true);
        Remove Element with respective ID from the set provided it exists.
        If the ID exists and there is no request to delete the attached
        element, then return its address.  If the ID does not exists or
        if it is deleted as requested, return nullptr.

    Element* remove (Element *e, bool delet = true);
        Remove the respective element from the set.  No elements must have
        been added or removed since the element's address was retrieved.

    [const] Element* first() [const];
    [const] Element* last() [const];
        Return the address of the first resp. last element in the set,
        or nullptr if the set is empty.

    enum RBStringMap::Traverse { TopDown, DownTop, LeftRight, RightLeft };
        Used in the function below.

    int traverse (Traverse how, traverse_[const_]func f, void *data) [const];
        Traverse all elements of the set.  The set is internally represented
        as a binary tree, and HOW determines how that tree will be traversed.
        DATA defaults to nullptr and is passed to each invocation of F.
        If F returns 0, the traversal stops.  In any case, the return value
        of F will be returned by traverse. Traversing an empty set returns 0.

        Example:
        
            using S = RBStringMap<int,false>;
            S s;
            s.add (nullptr, 0);
            s.add ("two", 2);
            s.add ("six", 6);
            s.add ("zero", 0);
            s.traverse (S::Traverse::LeftRight,
                        [] (const S::Element &e, void*) -> int
                        {
                            printf ("'%s' -> %d\n", e.id(), e.value());
                            return 0; // Don't stop.
                        });
        Output:
            '(null)' -> 0
            'six' -> 6
            'two' -> 2
            'zero' -> 0

    Time Consumption
    ================

    All functions execute in time  O(log N)  with N = size(), i.e. the number
    of elements currently held in the set.

    Exceptions to the rule:
    *   size() and empty() perform in constant time.
    *   Deletion takes time O(n), same for traversing the whole set.
*/

namespace
{
    template<bool b> struct _rbStringMap_IdType;
    template<> struct _rbStringMap_IdType<true> { typedef char* type; };
    template<> struct _rbStringMap_IdType<false> { typedef const char* type; };
} // ::anon


template<class T, bool copy_id>
class RBStringMap
{
public:
    struct Element;
    typedef RBTree<Element> Tree;
    typedef typename Tree::Node Node;
    typedef T value_type;
    typedef Element element_type;
    typedef int (*traverse_func) (Element&, void*);
    typedef int (*traverse_const_func) (const Element&, void*);

    typedef typename _rbStringMap_IdType<copy_id>::type IdType;
    static constexpr bool copy_id_ = copy_id;

    struct Element
    {
        friend RBStringMap;
        static constexpr bool copy_id_ = copy_id;

    private:
        IdType id_;
        Node* node_ = nullptr;
        T t_;

        Element() = delete;

        static char* cid (const char *id)
        {
            char *xid = const_cast<char*> (id);
            if (id && copy_id_)
            {
                xid = new char[1 + std::strlen (id)];
                std::strcpy (xid, id);
            }
            return xid;
        }

        Element (const Element& e) : Element (e.id_, e.node_, e.t_) {}

        Element& operator = (const Element&) = delete;
        Element& operator = (Element&&) = delete;

        Element (const char* id, Node *node, T const& t)
            : id_(cid(id)), node_(node), t_(t) {}

        Element (const char* id, Node *node, T&& t)
            : id_(cid(id)), node_(node), t_(std::move (t)) {}

    public:
        Element (Element&& e)
        {
            std::memcpy (this, &e, sizeof (Element));
            e.id_ = nullptr;
        }

        ~Element()
        {
            if (copy_id_)
                delete [] id_;
        }

        const char* id () const { return id_; };

        // non-const
        T& value() { return t_; }
        Element* next() { return node_->next(); }
        Element* prev() { return node_->prev(); }

        // Dito, as const
        const T& value() const { return t_; }
        const Element* next() const { return node_->next(); }
        const Element* prev() const { return node_->prev(); }
       
        bool operator < (Element const& x) const
        {
            return idcmp (id_, x.id_) < 0;
        }

    private:
        struct Search
        {
            int operator () (Element const& e, const char* const& id) const
            {
                return idcmp (e.id_, id);
            }
        };

        static int idcmp (const char *a, const char *b)
        {
            return 0?0
                : a && b ? std::strcmp (a, b)
                : a < b  ? -1
                : a > b;
        }
    }; // struct Element

private:
    Tree tree_;

public:
    Element* add (const char *id, T const& t, bool *added_p = nullptr)
    {
        int size = tree_.size();
        Node *n = tree_.add (std::move (Element (id, nullptr, t)),
                             true /* unique */);
        assert (n);

        if (added_p)
            *added_p = size != tree_.size();

        if (size != tree_.size())
            n->t_.node_ = n;

        return & n->t_;
    }

    Element* add (const char *id, T&& t, bool *added_p = nullptr)
    {
        int size = tree_.size();
        Element e = Element (id, nullptr, std::move (t));
        Node *n = tree_.add (std::move (e), true /* unique */);
        assert (n);

        if (added_p)
            *added_p = size != tree_.size();

        if (size != tree_.size())
            n->t_.node_ = n;

        return & n->t_;
    }

    const Element* find (const char *id) const
    {
        auto *n = tree_.template search<typename Element::Search, const char*>(id);
        return n ? & n->t_ : nullptr;
    }

    Element* find (const char *id)
    {
        auto *n = tree_.template search<typename Element::Search, const char*>(id);
        return n ? & n->t_ : nullptr;
    }

    Element* remove (const char *id, bool delet = true)
    {
        Node *n = tree_.template
            search_and_remove<typename Element::Search, const char*> (id, delet);

        return n ? nullptr : & n->t_;
    }

    Element* remove (Element *e, bool delet = true)
    {
        Node *n = e ? tree_.remove (e->node_, delet) : nullptr;
        return n ? & n->t_ : nullptr;
    }

    int size () const { return tree_.size(); }
    bool empty () const { return tree_.empty(); }

    Element* first() { Node *n = tree_.first(); return n ? & n->t_ : nullptr; }
    Element* last() { Node *n = tree_.last(); return n ? & n->t_ : nullptr; }
    const Element* first() const { const Node *n = tree_.first(); return n ? & n->t_ : nullptr; }
    const Element* last() const { const Node *n = tree_.last(); return n ? & n->t_ : nullptr; }

    typedef typename Tree::Traverse Traverse;

private:
    // Traverse helpers.
    
    struct TraverseData
    {
        traverse_const_func f_const;
        traverse_func f;
        void *data;
    };

    static int f_traverse_constTree (const typename Tree::Node& n, void *v)
    {
        auto td = (TraverseData*) v;
        td->f_const (n.t_, td->data);
        return 0;
    };
    static int f_traverse_Tree (typename Tree::Node& n, void *v)
    {
        auto td = (TraverseData*) v;
        td->f (n.t_, td->data);
        return 0;
    };

public:
    int traverse (Traverse how, traverse_const_func f, void *data = nullptr) const
    {
        TraverseData td;
        td.f_const = f;
        td.data = data;
        return tree_.traverse (how, f_traverse_constTree, & td);
    }
    int traverse (Traverse how, traverse_func f, void *data = nullptr)
    {
        TraverseData td;
        td.f = f;
        td.data = data;
        return tree_.traverse (how, f_traverse_Tree, & td);
    }
};

#endif // RB_STRINGMAP_H
