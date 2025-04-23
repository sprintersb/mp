#include "expr-parse.h"

#include "rb-tree.h"

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <cmath>
#include <cstdarg>

#include "diagnostic.h"

using EP = ExpressionParser;

#define NORETURN __attribute__((__noreturn__))
#define DEPRECATED __attribute__((__deprecated__))

// Payload in a RBTree node that maps a string (identifier name) to
// a value or a function pointer.  If the value / function for an
// identifier is looked up, the found Ident's address will be stored
// in the Token.ident_ field for swift future lookup.
// This uses the feature that RBTree::Node's never change their payload.
// If a function requires exactly 1 argument, then it's pointer is f_.
// If a function requires more than 1 argument or accepts a variable
// number of arguments like "min" and "max", then it's a Token method
// pointer.

struct EP::Ident
{
    typedef double (*Func)(double);
    typedef double (Token::* FuncN)() const;
    typedef Token* (*Maker)(Lexer*, char*, const char*, const char*, Ident*);

    const char *name_;
    FuncN fn_ = nullptr;
    Func fx_ = nullptr;
    Maker make_ = nullptr;

    double value_ = std::nan("");
    char n_args = -1;
    bool built_in_ = true;
    bool undefined_ = false;

    Ident (const char *str, double val, bool bi)
        : name_(str), value_(val), n_args(0), built_in_(bi) {}

    Ident (int n, const char *name, Maker m) : name_(name), make_(m), n_args(n) {}

    // How RBTree::Node's compare.
    bool operator < (const Ident& x) const
    {
        return std::strcmp (name_, x.name_) < 0;
    }

    // To find RBTree::Node with a given name.
    struct Search
    {
        int operator () (const Ident& i, const char* const& name) const
        {
            return std::strcmp (i.name_, name);
        }
    };
};

using Tree = RBTree<EP::Ident>;

class EP::Lexer
{
    friend EP;
    friend Token;
private:

    const char* const text_;
    const char* caret_;

    int flags_, id_ = 0;

    Token *sp = nullptr;
    bool no_hexfloat_;

public:
    struct
    {
        bool tokenize, parse, lookup, alloc, tree;
    } debug;

private:
    // Set of known identifiers.  This includes known constants like "pi",
    // functions like "sin" and (helper) variables like "a" that can be set
    // and used by "a = x^2, sin (pi * a)" which will be evaluated the same
    // like "sin (pi * a^2)".
    Tree ident_;

    // The Ident* for "x" to speed up as expressions usually
    // have "x" as input variable
    Ident *ident_x = nullptr;

    auto search (const char* name) const -> const Tree::Node*
    {
        return ident_.search<Ident::Search, const char*> (name);
    }
    auto search (const char* name) -> Tree::Node*
    {
        return ident_.search<Ident::Search, const char*> (name);
    }

    Tree::Node* get_node (const char *name, double val, bool is_const)
    {
        auto node = ident_.add (Ident (name, val, is_const), true /* unique */);
        assert (node);
        return node;
    }

    Ident* get_ident (const char *name, double val, bool is_const)
    {
        bool is_x = name[0] == 'x' && name[1] == '\0';
        if (is_x && ident_x)
            return ident_x;

        Tree::Node *node = get_node (name, val, is_const);

        Ident *ident = & node->t_;
        if (is_x)
            ident_x = ident;
        return ident;
    }

    Ident* set (const char *name, double);
    Ident* set (const char *name, Ident::Func);
    Ident* set (int n, const char *name, Ident::FuncN);
    Ident* set_var (const char *name, double);

    void set (const char *a, const char *b, Ident::Func f)
    {
        set (a, f);
        set (b, f);
    }

public:
    class Number;
    class Identifier;
    class Variable;
    class Operator;
    class Function;
    struct Location { int from, to; };

    // A singly linked list of tokens due to input string tokenization.
    Token *tokens_ = nullptr, *prev_ = nullptr;

public:
    Lexer (const char *text, int flags);
    ~Lexer ();
    void init_ident();
    NORETURN static void leave_va (const char*, int, int, const char*, va_list);
    __attribute__((__format__(printf,4,5)))
    NORETURN void leave (int, int, const char*, ...) const;
    __attribute__((__format__(printf,3,4)))
    NORETURN void leave (Location, const char*, ...) const;

    void add (Token*);
    Token* lex ();
    Token* tokenize();
    Token* parse ();
    bool expect (const char*, const Token* = nullptr) const;
    static int prio (const Token*);

    Location loc (const char *start, const char *end) const
    {
        return Location { (int)(start - text_), (int)(end - text_) };
    }
    // Stack of tokens used during parse.  If all goes well, after parsing
    // SP will contain exactly one token which represents the syntax tree.
    void push (Token*, const char* = "");
    Token* pop ();
    void diagnostic_after_parse () const;
    void diagnostic_after_push (const Token*) const;
    NORETURN void leave_missing_open_paren (const Token*) const;
    void dump_sp (const Token*) const;
};


class EP::Token
{
    friend EP;
    friend Lexer;
protected:
    Lexer* const lexer_;
    // Copy of a portion of the original string.
    char *text_ = nullptr;
    // Pointer to the 1st char and one past the last char in the original
    // string lexer_->text_.
    const char *const start_;
    const char *end_;

    // A unique number 0 ... lexer_->id_ - 1.
    const int id_;

    // For quick identify: text_[0] for operator, 'N' = number, 'X' = ident.
    char c_;

    // For operators...
    int prio = 0;

    // This is a n-ary node, n = 0...2.
    int n_ary = 0;

    // Sons, n_arg of them.
    Token *arg[2] = { nullptr, nullptr };

    // The tokens are held in a linked list.
    Token *prev_, *next_ = nullptr;

    // Only used temporarily during parse to build a stack (linked list).
    Token *sp = nullptr;

    void alloc_text (const char *end);

    // Multivariate function implementors...
    double min () const
    {
        return c_ == ',' ? std::fmin (arg[0]->min(), arg[1]->min()) : value();
    }
    double max () const
    {
        return c_ == ',' ? std::fmax (arg[0]->max(), arg[1]->max()) : value();
    }
    double pow () const;
    double sat () const;
    double poly () const;
    double atan2 () const { return std::atan2 (arg[0]->value(), arg[1]->value()); }
    double hypot () const { return std::sqrt (hypot2()); }

    // ...and some helpers.
    double hypot2 () const;
    double horner (double, double) const;
    const Token* first_arg () const { return c_ == ',' ? arg[0]->first_arg() : this; }
    const Token* leftmost () const { return n_ary && ! is_pre ? arg[0]->leftmost() : this; }
    const Token* rightmost () const { return n_ary && ! is_post ? arg[n_ary-1]->rightmost() : this; }

    const Token* nth_value (int, int) const;
    int n_values () const
    {
        return c_ == ',' ? arg[0]->n_values() + arg[1]->n_values() : 1;
    }

public:
    mutable Ident *ident_ = nullptr;

    NORETURN void leave (const char*, ...) const;
    virtual double value () const = 0;

    double value (double x) const
    {
        if (lexer_->ident_x)
        {
            lexer_->ident_x->value_ = x;
            lexer_->ident_x->undefined_ = false;
        }
        else
            lexer_->set_var ("x", x);

        return value ();
    }
protected:
    union
    {
        struct
        {
            // The 3 main categories.
            bool is_id  : 1;
            bool is_num : 1;
            bool is_op  : 1;

            bool is_var : 1; // is_id with n_ary == 0.
            bool is_func: 1; // is_id with n_ary != 0.
            bool is_end : 1;
            bool is_pre : 1;
            bool is_post: 1;
            bool is_infix: 1;
            bool is_const: 1; // built-in.
            bool is_value: 1; // Ready to be evaluated by means of value().
            bool is_list: 1;  // 'L , v' for list of expressions.
        };
        int is_;
    };

protected:
    Token (Lexer*, char *text, const char *start, const char *end);
    Token (Lexer *l, const char *start) : Token (l, nullptr, start, start) {}
    virtual ~Token ();
    void finish ();

    virtual bool check (int) const { return true; }
    void check_n_values () const;

    void print (bool paren = false) const;
    void dump () const;
    void dump_id (const char* = "") const;
    const char* syntax () const;
};

/*******************************************************************/
/*******************************************************************/
/*******************************************************************/

double EP::eval (const char *text, double x, int flags)
{
    Expression *ex = parse (text, flags);
    double val = ex->value (x);
    delete ex;
    return val;
}

bool EP::supports_hexfloat ()
{
    char *tail;
    (void) strtod ("0x1.2p0", &tail);
    return *tail == '\0';
}

// EP::Expression is just a thin wrapper around Token in order to keep
// the header clean from all the gory bells and whistles.

EP::Expression::Expression (Token *token)
    : expression_(token)
{}

EP::Expression::~Expression ()
{
    // FIXME: We MUST NOT "delete expression_->lexer_;" because that may crash
    // as the compiler may use "expression_" AFTER the respective Token has
    // been deleted by the Lexer.  This works around a GCC bug, see
    // https://gcc.gnu.org/PR52339
    // https://stackoverflow.com/a/76172573/1556746
    Lexer *l = expression_->lexer_;
    delete l;
}

double EP::Expression::value () const
{
    return expression_->value();
}

double EP::Expression::value (double x) const
{
    return expression_->value (x);
}

const char* EP::Expression::text () const
{
    return expression_->lexer_->text_;
}

void EP::Expression::set (const char *name, double val) const
{
    expression_->lexer_->set_var (name, val);
}

// Register a new token.
void EP::Lexer::add (Token *token)
{
    // Manage their administration in a light-weight list.

    if (! tokens_)
        tokens_ = token;

    if (prev_)
        prev_->next_ = token;

    token->prev_ = prev_;
    prev_ = token;
}

EP::Token::Token (Lexer *lexer, char *text, const char* start, const char *end)
    : lexer_(lexer), text_(text), start_(start), end_(end), id_(lexer->id_++)
{
    is_ = 0;
    lexer_->add (this);
}


// End points to intended position of the string end '\0'.
void EP::Token::alloc_text (const char *end)
{
    auto len = end - start_;
    text_ = new char[1 + len];
    std::strncpy (text_, start_, len);
    text_[len] = '\0';
    end_ = end;
}

class EP::Lexer::Number : public EP::Token
{
    int radix_ = 10;
    // Pointers into the original string or 0.
    const char *dot_ = nullptr;
    const char *ep_ = nullptr;
    double value_;

public:
    Number (EP::Lexer *lexer, const char*& str)
        : Token (lexer, str)
    {
        //info (lexer_->debug.tokenize, "...%s", str);
        is_num = is_value = true;
        n_ary = 0;
        c_ = 'N';
        while (consumes (str))
            ++str;

        alloc_text (str);
assert (text_ && *text_);
        value_ = strtod (text_, nullptr);
        check();
    }

    ~Number()
    {
        info (lexer_->debug.alloc, "~Number(%p) \"%s\"", this, text_);
    }

    double value() const override
    {
        return value_;
    }

    bool consumes (const char *str)
    {
        char c = *str;
        bool ok = true;

        if (c == '\0')
            ok = false;
        else if (radix_ == 10 && str == 1 + start_ && str[-1] == '0'
                 && strchr ("xX", c))
        {
            if (lexer_->no_hexfloat_)
                lexer_->leave (lexer_->loc (start_, str),
                               "'strtod' does not support hexadecimal float");
            radix_ = 16;
        }
        else if (c == '.' && !dot_ && !ep_
                 && (radix_ == 10 || str - start_ != 2))
        {
            dot_ = str;
        }
        else if (ep_ && str == 1 + ep_ && strchr ("+-", c))
        { /* ok */ }
        else if (ep_ && isdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 10 && isdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 16 && isxdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 10 && strchr ("eE", c))
            ep_ = str;
        else if (!ep_ && radix_ == 16 && strchr ("pP", c))
            ep_ = str;
        else
            ok = false;

        return ok;
    }

    void check ()
    {
        int len = (int) std::strlen (text_);

        if (ep_ && !isdigit (end_[-1]))
            leave ("bad exponent in number token near '%s'",
                   text_ + (ep_ - lexer_->text_));

        if (len == 2 && strchr ("xX", end_[-1]))
            leave ("stray hex prefix in number token");

        if (dot_ && len == 1)
            leave ("'%s' is not a valid number", text_);
    }
};

class EP::Lexer::Identifier : public EP::Token
{
protected:
    Identifier (Lexer *lexer, char *text, const char *start, const char *end,
                Ident *id)
        : Token (lexer, text, start, end)
    {
        //info (lexer_->debug.tokenize, "%s", text);
        is_id = true;
        c_ = 'X';

        if (id && ident_)
            assert (id == ident_);
        else if (id)
            ident_ = id;
        else if (! ident_)
            ident_ = lexer->get_ident (text, std::nan(""), false);

        assert (text_);
        assert (ident_);
    }
public:
    static Token* make (Lexer*, const char*&);
    ~Identifier()
    {
        info (lexer_->debug.alloc, "~Identifier(%p) \"%s\"", this, text_);
    }
};

class EP::Lexer::Function : public EP::Lexer::Identifier
{
public:
    Function (Lexer *lexer, char *text, const char *start, const char *end,
              Ident *id)
        : Identifier (lexer, text, start, end, id)
    {
        assert (id);
        is_func = true;
        n_ary = 1;
        assert (ident_);
        info (lexer_->debug.alloc, "id=%p Function \"%s\"", ident_, text);
    }

    double value() const override
    {
        assert (ident_);
        assert (ident_->fn_);
        assert (n_ary == 1 && arg[0]);
        return (arg[0] ->* ident_->fn_) ();
    }

    bool check (int phase) const override
    {
        if (phase == 2)
        {
            assert (ident_);
            if (is_value && ident_->fn_ == & Token::poly
                && arg[0]->n_values() < 2)
            {
                const Token *t1 = arg[0]->nth_value (0, 1);
                const Token *from = this;
                const Token *to = t1->rightmost();
                auto loc = lexer_->loc (from->start_, to->end_ - 1);
                lexer_->leave (loc, "polynomial without coefficient");
            }
        }
        return true;
    }
};

auto EP::Lexer::set (int n, const char *name, Ident::FuncN f) -> Ident*
{
    Ident *id = get_ident (name, std::nan(""), true);
    info (debug.lookup, "set id=%p in %d %s(#%d) = @%p", id, ident_.size(),
          name, n, &f);

    id->fn_ = f;
    id->n_args = n;
    id->make_ =
    [] (EP::Lexer *l, char *t, const char *s, const char *e, EP::Ident *i)
        -> EP::Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new Function (l, t, s, e, i);
    };
    return id;
}


class EP::Lexer::Variable : public Lexer::Identifier
{
public:
    Variable (Lexer *lexer, char *text, const char *start, const char *end, Ident *id)
        : Identifier (lexer, text, start, end, id)
    {
        assert (text_);
        assert (ident_);
        is_var = is_value = true;
        n_ary = 0;

        //info (lexer_->debug.tokenize, "id=%p Variable %s", ident_, text_);
    }

    ~Variable()
    {
        info (lexer_->debug.alloc, "~Variable(%p) \"%s\"", this, text_);
    }

    static Ident::Maker maker;

    double value() const override
    {
        assert (ident_);
        if (ident_->undefined_)
            leave ("undefined variable '%s'", text_);
        //info (3, "%s = %f", text_, ident_->value_);
        return ident_->value_;
    }
};


EP::Ident::Maker EP::Lexer::Variable::maker =
    [] (EP::Lexer *l, char *t, const char *s, const char *e, EP::Ident *i)
        -> EP::Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new EP::Lexer::Variable (l, t, s, e, i);
    };


class EP::Lexer::Operator : public EP::Token
{
public:
    static constexpr const char* operators_ = "+-*/%^!(,)<=>?:;#";
    static constexpr const char* priority_ = "#(),;=?:>+*^!";
    Operator (EP::Lexer *lexer, const char*& caret)
        : Token (lexer, caret)
    {
        c_ = *caret++;
        is_op = true;
        if (c_ == '\0') c_ = '#';
        is_end = c_ == '#';
        if (!std::strchr (Operator::operators_, c_))
            fatal (c_ >= 0x20
                   ? "unknown operator char %d = 0x%02x = '%c'"
                   : "unknown operator char %d = 0x%02x", c_, c_, c_);

        if (c_ == '*' && *caret == '*')
            ++caret, c_ = '^';

        if (is_end)
        {
            end_ = 1 + start_;
            text_ = new char[1 + std::strlen ("<string end>")];
            std::strcpy (text_, "<string end>");
        }
        else
            alloc_text (caret);

        set_n_args_and_prio ();

        info (lexer_->debug.tokenize, "Operator \"%s\" #%d", text_, n_ary);
    }

    ~Operator()
    {
        info (lexer_->debug.alloc, "~Operator(%p) \"%s\"", this, text_);
    }

    double value() const override
    {
        for (int i = 0; i < n_ary; ++i)
            assert (arg[i]);

        switch (c_)
        {
            default:
                break;

            case '^': return pow ();
            case '?': return conditional ();
            case '=': return assignment ();
            case '!': return factorial ();
        }

        double v[3] = { std::nan(""), std::nan(""), std::nan("") };
        const double& a = v[0];
        const double& b = v[1];

        for (int i = 0; i < n_ary; ++i)
            v[i] = arg[i]->value();

        switch (c_)
        {
            default:
                fatal ("todo: \"%s\"", text_);
                assert (0);

            case '>' : return a > b;
            case '<' : return a < b;
            case '*' : return a * b;
            case '/' : return a / b;
            case '%' : return std::fmod (a, b);
            case '+' : return n_ary == 1 ? +a : a + b;
            case '-' : return n_ary == 1 ? -a : a - b;
            case ';' : return ((void) a, b);
        }

        return std::nan ("");
    }

    void set_n_args_and_prio ()
    {
        char c = c_;

        assert (c);
        is_post = strchr ("!", c_);
        is_pre = strchr ("", c_); // No static prefix operators for now.

        bool maybe_pre = !prev_ || (prev_->is_op && ! prev_->is_post && prev_->c_ != ')');
        if (c == '-' && maybe_pre)
            c = '*', is_pre = 1; // Unary minus is same as *.
        if (c == '+' && maybe_pre)
            c = '*', is_pre = 1; // Unary plus is same as *.
        if (c == '-')   c = '+';
        if (c == '/')   c = '*';
        if (c == '%')   c = '*';
        if (c == '<')   c = '>';

        n_ary = 0?0
            : is_pre || is_post ? 1
            : strchr ("()#", c_) ? 0
            : 2;

        is_infix = n_ary == 2 && ! strchr ("()?:,=", c_);

        const char *cs = Operator::priority_;
        prio = 1 + (int) (std::strchr (cs, c) - cs);
        assert (prio > 0);
    }

    bool check (int phase) const override
    {
        if (phase != 2 || !is_value)
            return true;

        for (int i = 0; i < n_ary; ++i)
        {
            assert (arg[i]);
            if (arg[i]->c_ == ',' && !strchr (",#", c_))
                leave ("'%s' cannot be argument of '%s'", arg[i]->text_, text_);
        }

        if (c_ == '?' && arg[1]->c_ != ':')
            leave ("'?' conditional must be followed by two expressions"
                   " separated by ':', have '%s'", arg[1]->text_);

        if (c_ == '=')
        {
            Lexer::Location loc = lexer_->loc (arg[0]->start_, start_);
            if (! arg[0]->is_var)
                lexer_->leave (loc, "assignment to function '%s', expect a"
                               " variable", arg[0]->text_);

            if (arg[0]->is_id)
            {
                auto var = arg[0];
                assert (var->ident_);
                if (var->ident_->built_in_)
                    lexer_->leave (loc, "assignment to built-in variable"
                                   " '%s = %.15f'", var->text_, var->value());
            }
        }

        return true;
    }

private:
    double conditional () const
    {
        assert (c_ == '?' && arg[1]->c_ == ':');
        bool sel = arg[0]->value();
        return arg[1]->arg[! sel]->value();
    }

    double assignment () const
    {
        assert (arg[0]->is_id);
        assert (arg[0]->ident_);
        assert (arg[0]->ident_->n_args == 0);

        arg[0]->ident_->undefined_ = false;
        return arg[0]->ident_->value_ = arg[1]->value();
    }

    double factorial () const
    {
        return std::tgamma (1 + arg[0]->value());
    }
};


auto EP::Lexer::Identifier::make (Lexer *lexer, const char*& str) -> Token*
{
    const char *start = str;
    while (isalpha (*str) || *str == '_' || (str != start && isdigit (*str)))
        ++str;
    auto len = str - start;
    char *text = new char[1 + len];
    std::strncpy (text, start, len);
    text[len] = '\0';

    info (lexer->debug.tokenize, "search \"%s\"", text);
    Ident *id = lexer->get_ident (text, std::nan(""), false);

    bool undefined = !id->make_;
    if (undefined)
    {
        // Unknown stuff can only be some variable like in "x = 10" or
        // some syntax error.
        id->make_ = Lexer::Variable::maker;
        id->undefined_ = true;
    }

    return (* id->make_) (lexer, text, start, str /* end */, id);
}


namespace
{

class _F_d : public EP::Lexer::Function
{
public:
    _F_d (EP::Lexer *lexer, char *text, const char *start, const char *end, EP::Ident *id)
        : Function (lexer, text, start, end, id)
    {
        assert (1 == n_ary);
        assert (id);
        assert (1 == id->n_args);
        assert (id->fx_);
        //info (lexer_->debug.alloc, "\"%s\"", text);
    }

    double value() const override
    {
        assert (ident_->fx_ && arg[0]);
        return ident_->fx_ (arg[0]->value());
    }
};

}; // ::anon

EP::Lexer::Lexer (const char *t, int f)
    : text_(std::strcpy (new char[1 + std::strlen (t)], t))
    , caret_(text_), flags_(f)
{
    debug.tokenize = f & FlagDebugTokenize;
    debug.parse = f & FlagDebugParse;
    debug.lookup = f & FlagDebugLookup;
    debug.alloc = f & FlagDebugAlloc;
    debug.tree = f & FlagDebugTree;

    info (debug.alloc, "new Lexer(%p)", this);

    init_ident();
    no_hexfloat_ = ! EP::supports_hexfloat ();
}


auto EP::Lexer::set_var (const char *name, double val) -> Ident*
{
    info (debug.lookup, "set_var in %d %s = %e", ident_.size(), name, val);

    Ident *id = get_ident (name, val, false);

    if (id->built_in_)
        leave (1, 0, id->n_args ? "%sfunction '%s'" : "%svariable '%s = %.15f'",
               "you cannot set built-in ", id->name_, id->value_);
    id->undefined_ = false;
    id->value_ = val;
    id->make_ = Variable::maker;
    return id;
}

auto EP::Lexer::set (const char *name, double val) -> Ident*
{
    Ident *id = get_ident (name, val, true);
    info (debug.lookup, "set id=%p in %d %s=%e", id, ident_.size(), name, val);

    id->make_ = Variable::maker;
    return id;
}

auto EP::Lexer::set (const char *name, Ident::Func f) -> Ident*
{
    Ident *id = get_ident (name, std::nan(""), true);
    info (debug.lookup, "set id=%p in %d %s(#1) = @%p", id, ident_.size(),
          name, f);
    id->fx_ = f;
    id->n_args = 1;
    id->make_ =
    [] (Lexer *l, char *t, const char *s, const char *e, Ident *i)-> Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new _F_d (l, t, s, e, i);
    };
    return id;
}

void EP::Lexer::init_ident()
{
    set ("pi",   3.1415926535897932384626433832795028841972);
    set ("e",    2.7182818284590452353602874713526624977572);
    set ("ln2",  0.6931471805599453094172321214581765680755);
    set ("ln10", 2.3025850929940456840179914546843642076011);
    set ("inf",  HUGE_VAL);
    set ("InF",  HUGE_VAL);
    set ("nan",  std::nan(""));
    set ("NaN",  std::nan(""));

    set ("sin", std::sin);
    set ("cos", std::cos);
    set ("tan", std::tan);
    set ("asin", "arcsin", std::asin);
    set ("acos", "arccos", std::acos);
    set ("atan", "arctan", std::atan);

    set ("sinh", std::sinh);
    set ("cosh", std::cosh);
    set ("tanh", std::tanh);
    set ("arsinh", "asinh", std::asinh);
    set ("arcosh", "acosh", std::acosh);
    set ("artanh", "atanh", std::atanh);

    set ("exp", std::exp);
    set ("log", "ln", std::log);
    set ("log2", std::log2);
    set ("log10", std::log10);

    set ("sqrt", std::sqrt);
    set ("cbrt", std::cbrt);
    set ("abs", std::fabs);
    set ("ceil", std::ceil);
    set ("floor", std::floor);
    set ("round", std::round);
    set ("mod1", [] (double x) -> double // Modulo 1 in [0, 1)
    {
        x -= std::floor (x);
        return x < 0 ? x+1 : x >= 1 ? x-1 : x;
    });

    set ("Gamma", std::tgamma);
    set ("gd", [] (double x) -> double // Gudermann
    {
        x = std::ldexp (x, -1);
        x = std::atan (std::tanh (x));
        return std::ldexp (x, 1);
    });

    // These functions actually have exactly 1 argument: a "," pack.
    set (2, "pow",   &Token::pow);
    set (2, "atan2", &Token::atan2);
    set (3, "sat",   &Token::sat);

    // Have 1 argument: a "," pack or something else.
    set (-1, "min", &Token::min);
    set (-1, "max", &Token::max);
    set (-1, "hypot", &Token::hypot);
    set (-1, "poly",  &Token::poly);
}

auto EP::parse (const char *text, int flags) -> Expression*
{
    Lexer *lexer = new Lexer (text, flags);

    lexer->tokenize ();

    info (lexer->debug.alloc, "Lexer(%p) #Tokens = %d, #Ident = %d", lexer,
          lexer->id_, lexer->ident_.size());

    Token *token = lexer->parse ();

    return new Expression (token);
};


EP::Lexer::~Lexer ()
{
    info (debug.alloc, "~Lexer(%p)", this);

    for (Token *t = tokens_; t; )
    {
        Token *next = t->next_;
        info (debug.alloc, "... ~Lexer(%p) ~%p \"%s\"", this, t, t->text_);
        delete t;
        t = next;
    }

    delete[] text_;
}

void EP::Lexer::leave_va (const char *text, int from, int to,
                          const char *fmt, va_list args)
{
    FILE* stream = stderr;

    fprintf (stream, "\nerror: ");

    vfprintf (stream, fmt, args);

    if (text && to >= from && to <= 1 + (int) std::strlen (text))
    {
        fprintf (stream, "\n%s", text);
        fprintf (stream, "\n%*s", from, "");
        while (from++ <= to)
            fprintf (stream, "^");
    }
    fprintf (stream, "\n");

    std::exit (-1);
}

void EP::Token::leave (const char *fmt, ...) const
{
    int from = (int) (start_ - lexer_->text_);
    int to   = (int) (end_ - 1 - lexer_->text_);

    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (lexer_->text_, from, to, fmt, args);
    va_end (args);
}

void EP::Lexer::leave (Location loc, const char *fmt, ...) const
{
    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (text_, loc.from, loc.to, fmt, args);
    va_end (args);
}


void EP::Lexer::leave (int from, int to, const char *fmt, ...) const
{
    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (text_, from, to, fmt, args);
    va_end (args);
}

void EP::Lexer::dump_sp (const Token *t) const
{
    if (t)
    {
        dump_sp (t->sp);
        printf (" %s%s", t->syntax(), t == sp ? "\n" : "");
    }
}

// Display Token like in the syntax table above.
const char* EP::Token::syntax () const
{
    return 0?0
        : is_list ? "L"
        : is_num  ? text_
        : is_var  ? text_
        : is_value ? "v"
        : is_func  ? text_
        : text_;
}

void EP::Lexer::push (Token *t, const char *s)
{
    if (debug.parse)
        out ("%s push %s => ", s, t->syntax());
    t->check (2);
    t->sp = sp;
    sp = t;
    if (debug.parse)
        dump_sp (sp);
}

auto EP::Lexer::pop () -> Token*
{
    assert (sp);
    if (debug.parse)
        out (" pop%s ", sp->syntax());
    Token *t = sp;
    sp = t->sp;
    t->sp = nullptr;
    return t;
}

int EP::Lexer::prio (const Token *t)
{
    return 0?0
        : t == nullptr  ? 0
        : t->is_value   ? 0
        : t->is_op ? t->prio
        : 0;
}


// Return TRUE if the top of the token stack has a specific sequence
// of tokens, each specified by a single character in string S.

bool EP::Lexer::expect (const char *s0, const Token *t0) const
{
    // Special chars and operators that do *not*  v -> v o v.
    const char *special = "()?:,=";
    const char *s = s0 + std::strlen (s0);

    for (const Token *t = t0 ? t0 : sp; t && s != s0; t = t->sp)
    {
#define FILTER_OUT(COND) \
    if (COND) return false; else continue

        switch (* --s)
        {
            case 'v': FILTER_OUT (! t->is_value);
            case 'X': FILTER_OUT (! t->is_var);
            case 'L': FILTER_OUT (! t->is_list);

            case 'F': FILTER_OUT (t->is_value || ! t->is_func);
            case '~': FILTER_OUT (t->is_value || ! t->is_pre);
            case '!': FILTER_OUT (t->is_value || ! t->is_post);
            case 'o': FILTER_OUT (t->is_value || ! t->is_infix);
            case '.': FILTER_OUT (0 /* nothing*/);
            case '#': FILTER_OUT (t->sp /* String start */);

            default:
                if (! std::strchr (special, *s))
                    fatal ("unrecognized expect '%c'", *s);

                FILTER_OUT (t->c_ != *s || t->is_value || t->is_list);
        } // switch
#undef FILTER_OUT
    }

    return s == s0;
}


/*
    Implement the following productions starting at S:

    Terminals:
        N = Number
        X = Variable identifier
        F = Function identifier
        # = End of Input

        Operators:
        o = Binary Operations: '+', '-', '*', '/', ';', '^', '**', ...
        ~ = Unary prefix Operations: '-', '+'.
        ! = Unary postfix Operations: '!'.
        = = Assignment '=' as binary Operation
        , = Seperate expressions by ',' to build Argument Lists.
        Special charactes like '(', ')', '#', '?', ':', ...

    Productions:
        S -> v #            Start
        L -> v              List of ',' separated values
        L -> L , v          List of ',' separated values

        v -> F ( L )        Function call
        v -> N              Number Terminal
        v -> X              Variable Terminal

        // Operators
        v -> ~ v            Unary prefix
        v -> v !            Unary postfix
        v -> v o v          Binary non-assignment, non-list
        v -> X = v          Binary assignment
        v -> v ? v : v      Ternary conditional

        v -> ( v )          Parenthesis for priority
*/

auto EP::Lexer::parse () -> Token*
{
    if (debug.parse)
        out ("parsing...\n");

    for (Token *t = tokens_; t; t = t->next_)
    {
        /* In almost all cases, Identifiers and Numbers alternate with
           Operators.  Do some sanity checks to spot syntax errors early
           on and to give more specific diagnostics.  */

        // v -> N
        // v -> X
        // v -> F (...
        if (! t->is_op)
        {
            assert (t->is_id || t->is_num);
            if (t->next_)
            {
                if (t->is_func && t->next_->c_ != '(')
                    t->leave ("expect '(' after function '%s'", t->text_);
                if (t->is_var && t->next_->c_ == '(')
                    t->leave ("'%s' is not a function", t->text_);
            }
        }

        if (0)
        {
reduce:;
            if (debug.parse)
                out ("-- ...\n");
        }

        diagnostic_after_push (t);

        // v -> v o v
        if (expect ("vov")
            && prio (t) <= sp->sp->prio)
        {
            if (debug.parse)
                out ("-- reduce v -> v %s v", sp->sp->text_);

            Token *v = pop();
            Token *o = pop();
            o->arg[0] = pop();
            o->arg[1] = v;
            o->is_value = 1;
            push (o);
            goto reduce;
        }

        // v -> X = v
        if (expect ("X=v")
            && prio (t) <= sp->sp->prio)
        {
            if (debug.parse)
                out ("-- reduce v -> '%s' = v", sp->sp->sp->text_);
            Token *v = pop();
            Token *o = pop();
            o->arg[0] = pop();
            o->arg[1] = v;
            o->is_value = 1;
            push (o);
            goto reduce;
        }

        // v -> ~ v
        if (expect ("~v")
            && prio (t) <= sp->sp->prio)
        {
            if (debug.parse)
                out ("-- reduce v -> %s v", sp->sp->text_);

            Token *v = pop();
            Token *o = pop();
            o->arg[0] = v;
            o->is_value = 1;
            push (o);
            goto reduce;
        }

        // v -> v !
        if (expect ("v") && expect ("!", t)
            && prio (sp->sp) <= t->prio)
        {
            if (debug.parse)
                out ("-- reduce v -> v %s", t->text_);
            t->arg[0] = pop();
            t->is_value = true;
            push (t);
            continue;
        }

        // v -> v ? v : v
        if (expect ("v?v:v")
            && prio (t) <= sp->sp->prio
            && prio (t) <= sp->sp->sp->sp->prio)
        {
            if (debug.parse)
                out ("-- reduce v -> v ? v : v");

            Token *v = pop();
            Token *o = pop();
            o->arg[0] = pop();
            o->arg[1] = v;
            o->is_value = 1;
            sp->arg[1] = o;
            o = pop();
            o->arg[0] = pop();
            o->is_value = 1;
            push (o);
            goto reduce;
        }

        // L -> v , v
        // L -> L , v

        if ((expect ("v,v") || expect ("L,v"))
            && prio (t) <= sp->sp->prio)
        {
            if (debug.parse)
                out ("-- reduce L -> %s , v", sp->sp->sp->syntax());
            Token *v = pop();
            Token *o = pop();
            o->arg[0] = pop();
            o->arg[0]-> is_value = 1;
            o->arg[1] = v;
            o->is_list = 1;
            push (o);
            goto reduce;
        }

        if (t->c_ == ')')
        {
            // v -> F ( v )
            // v -> F ( L )
            // v -> ( v )
            if (expect ("(v") || expect ("(L"))
            {
                if (debug.parse)
                {
                    if (expect ("F.."))
                        out ("-- reduce v -> '%s' ( %s )",
                             sp->sp->sp->text_, sp->syntax());
                    else
                        out ("-- reduce v -> ( %s )", sp->syntax());
                }

                Token *vL = pop();
                Token *pa = pop();

                if (expect ("F"))
                {
                    // v -> F ( v | L )
                    vL->is_value = 1;
                    Token *f = pop();
                    f->arg[0] = vL;
                    f->check_n_values ();
                    f->is_value = 1;
                    push (f);
                    continue;
                }

                if (vL->c_ == ',') // "L"
                    pa->leave ("expect function before argument list");

                // v -> ( v )
                vL->is_value = 1;
                push (vL);
                continue;
            }

            leave_missing_open_paren (t);
        }

        // '#' is the lowest-priority operator.  If it didn't collapse
        // everything by now, then it's a syntax error.
        if (t->is_end)
            break;

        // No more reductions: throw next token into the game.
        push (t, "--");
    } // for tokens

    diagnostic_after_parse ();

    if (debug.parse)
    {
        printf ("\n");
        info (true, "--parse success");
        sp->print();
        printf ("\n");
    }

    sp->finish();

    if (debug.tree)
        sp->dump();

    return sp;
}

void EP::Lexer::leave_missing_open_paren (const Token *t) const
{
    if (expect ("("))
        sp->leave (expect ("F.")
                   ? "empty argument list"
                   : "expect expression between '(...)'");

    for (const Token *y = sp; y; y = y->sp)
        if (y->c_ == '(')
            leave (loc (y->next_->start_, t->prev_->end_ - 1),
                   expect ("F.", y)
                   ? "bad argument list"
                   : "expect expression between '(...)'");

    t->leave ("missing opening '(' for '%s'", t->text_);
}


void EP::Lexer::diagnostic_after_push (const Token *t) const
{
    char c = t->c_;
    bool no_v = c == ')' ||  c == ',' || t->is_end || expect ("!", t) || expect ("o", t);

    // First, point out errors "after".

    if (expect ("o") && no_v)
        t->leave ("expect expression after binary '%s'", sp->text_);

    if (expect ("~") && no_v)
        t->leave ("expect expression after unary '%s'", sp->text_);

    if (expect ("=") && no_v)
        t->leave ("expect expression after '=' assignment");

    if (expect (",") && no_v)
        t->leave ("expect expression after ',' in argument list");

    if (expect ("?") && no_v)
        t->leave ("expect expression after '?' conditional");

    if (expect (":") && no_v)
        t->leave ("expect expression after ':' in conditional");

    if (expect ("v("))
        sp->leave ("misplaced '(' after '%s'", sp->sp->rightmost()->text_);

    if (expect ("!."))
    {
        const char *p = sp->sp->text_;
        if (expect ("("))
            sp->leave ("misplaced '(' after postfix '%s'", p);
        if (expect ("~"))
            sp->leave ("misplaced prefix '%s' after postfix '%s'", p);
        if (expect ("v"))
            sp->leftmost()->leave ("misplaced expression after postfix '%s'", p);
        if (expect ("F"))
            sp->leave ("misplaced function call after postfix '%s'", p);
    }

    if (expect ("v."))
    {
        const char *p = sp->sp->text_;
        if (expect ("v"))
            sp->leftmost()->leave ("misplaced expression after '%s' expression", p);
        if (expect ("F"))
            sp->leave ("misplaced function call after '%s' expression", p);
    }

    // Point out errors "before".

    if (expect ("o") && ! expect ("v."))
        sp->leave ("expect expression before binary '%s'", sp->text_);

    if (expect ("=") && ! expect ("X."))
        sp->leave ("expect variable before '=' assignment");

    if (c == ',' && ! expect ("L") && ! expect ("v"))
        t->leave ("expect expression before ',' in argument list");

    if (expect ("?") && ! expect ("v."))
        sp->leave ("expect expression before '?' conditional");

    if (expect (":") && ! expect ("?v."))
        sp->leave ("expect '?' conditional and expression before ':'");

    if ((! sp || sp->c_ == '(') && expect ("!", t))
        t->leave ("expect expression before postfix '%s'", t->text_);
}

void EP::Lexer::diagnostic_after_parse () const
{
    if (!sp)
        leave (0, (int) std::strlen (text_) - 1, "empty string");

    for (const Token *t = sp; t; t = t->sp)
        if (t->is_op && t->n_ary > 1 && ! t->is_value && t->c_ == ',')
            t->leave ("incomplete function call");
        else if (t->is_op && t->n_ary > 1 && ! t->is_value)
            t->leave ("incomplete %s'%s' expression",
                      t->c_ == '?'  ? "conditional "
                      : t->is_pre   ? "prefix "
                      : t->is_post  ? "postfix "
                      : t->is_infix ? "binary "
                      : "", t->text_);
        else if (std::strchr ("()", t->c_))
            t->leave ("unmatched '%s'", t->text_);

    if (! expect ("v"))
        sp->leave ("incomplete expression");

    if (sp->sp)
        sp->leave ("syntax error");
}


auto EP::Lexer::tokenize () -> Token*
{
    info (debug.tokenize, "tokenizing \"%s\" ...", text_);

    while (Token *token = lex())
    {
        if (debug.tokenize)
        {
            info (true, "new token: ");
            const char *n = 0?0
                : token->is_num ? "num"
                : token->is_id  ? "var"
                : token->is_end ? "#"
                : "op";
            (void) n;
            if (token->is_op)
                info (true, "%s[%d]!%d = '%s' ", n, token->n_ary,
                      token->prio, token->text_);
            else
                info (true, "%s = '%s' ", n, token->text_);

            if (debug.alloc)
                info (true, "= @%p", token);
            else
                info (true, "%s", "");
        }

        if (token->is_end)
            break;
    };

    if (debug.tokenize)
    {
        info (true, "done tokenize %d tokens:", id_);
        for (const Token *t = tokens_; t; t = t->next_)
            out ("  %s", t->text_);
        out ("\n\n");
    }

    if (debug.lookup)
    {
        info (true, "Identifiers %d:", ident_.size());
        for (const Tree::Node *n = ident_.first(); n; n = n->next())
        {
            const Ident& id = n->t_;
            const char* const f = 0?0
                : id.n_args == 1 ? "f(x)"
                : id.n_args == 2 ? "f(x,y)"
                : id.n_args == 3 ? "f(x,y,z)"
                : "f(...)";
            (void) f;
            if (id.n_args == 0)
                info (true, "id=%p '%s' = %e", &id, id.name_, id.value_);
            else if (id.n_args >= -1 && id.n_args <= 3)
                info (true, "id=%p '%s' = %s", &id, id.name_, f);
            else
                assert (0);
        }
    }

    if (! tokens_ || tokens_->is_end)
        leave (0, (int) (caret_ - text_ - 1), "no tokens in empty string");

    return tokens_;
}


auto EP::Lexer::lex () -> Token*
{
    while (std::isspace (*caret_))
        ++caret_;

    char c = *caret_;

    if (strchr (Operator::operators_, c))
        return new Operator (this, caret_);
    if (isalpha (c) || c == '_')
        return Identifier::make (this, caret_);
    if (isdigit (c) || c == '.')
        return new Number (this, caret_);

    leave (loc (caret_, caret_),
           c >= ' ' ? "%s '%c' = %d = 0x%02x" : "%s %d = 0x%02x",
           "spurious character",  c, c, c);
}


EP::Token::~Token ()
{
    info (lexer_->debug.alloc, "~Token(%p)", this);
    delete[] text_;
}


// Parsing is finished, check that whole tree.
void EP::Token::finish ()
{
    for (int i = 0; i < 2; ++i)
        assert (i < n_ary ? !! arg[i] : ! arg[i]);

    if (is_id && n_ary)
        assert (n_ary == 1);

    for (int i = 0; i < n_ary; ++i)
        arg[i]->finish();
}

void EP::Token::dump_id (const char *str) const
{
    (void) str;
    out ("%s{%d=", str, id_);
    if (!is_op)
        out ("%c", c_);
    else if (n_ary == 1 && is_pre)// && c_ == '(')
        out ("%c%c", c_, arg[0]->c_);
    else if (n_ary == 1 && is_post)// && c_ == '(')
        out ("%c%c", arg[0]->c_, c_);
    else if (n_ary == 2)
        out ("%c%c%c", arg[0]->c_, c_, arg[1]->c_);
    else
        assert (0);

    out ("}");
}

void EP::Token::dump () const
{
    out ("@%p:", this);
    dump_id();

    if (n_ary > 0)
        out (" [%d] ->", n_ary);
    for (int i = 0; i < n_ary; ++i)
        if (arg[i])
            arg[i]->dump_id (" ");
        else
            out (" (null)");
    out ("  |  \"%s\"\n", text_);

    for (int i = 0; i < n_ary; ++i)
        if (arg[i])
            arg[i]->dump();
}

void EP::Token::print (bool paren) const
{
    assert (is_value);

    const Token* const a = arg[0];
    const Token* const b = arg[1];

    if (paren) out ("(");

    if (is_num)
        out ("%s", text_);
    else if (is_id)
    {
        out ("%s", text_);
        if (a)
            a->print (1);
    }
    else if (is_op)
    {
        if (n_ary == 2)
        {
            bool pb = b->prio <= prio;
            pb |= next_ && next_->is_pre;
            assert (a && b);
            a->print (a->is_op && a->prio < prio);
            out (" %s ", text_);
            b->print (b->is_op & pb);
        }
        else if (n_ary == 1)
        {
            assert (a);
            if (is_pre) out ("%s%s", text_, a->prio > prio ? " " : "");
            a->print (a->is_op && (a->prio < prio || (next_ && next_->is_pre)));
            if (is_post) out ("%s", text_);
        }
        else
            assert (0);
    }
    else
        assert (0);

    if (paren) out (")");
}

double EP::Token::hypot2 () const
{
    double x;
    return c_ == ','
        ? arg[0]->hypot2() + arg[1]->hypot2()
        : (x = value(), x * x);
}

double EP::Token::pow () const
{
    auto& a = arg[0];
    auto& b = arg[1];
    assert (a && b);
    double x = a->value();
    double y = b->value();

    if (y == 2)
        return x * x;

    if (y == 0.5)
        return std::sqrt (x);

    if (a->is_var && 0 == strcmp ("e", a->text_))
        return std::exp (y);

    if (x == 2)
        return std::exp2 (y);

    return std::pow (x, y);
}

double EP::Token::sat () const
{
    assert (',' == c_);
    assert (',' == arg[0]->c_);
    double xx = arg[0]->arg[0]->value();
    double lo = arg[0]->arg[1]->value();
    double hi = arg[1]->value();
    return std::fmin (hi, std::fmax (lo, xx));
}

double EP::Token::poly () const
{
    double x = first_arg()->value();
    return horner (x, 0.0);
}

inline double EP::Token::horner (double x, double val) const
{
    return c_ == ','
        ? arg[0]->horner (x, x * val + arg[1]->value())
        : val;
}


auto EP::Token::nth_value (int n, int n_vals) const -> const Token*
{
    assert (n >= 0 && n < n_vals);

/*    if (c_ != ',')
    {
        assert (n == 0);
        return this;
    }*/
out ("\n\n");
    for (const Token *t = this; ; t = t->arg[0], --n_vals)
    {
        out ("%c %d/%d\n", t->c_, n,  n_vals);
        if (t->c_ != ',')
        {
            assert (n_vals == 1);
            return this;
        }
        else
        {
            assert (t->arg[0] && t->arg[1]);
            // ',' is left-associative.
            assert (t->arg[1]->c_ != ',');
            if (n == n_vals - 1)
                return t->arg[1];
        }
    }

    assert (0);
    return nullptr;
}

// Check number or argument for something that looks like a function call.
void EP::Token::check_n_values () const
{
    int n = arg[0]->n_values ();
    assert (ident_);
    if (n_ary == 0)
        leave ("variable \"%s\" cannot have arguments but has %d", text_, n);

    assert (n_ary == 1);

    if (ident_->n_args != -1 && ident_->n_args != n)
    {
        const Token *from, *to;

        if (n <= ident_->n_args)
        {
            const Token *t1 = arg[0]->nth_value (n - 1, n);
            from = this;
            to   = t1->rightmost();
        }
        else
        {
            const Token *t0 = arg[0]->nth_value (ident_->n_args, n);
            const Token *t1 = arg[0]->nth_value (n - 1, n);
            from = t0->leftmost()->prev_;
            to   = t1->rightmost();
        }

        auto loc = lexer_->loc (from->start_, to->end_ - 1);
        lexer_->leave (loc, "\"%s\" needs %d arguments but has %d", text_,
                       ident_->n_args, n);
    }
}
