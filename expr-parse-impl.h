#ifndef EXPR_PARSE_IMPL_H
#error Don't include this file, it's used by expr-parse-double.cpp and expr-parse-Float.cpp
#endif

#include "rb-tree.h"

#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cassert>
#include <cmath>
#include <cstdarg>

#include "diagnostic.h"

#define NORETURN __attribute__((__noreturn__))
#define DEPRECATED __attribute__((__deprecated__))

// Some templates to factor out difference between 'double' and 'Float'.
namespace
{
    // For diagnostics
    template<typename T> const char* _float_name () { return "double"; }

    // NaN and Inf
    template <typename F> F _nan() { return std::nan(""); }
    template <typename F> F _inf() { return HUGE_VAL; }

    // Gamma and factorial
    template <typename F> F _gamma (typename ExpressionParser<F>::Farg x)
    {
        return std::tgamma (x);
    }
    template <typename F>
    F _factorial (int, typename ExpressionParser<F>::Farg x)
    {
        return std::tgamma (1.0 + x);
    }

    double _m_e    = 2.7182818284590452353602874713526624977572;
    double _m_pi   = 3.1415926535897932384626433832795028841972;
    double _m_ln2  = 0.6931471805599453094172321214581765680755;
    double _m_ln10 = 2.3025850929940456840179914546843642076011;
    template <typename F> F _f_e()    { return _m_e; }
    template <typename F> F _f_pi()   { return _m_pi; }
    template <typename F> F _f_ln2()  { return _m_ln2; }
    template <typename F> F _f_ln10() { return _m_ln10; }
    template <typename F> F _f_prec() { return (F) __DBL_MANT_DIG__; }

    template<typename F> bool _supports_hexfloat ()
    {
        char *tail;
        (void) std::strtod ("0x1.2p0", &tail);
        return *tail == '\0';
    }

    template<typename F> bool _fixed_precision () { return true; }

    template<typename F> F _from_str (const char *str)
    {
        return std::strtod (str, nullptr);
    }
} // anon

// Other stuff.
namespace
{
    const char *const _nullstr = "<nullptr>";

    // Return int if s >= 0 and s < 1000, else -1.
    int _get_natural_value (const char *s)
    {
        int val = 0;
        for (int i = 0; i <= 3; ++i, ++s)
        {
            if (std::isdigit (*s))
                val = val * 10 + s[0] - '0';
            else
                return *s == '\0' ? val : -1;
        }
        return -1;
    }
} // anon

// Payload in a RBTree node that maps a string (identifier name) to
// a value or a function pointer.  If the value / function for an
// identifier is looked up, the found Ident's address will be stored
// in the Token.ident_ field for swift future lookup.
// This uses the feature that RBTree::Node's never change their payload.
// If a function requires exactly 1 argument, then it's pointer is f_.
// If a function requires more than 1 argument or accepts a variable
// number of arguments like "min" and "max", then it's a Token method
// pointer.

template<typename F>
struct ExpressionParser<F>::Ident
{
    typedef F (*FuncVoid)();
    typedef F (*Func)(Farg);
    typedef F (Token::* FuncN)() const;
    typedef Token* (*Maker)(Lexer*, char*, const char*, const char*, Ident*);

    const char *name_;
    FuncN fn_ = nullptr;
    Func fx_ = nullptr;
    FuncVoid f0_ = nullptr;
    Maker make_ = nullptr;

    F value_ { _nan<F>() };
    char n_args = -1;
    bool built_in_ = true;
    bool undefined_ = false;

    Ident (const char *str, Farg val, bool bi)
        : name_(str), value_(val), n_args(0), built_in_(bi) {}

    Ident (int n, const char *name, Maker m) : name_(name), make_(m), n_args(n)
    {}

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


template<typename F>
class ExpressionParser<F>::Lexer
{
    friend ExpressionParser<F>;
    friend Token;
    typedef RBTree<ExpressionParser<F>::Ident> Tree;

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

    auto search (const char* name) const -> const typename Tree::Node*
    {
        using S = typename Ident::Search;
        return ident_.template search< S, const char*> (name);
    }
    auto search (const char* name) -> typename Tree::Node*
    {
        using S = typename Ident::Search;
        return ident_.template search<S, const char*> (name);
    }

    typename Tree::Node* get_node (const char *name, Farg val, bool is_const)
    {
        auto node = ident_.add (Ident (name, val, is_const), true /* unique */);
        assert (node);
        return node;
    }

    Ident* get_ident (const char *name, Farg val, bool is_const)
    {
        bool is_x = name[0] == 'x' && name[1] == '\0';
        if (is_x && ident_x)
            return ident_x;

        typename Tree::Node *node = get_node (name, val, is_const);

        Ident *ident = & node->t_;
        if (is_x)
            ident_x = ident;

        return ident;
    }

    Ident* get_ident (const char *name, typename Ident::FuncVoid f0)
    {
        typename Tree::Node *node = get_node (name, _nan<F>(), true);

        Ident *ident = & node->t_;
        ident->f0_ = f0;

        return ident;
    }

    Ident* set (const char *name, Farg);
    Ident* set (const char *name, typename Ident::FuncVoid);
    Ident* set (const char *name, typename Ident::Func);
    Ident* set (int n, const char *name, typename Ident::FuncN);
    Ident* set_var (const char *name, Farg);

    void set (const char *a, const char *b, typename Ident::Func f)
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
    virtual ~Lexer ();
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


template<typename F>
class ExpressionParser<F>::Token
{
    friend ExpressionParser;
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
    F min () const
    {
        return c_ == ',' ? ::fmin (arg[0]->min(), arg[1]->min()) : value();
    }
    F max () const
    {
        return c_ == ',' ? ::fmax (arg[0]->max(), arg[1]->max()) : value();
    }
    F pow () const;
    F sat () const;
    F poly () const;
    F atan2 () const { return ::atan2 (arg[0]->value(), arg[1]->value()); }
    F hypot () const { return ::sqrt (hypot2()); }

    // ...and some helpers.
    F hypot2 () const;
    F horner (Farg, Farg) const;
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
    virtual F value () const = 0;
    virtual int get_natural () const { return -1; }

    F value (Farg x) const
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
            bool is_const: 1; // built-in, read-only.
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


template<typename F>
F ExpressionParser<F>::eval (const char *text, Farg x, int flags)
{
    Expression *ex = parse (text, flags);
    F val = ex->value (x);
    delete ex;
    return val;
}

// ExpressionParser::Expression is just a thin wrapper around Token in order
// to keep the header clean from all the gory bells and whistles.

template<typename F>
ExpressionParser<F>::Expression::Expression (Token *token)
    : expression_(token)
{}

template<typename F>
ExpressionParser<F>::Expression::~Expression ()
{
    // FIXME: We MUST NOT "delete expression_->lexer_;" because that may crash
    // as the compiler may use "expression_" AFTER the respective Token has
    // been deleted by the Lexer.  This works around a GCC bug, see
    // https://gcc.gnu.org/PR52339
    // https://stackoverflow.com/a/76172573/1556746
    Lexer *l = expression_->lexer_;
    delete l;
}

template<typename F>
F ExpressionParser<F>::Expression::value () const
{
    return expression_->value();
}

template<typename F>
F ExpressionParser<F>::Expression::value (Farg x) const
{
    return expression_->value (x);
}

template<typename F>
const char* ExpressionParser<F>::Expression::text () const
{
    return expression_->lexer_->text_;
}

template<typename F>
void ExpressionParser<F>::Expression::set (const char *name, Farg val) const
{
    expression_->lexer_->set_var (name, val);
}

// Register a new token.
template<typename F>
void ExpressionParser<F>::Lexer::add (Token *token)
{
    // Manage their administration in a light-weight list.

    if (! tokens_)
        tokens_ = token;

    if (prev_)
        prev_->next_ = token;

    token->prev_ = prev_;
    prev_ = token;
}

template<typename F>
ExpressionParser<F>::Token::Token (Lexer *lexer, char *text,
                                   const char* start, const char *end)
    : lexer_(lexer), text_(text), start_(start), end_(end), id_(lexer->id_++)
{
    is_ = 0;
    lexer_->add (this);
}


// End points to intended position of the string end '\0'.
template<typename F>
void ExpressionParser<F>::Token::alloc_text (const char *end)
{
    auto len = end - start_;
    text_ = new char[1 + len];
    std::strncpy (text_, start_, len);
    text_[len] = '\0';
    end_ = end;
}

template<typename F>
class ExpressionParser<F>::Lexer::Number : public ExpressionParser<F>::Token
{
    int radix_ = 10;
    int natural_value_ = -1;
    // Pointers into the original string or 0.
    const char *dot_ = nullptr;
    const char *ep_ = nullptr;
    F value_;

public:
    Number (ExpressionParser<F>::Lexer *lexer, const char*& str)
        : Token (lexer, str)
    {
        //info (lexer_->debug.tokenize, "...%s", str);
        // For the use of "this->", see
        // https://stackoverflow.com/q/1120833/1556746
        this->is_num = this->is_value = true;
        this->n_ary = 0;
        this->c_ = 'N';
        while (consumes (str))
            ++str;

        this->alloc_text (str);
        assert (this->text_ && * this->text_);
        value_ = _from_str<F> (this->text_);
        natural_value_ = _get_natural_value (this->text_);
        check();
    }

    int get_natural () const override
    {
        return natural_value_;
    }

    ~Number()
    {
        info (this->lexer_->debug.alloc, "~Number(%p) \"%s\"", this,
              this->text_);
    }

    F value() const override
    {
        return _fixed_precision<F>()
            ? value_
            : _from_str<F> (this->text_);
    }

    bool consumes (const char *str)
    {
        char c = *str;
        bool ok = true;

        if (c == '\0')
            ok = false;
        else if (radix_ == 10 && str == 1 + this->start_ && str[-1] == '0'
                 && std::strchr ("xX", c))
        {
            if (this->lexer_->no_hexfloat_)
                this->lexer_->leave (this->lexer_->loc (this->start_, str),
                                     "'strtod' does not support hexadecimal"
                                     " float");
            radix_ = 16;
        }
        else if (c == '.' && !dot_ && !ep_
                 && (radix_ == 10 || str - this->start_ != 2))
        {
            dot_ = str;
        }
        else if (ep_ && str == 1 + ep_ && std::strchr ("+-", c))
        { /* ok */ }
        else if (ep_ && std::isdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 10 && std::isdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 16 && std::isxdigit (c))
        { /* ok */ }
        else if (!ep_ && radix_ == 10 && std::strchr ("eE", c))
            ep_ = str;
        else if (!ep_ && radix_ == 16 && std::strchr ("pP", c))
            ep_ = str;
        else
            ok = false;

        return ok;
    }

    void check ()
    {
        int len = (int) std::strlen (this->text_);

        if (ep_ && ! std::isdigit (this->end_[-1]))
            this->leave ("bad exponent in number token near '%s'",
                         this->text_ + (ep_ - this->lexer_->text_));

        if (len == 2 && std::strchr ("xX", this->end_[-1]))
            this->leave ("stray hex prefix in number token");

        if (dot_ && len == 1)
            this->leave ("'%s' is not a valid number", this->text_);
    }
};

template<typename F>
class ExpressionParser<F>::Lexer::Identifier : public ExpressionParser<F>::Token
{
protected:
    Identifier (Lexer *lexer, char *text, const char *start, const char *end,
                Ident *id)
        : Token (lexer, text, start, end)
    {
        //info (lexer_->debug.tokenize, "%s", text);
        this->is_id = true;
        this->c_ = 'X';

        if (id && this->ident_)
            assert (id == this->ident_);
        else if (id)
            this->ident_ = id;
        else if (! this->ident_)
            this->ident_ = lexer->get_ident (text, _nan<F>(), false);

        assert (this->text_);
        assert (this->ident_);
    }
public:
    static Token* make (Lexer*, const char*&);
    ~Identifier()
    {
        info (this->lexer_->debug.alloc, "~Identifier(%p) \"%s\"", this,
              this->text_);
    }
};

template<typename F>
class ExpressionParser<F>::Lexer::Function : public ExpressionParser::Lexer::Identifier
{
public:
    Function (Lexer *lexer, char *text, const char *start, const char *end,
              Ident *id)
        : Identifier (lexer, text, start, end, id)
    {
        assert (id);
        this->is_func = true;
        this->n_ary = 1;
        assert (this->ident_);
        info (this->lexer_->debug.alloc, "id=%p Function \"%s\"",
              this->ident_, text);
    }

    F value() const override
    {
        assert (this->ident_);
        assert (this->ident_->fn_);
        assert (this->n_ary == 1 && this->arg[0]);
        return (this->arg[0] ->* this->ident_->fn_) ();
    }

    bool check (int phase) const override
    {
        if (phase == 2)
        {
            assert (this->ident_);
            if (this->is_value && this->ident_->fn_ == & Token::poly
                && this->arg[0]->n_values() < 2)
            {
                const Token *t1 = this->arg[0]->nth_value (0, 1);
                const Token *from = this;
                const Token *to = t1->rightmost();
                auto loc = this->lexer_->loc (from->start_, to->end_ - 1);
                this->lexer_->leave (loc, "polynomial without coefficient");
            }
        }
        return true;
    }
};

template<typename F>
auto ExpressionParser<F>::Lexer::set (int n, const char *name,
                                      typename Ident::FuncN f) -> Ident*
{
    Ident *id = get_ident (name, _nan<F>(), true);
    info (debug.lookup, "set id=%p in %d %s(#%d) = @%p", id, ident_.size(),
          name, n, &f);

    id->fn_ = f;
    id->n_args = n;
    id->make_ =
    [] (ExpressionParser::Lexer *l, char *t, const char *s, const char *e,
        ExpressionParser::Ident *i) -> ExpressionParser::Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new Function (l, t, s, e, i);
    };
    return id;
}


template<typename F>
class ExpressionParser<F>::Lexer::Variable : public Lexer::Identifier
{
public:
    Variable (Lexer *lexer, char *text, const char *start, const char *end,
              Ident *id)
        : Identifier (lexer, text, start, end, id)
    {
        assert (this->text_);
        assert (this->ident_);
        this->is_var = this->is_value = true;
        this->n_ary = 0;

        //info (lexer_->debug.tokenize, "id=%p Variable %s", ident_, text_);
    }

    ~Variable()
    {
        info (this->lexer_->debug.alloc, "~Variable(%p) \"%s\"", this,
              this->text_);
    }

    static typename Ident::Maker maker;

    F value() const override
    {
        assert (this->ident_);
        if (this->ident_->undefined_)
            this->leave ("undefined variable '%s'", this->text_);
        //info (3, "%s = %f", text_, ident_->value_);
        if (! _fixed_precision<F>() && this->ident_->f0_)
            return this->ident_->f0_();

        return this->ident_->value_;
    }
};


template<typename F>
typename ExpressionParser<F>::Ident::Maker
ExpressionParser<F>::Lexer::Variable::maker =
    [] (ExpressionParser::Lexer *l, char *t, const char *s, const char *e,
        ExpressionParser::Ident *i) -> ExpressionParser::Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new ExpressionParser::Lexer::Variable (l, t, s, e, i);
    };


template<typename F>
class ExpressionParser<F>::Lexer::Operator : public ExpressionParser::Token
{
public:
    static constexpr const char* operators_ = "+-*/%^!(,)<=>?:;#";
    static constexpr const char* priority_ = "#(),;=?:>+*^!";
    Operator (ExpressionParser::Lexer *lexer, const char*& caret)
        : Token (lexer, caret)
    {
        char& c = this->c_;
        c = *caret++;
        this->is_op = true;
        if (c == '\0') c = '#';
        this->is_end = c == '#';
        if (! std::strchr (Operator::operators_, c))
            fatal (c >= 0x20
                   ? "unknown operator char %d = 0x%02x = '%c'"
                   : "unknown operator char %d = 0x%02x", c, c, c);

        if (c == '*' && *caret == '*')
            ++caret, c = '^';

        if (this->is_end)
        {
            this->end_ = 1 + this->start_;
            this->text_ = new char[1 + std::strlen ("<string end>")];
            std::strcpy (this->text_, "<string end>");
        }
        else
            this->alloc_text (caret);

        set_n_args_and_prio ();

        info (this->lexer_->debug.tokenize, "Operator \"%s\" #%d", this->text_,
              this->n_ary);
    }

    ~Operator()
    {
        info (this->lexer_->debug.alloc, "~Operator(%p) \"%s\"", this,
              this->text_);
    }

    F value() const override
    {
        for (int i = 0; i < this->n_ary; ++i)
            assert (this->arg[i]);

        switch (this->c_)
        {
            default:
                break;

            case '^': return this->pow ();
            case '?': return this->conditional ();
            case '=': return this->assignment ();
            case '!': return this->factorial ();
        }

        F v[3] = { _nan<F>(), _nan<F>(), _nan<F>() };
        const F& a = v[0];
        const F& b = v[1];

        for (int i = 0; i < this->n_ary; ++i)
            v[i] = this->arg[i]->value();

        switch (this->c_)
        {
            default:
                fatal ("todo: \"%s\"", this->text_);
                assert (0);

            case '>' : return a > b;
            case '<' : return a < b;
            case '*' : return a * b;
            case '/' : return a / b;
            case '%' : return ::fmod (a, b);
            case '+' : return this->n_ary == 1 ? +a : a + b;
            case '-' : return this->n_ary == 1 ? -a : a - b;
            case ';' : return ((void) a, b);
        }

        return _nan<F>();
    }

    void set_n_args_and_prio ()
    {
        char c = this->c_;

        assert (c);
        this->is_post = std::strchr ("!", this->c_);
        // No static prefix operators for now.
        this->is_pre = std::strchr ("", this->c_);

        bool maybe_pre = ! this->prev_ || (this->prev_->is_op
                                           && ! this->prev_->is_post
                                           && this->prev_->c_ != ')');
        if (c == '-' && maybe_pre)
            c = '*', this->is_pre = 1; // Unary minus is same as *.
        if (c == '+' && maybe_pre)
            c = '*', this->is_pre = 1; // Unary plus is same as *.
        if (c == '-')   c = '+';
        if (c == '/')   c = '*';
        if (c == '%')   c = '*';
        if (c == '<')   c = '>';

        this->n_ary = 0?0
            : this->is_pre || this->is_post ? 1
            : std::strchr ("()#", this->c_) ? 0
            : 2;

        this->is_infix = this->n_ary == 2 && ! std::strchr ("()?:,=", this->c_);

        const char *cs = Operator::priority_;
        this->prio = 1 + (int) (std::strchr (cs, c) - cs);
        assert (this->prio > 0);
    }

    bool check (int phase) const override
    {
        if (phase != 2 || ! this->is_value)
            return true;

        for (int i = 0; i < this->n_ary; ++i)
        {
            assert (this->arg[i]);
            if (this->arg[i]->c_ == ',' && ! std::strchr (",#", this->c_))
                this->leave ("'%s' cannot be argument of '%s'",
                             this->arg[i]->text_, this->text_);
        }

        if (this->c_ == '?' && this->arg[1]->c_ != ':')
            this->leave ("'?' conditional must be followed by two expressions"
                         " separated by ':', have '%s'", this->arg[1]->text_);

        if (this->c_ == '=')
        {
            Lexer::Location loc = this->lexer_->loc (this->arg[0]->start_,
                                                     this->start_);
            if (! this->arg[0]->is_var)
                this->lexer_->leave (loc, "assignment to function '%s', expe"
                                     "ct a variable", this->arg[0]->text_);

            if (this->arg[0]->is_id)
            {
                auto var = this->arg[0];
                assert (var->ident_);
                if (var->ident_->built_in_)
                    this->lexer_->leave (loc, "assignment to built-in variable"
                                         " '%s = %.15f'", var->text_,
                                         (double) var->value());
            }
        }

        return true;
    }

private:
    F conditional () const
    {
        assert (this->c_ == '?' && this->arg[1]->c_ == ':');
        bool sel = this->arg[0]->value() != 0;
        return this->arg[1]->arg[! sel]->value();
    }

    F assignment () const
    {
        assert (this->arg[0]->is_id);
        assert (this->arg[0]->ident_);
        assert (this->arg[0]->ident_->n_args == 0);

        this->arg[0]->ident_->undefined_ = false;
        return this->arg[0]->ident_->value_ = this->arg[1]->value();
    }

    F factorial () const
    {
        return _factorial<F> (this->arg[0]->get_natural(),
                              this->arg[0]->value());
    }
};

template<typename F>
auto ExpressionParser<F>::Lexer::Identifier::make (Lexer *lexer,
                                                   const char*& str) -> Token*
{
    const char *start = str;
    while (isalpha (*str)
           || *str == '_'
           || (str != start && std::isdigit (*str)))
        ++str;
    auto len = str - start;
    char *text = new char[1 + len];
    std::strncpy (text, start, len);
    text[len] = '\0';

    info (lexer->debug.tokenize, "search \"%s\"", text);
    Ident *id = lexer->get_ident (text, _nan<F>(), false);

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

template<typename F>
class _F_d : public ExpressionParser<F>::Lexer::Function
{
public:
    _F_d (typename ExpressionParser<F>::Lexer *lexer, char *text,
          const char *start, const char *end,
          typename ExpressionParser<F>::Ident *id)
        : ExpressionParser<F>::Lexer::Function (lexer, text, start, end, id)
    {
        assert (1 == this->n_ary);
        assert (id);
        assert (1 == id->n_args);
        assert (id->fx_);
        //info (lexer_->debug.alloc, "\"%s\"", text);
    }

    F value() const override
    {
        assert (this->ident_->fx_ && this->arg[0]);
        return this->ident_->fx_ (this->arg[0]->value());
    }
};

}; // ::anon

template<typename F>
ExpressionParser<F>::Lexer::Lexer (const char *t, int f)
    : text_(std::strcpy (new char[1 + std::strlen (t)], t))
    , caret_(text_), flags_(f)
{
    if (t == _nullstr)
        leave (0, 0, "no tokens in %s", _nullstr);
    debug.tokenize = f & FlagDebugTokenize;
    debug.parse = f & FlagDebugParse;
    debug.lookup = f & FlagDebugLookup;
    debug.alloc = f & FlagDebugAlloc;
    debug.tree = f & FlagDebugTree;

    info (debug.alloc, "new Lexer(%p)", this);

    init_ident();
    no_hexfloat_ = ! _supports_hexfloat<F> ();
}


template<typename F>
auto ExpressionParser<F>::Lexer::set_var (const char *name, Farg val)
    -> Ident*
{
    info (debug.lookup, "set_var in %d %s = %e", ident_.size(), name,
          (double) val);

    Ident *id = get_ident (name, val, false);

    if (id->built_in_)
        leave (1, 0, id->n_args ? "%sfunction '%s'" : "%svariable '%s = %.15f'",
               "you cannot set built-in ", id->name_, (double) id->value_);
    id->undefined_ = false;
    id->value_ = val;
    id->make_ = Variable::maker;
    return id;
}

template<typename F>
auto ExpressionParser<F>::Lexer::set (const char *name, Farg val) -> Ident*
{
    Ident *id = get_ident (name, val, true);
    info (debug.lookup, "set id=%p in %d %s=%e", id, ident_.size(), name,
          (double) val);

    id->make_ = Variable::maker;
    return id;
}


template<typename F>
auto ExpressionParser<F>::Lexer::set (const char *name,
                                      typename Ident::FuncVoid f) -> Ident*
{
    Ident *id = _fixed_precision<F>()
        ? get_ident (name, f(), true /* is_const */)
        : get_ident (name, f);

    info (debug.lookup, "set id=%p in %d %s%s%e", id, ident_.size(), name,
          _fixed_precision<F>() ? "=" : "()~", (double) f());

    id->make_ = Variable::maker;
    return id;
}

template<typename F>
auto ExpressionParser<F>::Lexer::set (const char *name,
                                      typename Ident::Func f) -> Ident*
{
    Ident *id = get_ident (name, _nan<F>(), true);
    info (debug.lookup, "set id=%p in %d %s(#1) = @%p", id, ident_.size(),
          name, f);
    id->fx_ = f;
    id->n_args = 1;
    id->make_ =
    [] (Lexer *l, char *t, const char *s, const char *e, Ident *i)-> Token*
    {
        info (l->debug.alloc, "make \"%s\"", t);
        return new _F_d<F> (l, t, s, e, i);
    };
    return id;
}

template<typename F>
void ExpressionParser<F>::Lexer::init_ident()
{
    set ("pi",   _f_pi<F>);   // 3.141...
    set ("e",    _f_e<F>);    // 2.718...
    set ("ln2",  _f_ln2<F>);  // 0.693...
    set ("ln10", _f_ln10<F>); // 2.302...
    set ("inf",  _inf<F>());
    set ("InF",  _inf<F>());
    set ("nan",  _nan<F>());
    set ("NaN",  _nan<F>());

    static auto f_eps = []() { return ::nextafter (F{1.0}, F{2.0}) - F{1.0}; };
    set ("EPS",  f_eps);
    set ("PREC", _f_prec<F>);

    set ("sin", ::sin);
    set ("cos", ::cos);
    set ("tan", ::tan);
    set ("asin", "arcsin", ::asin);
    set ("acos", "arccos", ::acos);
    set ("atan", "arctan", ::atan);

    set ("sinh", ::sinh);
    set ("cosh", ::cosh);
    set ("tanh", ::tanh);
    set ("arsinh", "asinh", ::asinh);
    set ("arcosh", "acosh", ::acosh);
    set ("artanh", "atanh", ::atanh);

    set ("exp", ::exp);
    set ("exp2", ::exp2);
    set ("log", "ln", ::log);
    set ("log2", ::log2);
    set ("log10", ::log10);

    set ("sqrt", ::sqrt);
    set ("cbrt", ::cbrt);
    set ("abs", ::fabs);
    set ("ceil", ::ceil);
    set ("floor", ::floor);
    set ("round", ::round);
    set ("mod1", [] (Farg x) -> F // Modulo 1 in [0, 1)
    {
        F y = x - ::floor (x);
        return y < 0 ? y + F{1} : y >= 1 ? y - F{1} : y;
    });

    set ("Gamma", _gamma<F>);
    set ("gd", [] (Farg x) -> F // Gudermann
    {
        F y = ::ldexp (x, -1);
        y = ::atan (::tanh (y));
        return ::ldexp (y, 1);
    });

    set ("sinq", [] (Farg x) -> F // sin(sqrt(x)) / sqrt(x)
    {
        if (x * x == 0)
            return F{1} - x / F{6};
        F q = ::sqrt (::fabs (x));
        return x > 0
            ? ::sin (q) / q
            : ::sinh (q) / q;
    });

    set ("cosq", [] (Farg x) -> F // cos(sqrt(x))
    {
        F q = ::sqrt (::fabs (x));
        return x >= 0 ? ::cos (q) : ::cosh (q);
    });

    set ("asinq", [] (Farg x) -> F // asin(sqrt(x)) / sqrt(x)
    {
        if (x * x == 0)
            return F{1} + x / F{6};
        F q = ::sqrt (::fabs (x));
        return x > 0
            ? ::asin (q) / q
            : ::asinh (q) / q;
    });

    set ("atanq", [] (Farg x) -> F // atan(sqrt(x)) / sqrt(x)
    {
        if (x * x == 0)
            return F{1} - x / F{3};
        F q = ::sqrt (::fabs (x));
        return x > 0
            ? ::atan (q) / q
            : ::atanh (q) / q;
    });

    set ("artanhq", [] (Farg x) -> F // artanh(sqrt(x)) / sqrt(x)
    {
        if (x * x == 0)
            return F{1} + x / F{3};
        F q = ::sqrt (::fabs (x));
        return x > 0
            ? ::atanh (q) / q
            : ::atan (q) / q;
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

template<typename F>
auto ExpressionParser<F>::parse (const char *text, int flags) -> Expression*
{
    if (text == nullptr)
        text = _nullstr;

    Lexer *lexer = new Lexer (text, flags);

    lexer->tokenize ();

    info (lexer->debug.alloc, "Lexer(%p) #Tokens = %d, #Ident = %d", lexer,
          lexer->id_, lexer->ident_.size());

    Token *token = lexer->parse ();

    return new Expression (token);
};


template<typename F>
ExpressionParser<F>::Lexer::~Lexer ()
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

template<typename F>
void ExpressionParser<F>::Lexer::leave_va (const char *text, int from, int to,
                                           const char *fmt, va_list args)
{
    FILE* stream = stderr;

    fprintf (stream, "\nExpressionParser<%s>: error: ", _float_name<F>());

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

template<typename F>
void ExpressionParser<F>::Token::leave (const char *fmt, ...) const
{
    int from = (int) (start_ - lexer_->text_);
    int to   = (int) (end_ - 1 - lexer_->text_);

    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (lexer_->text_, from, to, fmt, args);
    va_end (args);
}

template<typename F>
void ExpressionParser<F>::Lexer::leave (Location loc,
                                        const char *fmt, ...) const
{
    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (text_, loc.from, loc.to, fmt, args);
    va_end (args);
}


template<typename F>
void ExpressionParser<F>::Lexer::leave (int from, int to,
                                        const char *fmt, ...) const
{
    va_list args;
    va_start (args, fmt);
    Lexer::leave_va (text_, from, to, fmt, args);
    va_end (args);
}

template<typename F>
void ExpressionParser<F>::Lexer::dump_sp (const Token *t) const
{
    if (t)
    {
        dump_sp (t->sp);
        printf (" %s%s", t->syntax(), t == sp ? "\n" : "");
    }
}

// Display Token like in the syntax table above.
template<typename F>
const char* ExpressionParser<F>::Token::syntax () const
{
    return 0?0
        : is_list ? "L"
        : is_num  ? text_
        : is_var  ? text_
        : is_value ? "v"
        : is_func  ? text_
        : text_;
}

template<typename F>
void ExpressionParser<F>::Lexer::push (Token *t, const char *s)
{
    if (debug.parse)
        out ("%s push %s => ", s, t->syntax());
    t->check (2);
    t->sp = sp;
    sp = t;
    if (debug.parse)
        dump_sp (sp);
}

template<typename F>
auto ExpressionParser<F>::Lexer::pop () -> Token*
{
    assert (sp);
    if (debug.parse)
        out (" pop%s ", sp->syntax());
    Token *t = sp;
    sp = t->sp;
    t->sp = nullptr;
    return t;
}

template<typename F>
int ExpressionParser<F>::Lexer::prio (const Token *t)
{
    return 0?0
        : t == nullptr  ? 0
        : t->is_value   ? 0
        : t->is_op ? t->prio
        : 0;
}


// Return TRUE if the top of the token stack has a specific sequence
// of tokens, each specified by a single character in string S.

template<typename F>
bool ExpressionParser<F>::Lexer::expect (const char *s0, const Token *t0) const
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

template<typename F>
auto ExpressionParser<F>::Lexer::parse () -> Token*
{
    if (debug.parse)
    {
        out ("parsing...\n");
        out ("-- %d tokens:", id_);
        for (const Token *t = tokens_; t; t = t->next_)
            out (" %s", t->text_);
        out ("\n");
    }

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
        out ("-- parse success: ");
        sp->print();
        out ("\n");
    }

    sp->finish();

    if (debug.tree)
        sp->dump();

    return sp;
}

template<typename F>
void ExpressionParser<F>::Lexer::leave_missing_open_paren (const Token *t) const
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


template<typename F>
void ExpressionParser<F>::Lexer::diagnostic_after_push (const Token *t) const
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

template<typename F>
void ExpressionParser<F>::Lexer::diagnostic_after_parse () const
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


template<typename F>
auto ExpressionParser<F>::Lexer::tokenize () -> Token*
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
        //bool fip = _fixed_precision<F>();
        info (true, "Identifiers %d:", ident_.size());
        for (const typename Tree::Node* n = ident_.first(); n; n = n->next())
        {
            const Ident& id = n->t_;
            const char* const f = 0?0
                : id.n_args == 1 ? "f(x)"
                : id.n_args == 2 ? "f(x,y)"
                : id.n_args == 3 ? "f(x,y,z)"
                : "f(...)";
            (void) f;
            if (id.n_args == 0)
            {
                double val = (double) (id.f0_ ? id.f0_() : id.value_);
                if (! id.f0_)
                    info (true, "id=%p '%s' = %e", &id, id.name_, val);
                else
                    info (true, "id=%p '%s()' ~ %e", &id, id.name_, val);
            }
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


template<typename F>
auto ExpressionParser<F>::Lexer::lex () -> Token*
{
    while (std::isspace (*caret_))
        ++caret_;

    char c = *caret_;

    if (strchr (Operator::operators_, c))
        return new Operator (this, caret_);
    if (isalpha (c) || c == '_')
        return Identifier::make (this, caret_);
    if (std::isdigit (c) || c == '.')
        return new Number (this, caret_);

    leave (loc (caret_, caret_),
           c >= ' ' ? "%s '%c' = %d = 0x%02x" : "%s %d = 0x%02x",
           "spurious character",  c, c, c);
}


template<typename F>
ExpressionParser<F>::Token::~Token ()
{
    info (lexer_->debug.alloc, "~Token(%p)", this);
    delete[] text_;
}


// Parsing is finished, check that whole tree.
template<typename F>
void ExpressionParser<F>::Token::finish ()
{
    for (int i = 0; i < 2; ++i)
        assert (i < n_ary ? !! arg[i] : ! arg[i]);

    if (is_id && n_ary)
        assert (n_ary == 1);

    for (int i = 0; i < n_ary; ++i)
        arg[i]->finish();
}

template<typename F>
void ExpressionParser<F>::Token::dump_id (const char *str) const
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

template<typename F>
void ExpressionParser<F>::Token::dump () const
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

template<typename F>
void ExpressionParser<F>::Token::print (bool paren) const
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

template<typename F>
F ExpressionParser<F>::Token::hypot2 () const
{
    F x;
    return c_ == ','
        ? arg[0]->hypot2() + arg[1]->hypot2()
        : (x = value(), x * x);
}

template<typename F>
F ExpressionParser<F>::Token::pow () const
{
    auto& a = arg[0];
    auto& b = arg[1];
    assert (a && b);
    F x = a->value();
    F y = b->value();

    if (y == 2)
        return x * x;

    if (y == 0.5)
        return ::sqrt (x);

    if (a->is_var && 0 == std::strcmp ("e", a->text_))
        return ::exp (y);

    if (x == 2)
        return ::exp2 (y);

    return ::pow (x, y);
}

template<typename F>
F ExpressionParser<F>::Token::sat () const
{
    assert (',' == c_);
    assert (',' == arg[0]->c_);
    F xx = arg[0]->arg[0]->value();
    F lo = arg[0]->arg[1]->value();
    F hi = arg[1]->value();
    return ::fmin (hi, ::fmax (lo, xx));
}

template<typename F>
F ExpressionParser<F>::Token::poly () const
{
    F x = first_arg()->value();
    return horner (x, F{0});
}

template<typename F>
inline F ExpressionParser<F>::Token::horner (Farg x, Farg val) const
{
    return c_ == ','
        ? arg[0]->horner (x, ::fma (x, val, arg[1]->value()))
        : val;
}


template<typename F>
auto ExpressionParser<F>::Token::nth_value (int n, int n_vals) const
    -> const Token*
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
template<typename F>
void ExpressionParser<F>::Token::check_n_values () const
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
