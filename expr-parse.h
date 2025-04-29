#ifndef EXPR_PARSE_H
#define EXPR_PARSE_H

class ExpressionParser
{
public:
    struct Ident;
    class Token;
    class Lexer;
    class Expression
    {
        friend ExpressionParser;
        Expression (Token*);
        Expression (const Expression&) = delete;
        Expression (Expression&&) = delete;
        Token *expression_;
    public:
        // A copy of the text as passed to ExpressionParser::parse().
        const char* text() const;

        // Set the value of variable VAR to VALUE.  Free variables must have
        // been set at least once before value() is called to evaluate.
        void set (const char *var, double value) const;

        // Evaluate the expression.
        double value () const;

        // Short for:
        //      set ("x", x);
        //      value();
        double value (double x) const;

        // Decomission the expression and all resources attached to it.
        ~Expression();
    };

    enum Flag
    {
        FlagDebugTokenize = 1 << 0,
        FlagDebugParse = 1 << 1,
        FlagDebugLookup = 1 << 2,
        FlagDebugTree = 1 << 3,
        FlagDebugAlloc = 1 << 4,
        None = 0
    };

    /* Parse string TEXT and return an abstract representation that can be
       evaluated by means of expression->value() etc.  Supported operators
       are:

       Comment
            #               Ignore this and all following characters.

       Number literals:  Recognized numbers are the ones which the C languange
       would recognize like, for example
            1
            1.0
            2.34e-2         2.34 * 10^-2 = 0.0234
            0x0.ap5         Hexadecimal floating point.  Whether or not this is
                            actually supported depends on the capabilities of
                            the host's strtod implementation.  If not available,
                            the parser will terminate with an error.

       Binary infix:
            x + y           Addition
            x - y           Subtraction
            x * y           Multiplication
            x / y           Division
            x % y           Modulo
            x ^ y           Power, TeX style
            x ** y          Power, Python / GNUplot style
            x < y           Return 1 if x < y else 0
            x > y           Return 1 if x > y else 0

       Assignment:
            x = y           Set x to the value of y.  This can be used to
                            write expressions more neatly, for example
                                y = sin(x) + x,  y + 1 / y
                            is equivalent to
                                sin(x) + x + 1 / (sin(x) + x)

       Unary prefix:
            - x             Negation
            + x             Identity

       Unary postfix:
            x !             x Faktorial resp. Gamma (1 + x)

       Ternary:
            x ? y : z       Conditional return y if x != 0 else z.

       Function call:
            F (a, b, ...)   Number of accepted arguments depends on function

       Expression grouping:
            ( x )           For grouping operators like in  (1 + 2) * 3
                            which evaluates to 3 * 3 instead of 1 + 2 * 3
                            which evaluates to 1 + 6 because  *  has higher
                            priority than +.

            x; y; z; ...    Evaluates to the rightmost expression.
                            Only useful to set variable as side effects of
                            previous expression like in
                                a = x^2; y = x - a^2; a * y
                            which evaluates to x^2 * (x - (x^2)^2).

       All operators are left-associative and are evaluated from left to right.
       Operator priority, from lowest to highest:
            #               Comment
            ,               Separate arguments in function calls like pow
            ;               Sequence of expressions
            =               Assignment
            ? :             Ternary conditional
            < >             Comparison
            + -             Binary + and -
            + - * / %       Unary + and -, multiply, divide, modulo
            ^               Power
            !               Factorial
            ( )             Matching parenthesis

       Known constants
            pi              3.141592...  Ludolph's constant
            e               2.718281...  Euler's number
            ln2             0.693147...  Natural logarithm of 2
            ln10            2.302585...  Natural logarithm of 10
            Inf, inf        Infinity
            NaN, nan        Not-a-number

       Known functions
            sin (x)                 Sine
            cos (x)                 Cosine
            tan (x)                 Tangent
            asin (x), arcsin (x)    Arcus Sine
            acos (x), arccos (x)    Arcus Cosine
            atan (x), arctan (x)    Arrcus Tangent

            sinh (x)                Hyperbolic Sine
            cosh (x)                Hyperbolic Cosine
            tanh (x)                Hyperbolic Tangent
            arsinh (x), asinh (x)   Area Hyperbolic Sine
            arcosh (x), acosh (x)   Area Hyperbolic Cosine
            artanh (x), atanh (x)   Area Hyperbolic Tangent

            exp (x)                 Exponential Function, e^x
            log (x), ln (x)         Natural Logarithm (to base e)
            log2 (x)                Logarithm to Base 2
            log10 (x)               Logarithm to Base 10

            sqrt (x)                Square Root
            cbrt (x)                Cubic Root
            abs (x)                 Absolute Value
            ceil (x)                Round up to next integer
            floor (x)               Round down to next integer
            round (x)               Round to nearest integer
            mod1 (x)                Fractional Part (modulo 1 in [0, 1))

            Gamma (x)               Gamma Function
            gd (x)                  Gudermann's Function

            asinq (x)               asin(sqrt(x)) / sqrt(x)

       Two or more Arguments
            pow (x, y)              Power x ^ y
            atan2 (y, x)            The argument of complex number x + y * I
            sat (x, a, b)           Saturate x to the interval [a,b]

       Variable number of Arguments
            min (x, ...)            Minimum
            max (x, ...)            Maximum
            hypot (x, ...)          Hypothenuse of n-dimensional box
            poly (x, a0, a1, ...)   Polynomial in x:  a0 + a1 x + a2 x^2 + ...
    */
    static Expression* parse (const char *text, int flags = 0);

    // Parse TEXT, evaluate it with "x" = X, decommission all resources and
    // finally return the value.  For example,
    //
    //      ExpressionParser::eval ("y = 1 + pi; x * y", 123.0)
    //
    // will evaluate to  123.0 * (1 + pi).
    static double eval (const char *text, double x = 0.0, int flags = 0);

    // Whether or not the host C runtime's strtod function supports hexadecimal
    // float literals with base-2 exponent like 0x0.ap2 = 10*2^{-4} * 2^2 = 2.5.
    static bool supports_hexfloat ();
};

#endif // EXPR_PARSE_H
