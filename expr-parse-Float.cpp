#include "expr-parse-Float.h"

#define EXPR_PARSE_IMPL_H
#include "expr-parse-impl.h"
#undef EXPR_PARSE_IMPL_H

// Specializations for Float.
namespace
{
    // For diagnostics
    template<> const char* _float_name<Float> () { return "Float"; }

    // NaN and Inf
    template<> Float _nan<Float>() { return Float(); }
    template<> Float _inf<Float>() { return "inf"_Float; }

    // Gamma and factorial
    template<> Float _gamma<Float> (const Float& x) { return x.gamma(); }
    template<> Float _factorial<Float> (int n, const Float& x)
    {
        return n >= 0 && n <= 200
            ? Float::factorial (n)
            : (1_R + x).gamma();
    }

    template<> Float _f_e<Float>()    { return Float::e(); }
    template<> Float _f_pi<Float>()   { return Float::pi(); }
    template<> Float _f_ln2<Float>()  { return Float::ln2(); }
    template<> Float _f_ln10<Float>() { return Float::ln10(); }
    template<> Float _f_prec<Float>()
    {
        return Float { (double) Float::GetPrecision (2) };
    }

    template<> bool _supports_hexfloat<Float> () { return true; }

    template<> bool _fixed_precision<Float> () { return false; }

    template<> Float _from_str<Float> (const char *str)
    {
        return Float::from_str (str);
    }
} // anon

// Explicit instanciation.
template class ExpressionParser<Float>;

#include <type_traits>

static_assert (std::is_lvalue_reference<ExpressionParser<Float>::Farg>::value,
               "");
