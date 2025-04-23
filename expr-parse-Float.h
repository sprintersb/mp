#ifndef EXPR_PARSE_FLOAT_H
#define EXPR_PARSE_FLOAT_H

#include "mp-float.h"

#include "expr-parse-double.h"

namespace
{
    template<> struct ExpressionParserFarg<Float> { using T = const Float&; };
};

extern template class ExpressionParser<Float>;

#endif // EXPR_PARSE_FLOAT_H
