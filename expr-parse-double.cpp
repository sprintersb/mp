#include "expr-parse-double.h"

#define EXPR_PARSE_IMPL_H
#include "expr-parse-impl.h"
#undef EXPR_PARSE_IMPL_H

// Explicit instanciation.
template class ExpressionParser<double>;
