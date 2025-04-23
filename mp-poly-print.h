#if !defined(MP_POLY_PRINT_H) || !defined(MP_POLY_CPP)
#error "Do not compile this file. It's included in mp-poly.cpp."
#endif

#include <cctype>

// Provide functionality of gmp_asprintf etc. but with additional
// format-checking.  We use gmp_asprintf because in some implementations,
// asprintf is not available or broken.  If we actually need %Zd and such
// formats that are GMP extensions, use plain `gmp_asprintf` etc. which
// does not check the format.
namespace
{
    __attribute__((__format__(__printf__,2,3)))
    int my_asprintf (char **pstr, const char *fmt, ...)
    {
        va_list args;
        va_start (args, fmt);
        int ret = gmp_vasprintf (pstr, fmt, args);
        va_end (args);
        return ret;
    }

    __attribute__((__format__(__printf__,2,3)))
    int my_sprintf (char *str, const char *fmt, ...)
    {
        va_list args;
        va_start (args, fmt);
        int ret = gmp_vsprintf (str, fmt, args);
        va_end (args);
        return ret;
    }
} // anon

namespace
{
    char* _E_to_pow10 (char*, const char* s_mul);

    template<class A>
    A _split_sign (const A& a, int&)
    {
        return a;
    }

    template<> Float _split_sign<Float> (const Float& a, int& s)
    {
        s = a.sgn();
        return a.abs();
    }

    template<> Ratio _split_sign<Ratio> (const Ratio& a, int& s)
    {
        s = a.sgn();
        return a.abs();
    }

    template<> double _split_sign<double> (const double& a, int& s)
    {
        s = a < 0 ? -1 : a > 0 ? 1 : 0;
        return std::abs (a);
    }

    const double _nan { std::nan("") };

    char* _double_as_uint64_t (char *str, double d)
    {
        static_assert (sizeof (double) == sizeof (uint64_t), "");
        static_assert (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__, "");

        uint64_t u64 = 0;
        if (d)
            std::memcpy (&u64, &d, sizeof (u64));
        if (str)
            my_sprintf (str, "0x%0*" PRIx64, u64 ? 16 : 1, u64);
        else
            my_asprintf (&str, "0x%0*" PRIx64, u64 ? 16 : 1, u64);

        return str;
    }

    uint32_t _to_uint32_t (float f)
    {
        static_assert (sizeof (float) == sizeof (uint32_t), "");
        static_assert (__FLOAT_WORD_ORDER__ == __ORDER_LITTLE_ENDIAN__, "");
        uint32_t u32 = 0;
        if (f)
            std::memcpy (&u32, &f, sizeof (u32));
        return u32;
    }

    char* _float_as_uint32_t (char *str, float f)
    {
        uint32_t u32 = _to_uint32_t (f);
        if (str)
            my_sprintf (str, "0x%0*" PRIx32, u32 ? 8 : 1, u32);
        else
            my_asprintf (&str, "0x%0*" PRIx32, u32 ? 8 : 1, u32);

        return str;
    }

    template<class A>
    void _print (std::ostream& ost, const A& a, int, typename Poly<A>::Style)
    {
        ost << a;
    }

    template<>
    void _print<double> (std::ostream& ost, const double& d, int,
                         Poly<double>::Style style)
    {
        char *str = nullptr;

        using S = Poly<double>::Style;
        switch (style)
        {
            default:
                ost << d;
                break;

            case S::TeX:        case S::Desmos:
            case S::TeXDown:    case S::DesmosDown:
                my_asprintf (&str, "%e", d);
                str = _E_to_pow10 (str, "\\cdot ");
                break;

            case S::ListFloat:
                my_asprintf (&str, "%e", d);
                break;

            case S::ListDouble:
            case S::ListDoubleDown:
                my_asprintf (&str, "<double:floatf>%f", d);
                break;

            case S::ListFloatHex:
                my_asprintf (&str, "%a", d);
                break;

            case S::ListDoubleX64:
                str = _double_as_uint64_t (nullptr, d);
                break;

            case S::ListFloatX32:
                str = _float_as_uint32_t (nullptr, (float) d);
                break;

            case S::ListVHDLFloatX32:
            {
                uint32_t u32 = _to_uint32_t ((float) d);
                my_asprintf (&str, "X\"%08" PRIx32 "\"", u32);
                break;
            }
        } // switch

        if (str)
            ost << str;
        std::free (str);
    }

    char* _ditch_trailing_0s (char* const str, const char* const e)
    {
        // FIXME: For C/C++ values like in Style::CHorner we do not want to
        // ditch "." and "e+00" at the same time because "1" is not a float
        // type in C/C++ whereas "1.0" and 1E0" are.

        char *dot = std::strchr (str, '.');
        char *pex = std::strchr (str, *e);
        if (dot && pex && pex > dot)
        {
            char *n = pex;
            while (n[-1] == '0') --n;
            n -= n[-1] == '.';
            if (n != pex)
                std::memmove (n, pex, 1 + std::strlen (pex));
        }

        pex = std::strstr (str, e);
        if (pex && pex[std::strlen (e)] == '\0')
            *pex = '\0';

        return str;
    }

    char* _best_print (const Float& f, const Float& f0, bool is_C)
    {
        char *str;

        if (f.is_integer()
            && mpfr_fits_sint_p (f.mpfr(), MPFR_RNDN))
        {
            my_asprintf (&str, "%d%s", (int) mpfr_get_si (f.mpfr(), MPFR_RNDN),
                         is_C ? ".0" : "");
            return str;
        }
    
        mpfr_asprintf (&str, "%RNE", f.mpfr());

        char *dot = std::strchr (str, '.');
        char *e   = std::strchr (str, 'E');
        int digs = (int) (e - dot) - 1;

        if (f0.get_precision() > f.get_precision()
            && dot && e && e > dot)
        {
            std::free (str);
            mpfr_asprintf (&str, "%.*RNE", digs, f0.mpfr());
        }

        // Ditch suffix "e+00" and trailing 0's that might have been produced.
        return _ditch_trailing_0s (str, "E+00");
    }

    // Transform powers of 10 from "E+02" to "10^2".
    // We perform this on the *string* representation.
    char* _E_to_pow10 (char * const str, const char* s_mul)
    {
        assert (std::isdigit (str[0]) || str[0] == '-');

        int pow10 = 0;
        char *e    = std::strchr (str, 'E');
        if (! e) e = std::strchr (str, 'e');
        if (e)
        {
            // Evaluate the E+09 or whatever exponent and snip it from str.
            pow10 = (int) std::strtol (e + 1, nullptr, 10);
            *e = '\0';
        }

        // Split the mantissa at '.' into head and tail, and snip '.' from str,
        char *dot = std::strchr (str, '.');
        const char *tail = dot ? dot + 1 : "";
        const char *head = str;
        if (dot)
            *dot = '\0';
        const char *sign = *head == '-' ? "-" : "";
        head += *head == '-';

        // Fuse head and tail to an integral mantissa.  As the decimal point
        // is after the tail now, we have to adjust pow10 accordingly.
        pow10 -= (int) std::strlen (tail);

        // Remove leading 0's from i0mant.
        char *s, *i0mant, *imant;
        my_asprintf (&i0mant, "%s%s", head, tail);
        for (imant = i0mant; imant[0] == '0' && imant[1]; ++imant) {}

        // Remove trailing 0's from imant.
        for (char *end = imant + std::strlen (imant) - 1;
             end > imant && *end == '0';
             --end)
        {
            *end = '\0';
            ++pow10;
        }

        const int l_imant = (int) std::strlen (imant);

        char *s_pow10 = nullptr;
        if (pow10 > 10)
            my_asprintf (&s_pow10, "10^{%d}", pow10);
        else if (pow10 > 0)
            my_asprintf (&s_pow10, "10^%d", pow10);
        else
            my_asprintf (&s_pow10, "10^{%d}", pow10 + l_imant - 1);

        // Now insert the '.' at a convenient place.

        if (0 == std::strcmp (imant, "0"))
            my_asprintf (&s, "0");
        else if (pow10 >= 4 && 0 == std::strcmp (imant, "1"))
            // "Large" positive pow10 and imant == "1".
            my_asprintf (&s, "%s%s", sign, s_pow10);
        else if (pow10 >= 4)
            // "Large" positive pow10.
            my_asprintf (&s, "%s%s%s%s", sign, imant, s_mul, s_pow10);
        else if (pow10 >= 0)
            // No '.' after imant, filling with few 0's.
            my_asprintf (&s, "%s%s%.*s", sign, imant, pow10, "00000");
        else if (pow10 > - l_imant)
            // Place '.' somewhere in imant.
            my_asprintf (&s, "%s%.*s.%s", sign, l_imant + pow10, imant,
                          imant + l_imant + pow10);
        else if (pow10 > -3 - l_imant)
            // Place '.' before imant, filling with few 0's.
            my_asprintf (&s, "%s0.%.*s%s", sign, -l_imant - pow10, "00000",
                          imant);
        // "Large" negative pow10 and ...
        else if (0 == std::strcmp (imant, "1"))
            // ... imant = 1.
            my_asprintf (&s, "%s%s", sign, s_pow10);
        else if (imant[1] == '\0')
            // ... no tail -> no '.'.
            my_asprintf (&s, "%s%c%s%s", sign, imant[0], s_mul, s_pow10);
        else
            // ... everything else.
            my_asprintf (&s, "%s%c.%s%s%s", sign, imant[0], imant + 1,
                         s_mul, s_pow10);

        //my_asprintf (&s, "<<%s[.]%sP{%d}=%sP{%d} -> %s>>",
        //              head, tail, pow10+(int)std::strlen(tail), imant, pow10, s);
        std::free (str);
        std::free (i0mant);
        std::free (s_pow10);

        return s;
    }

    int _f7_extra_flag = 0;

    template<>
    void _print<Float> (std::ostream& ost, const Float& a0, int prec2,
                        Poly<Float>::Style style)
    {
        char *str = nullptr;

        Float a;
        a.copy (a0);
        if (prec2 > 0)
            a.set_precision_round (prec2, MPFR_RNDN);

        using S = Poly<Float>::Style;
        switch (style)
        {
            default:
                str = _best_print (a, a0, false);
                break;

            case S::CHorner:
                str = _best_print (a, a0, true);
                break;

            case S::None:
            //case S::CHorner:
                a.print (ost, "%RNe");
                break;

            case S::TeX:      case S::Desmos:
            case S::TeXDown:  case S::DesmosDown:
                str = _E_to_pow10 (_best_print (a, a0, false), "\\cdot ");
                break;

            case S::Text:
            case S::TextDown:
                a.print (ost, "%RNg");
                break;

            case S::ListDouble:
            {
                double d = (double) a;
                int prec10 = 10; // 1 + prec2 * log2 (10);
                char fmt[10];
                sprintf (fmt, "%%.%df", prec10);
                my_asprintf (&str, fmt, d);
                break;
            }

            case S::ListFloatHex:
                mpfr_asprintf (&str, "%RNa", a.mpfr());
                break;

            case S::ListDoubleX64:
                str = _double_as_uint64_t (nullptr, (double) a);
                break;

            case S::ListFloatX32:
                str = _float_as_uint32_t (nullptr, (float) (double) a);
                break;

            case S::ListVHDLFloatX32:
            {
                uint32_t u32 = _to_uint32_t ((float) (double) a);
                my_asprintf (&str, "X\"%08" PRIx32 "\"", u32);
                break;
            }

            case S::ListF7:
            case S::ListF7Normalized:
            {
                char *str2;
                int n_bytes = 7;
                int expo, bytes[n_bytes];
                int flag = a.as_bytes (n_bytes, bytes, &expo);
                flag |= _f7_extra_flag;
                my_asprintf (&str, flag > 9 ? "0x%02x, " : "%d, ", flag);
                for (int i = 0; i < n_bytes; ++i)
                {
                    my_asprintf (&str2, "%s0x%02x,", str, bytes[i]);
                    std::free (str);
                    str = str2;
                }
                my_asprintf (&str, "%s %d", str2, expo);
                std::free (str2);
                break;
            }

            // Print like IEEE float in hex bytes, LSBs first, but with
            // an encoded mantissa of prec2 bits.  This is used in
            // math lookup tables in AVR-LibC.
            case S::FLT40tab:
            {
                const int n_expo = 8;
                const int n_bits = 1 + n_expo + prec2;
                const int n_bytes = n_bits / 8 + !!(n_bits & 7);
                int *bytes = a.as_IEEE (n_bytes, n_expo);

                for (int n = 0; n < n_bytes; ++n)
                {
                    char txt[20];
                    sprintf (txt, "%s0x%02x", n ? "," : "", bytes[n]);
                    ost << txt;
                }

                delete[] bytes;
                str = nullptr;
                break;
            }
        } // switch

        if (str)
            ost << str;
        std::free (str);
    }

    // xxx_printf ("%RNe", ... ) choses a number of figures to print after
    // the decimal point.  However, for that same number of figures, we can get
    // a number with less error.  For example, printing mpq_t with intermediate
    // mpfr_t will print for 2/3 and with 30 bits
    //      "6.6666666698e-01"
    // where with the same number of figures we could get
    //      "6.6666666667e-01"
    // which is closer.
    //
    // The returned char* should be released with [std::]free().

    char* _best_print (const Ratio& q, int bits, bool)
    {
        char *str;

        if (q.is_integer() && mpz_fits_sint_p (q.numref()))
        {
            my_asprintf (&str, "%d", (int) mpz_get_si (q.numref()));
            return str;
        }

        mpfr_asprintf (&str, "%RNE", Float { q, 1 + bits }.mpfr());

        char *dot = std::strchr (str, '.');
        char *e   = std::strchr (str, 'E');
        int digs = (int) (e - dot) - 1;

        if (dot && e && e > dot)
        {
            std::free (str);
            mpfr_asprintf (&str, "%.*RNE", digs, Float { q, 18 + bits }.mpfr());
        }

        // Ditch suffix "e+00" and trailing 0's that might have been produced.
        return _ditch_trailing_0s (str, "E+00");
    }

    template<>
    void _print<Ratio> (std::ostream& ost, const Ratio& a, int prec2,
                        Poly<Ratio>::Style style)
    {
        prec2 = prec2 > 0 ? prec2 : 53;

        using S = Poly<Ratio>::Style;
        auto float_style = (Poly<Float>::Style) style;
        char *str = nullptr;

        switch (style)
        {
            default:
                ost << a;
                break;

            case S::TeXFrac:
            case S::TeXFracDown:
                if (a.is_integer())
                    ost << a;
                else
                    gmp_asprintf (&str, "%s\\frac{%Zd}{%Zd}",
                                  a < 0_Q ? "-" : "",
                                  a.abs().numref(), a.denref());
                break;

            case S::List:
                ost << a;
                break;

            case S::ListDouble:
            case S::ListDoubleDown:
                ost << "<ListF>" << a;
                break;

            case S::ListFloat:
                str = _best_print (a, prec2, false);
                break;

            case S::ListFloatHex:
                mpfr_asprintf (&str, "%RNa", Float { a, prec2 }.mpfr());
                break;

            case S::ListDoubleX64:
            {
                // Converting mpq_t to double truncates to zero,
                // therefore use mpfr_t as intermediate.
                _print<Float> (ost, Float { a, 8 + 53 }, 0, float_style);
                break;
            }
        } // switch

        if (str)
            ost << str;
        std::free (str);
    }

    struct Horner
    {
        int deg = 0;
        int n_paren = 0;
    };

} // ::anon


template<class A>
std::ostream& Poly<A>::print_horner (std::ostream& ost, Style style,
                                     const char* var, int prec2,
                                     const char *l, const char *r) const
{
    Horner horn;

    auto _ = [&ost, prec2, style, l, r](const A& a) -> const char*
    {
        ost << l;
        _print<A> (ost, a, prec2, style);
        ost << r;
        return "";
    };

    if (deg_ < 1)
        return ost << _(at(0));

    int n_out = 0;
    for (int i = 0; i <= deg_; ++i, ++n_out)
    {
        A ai = at(i);

        if (_is_zero<A>(ai))
        {
            --n_out;
            continue;
        }

        bool paren = false;
        bool is_one = false;
        int sign = 1;

        if (i == deg_)
        {
            ai = _split_sign<A> (ai, sign);
            is_one = ai == Scalar::one();
        }

        if (n_out)
            ost << (sign < 0 ? " - " : " + ");
        else if (sign < 0 && is_one)
            ost << "-";

        bool mul_neg = !n_out && sign < 0 && !is_one;
        bool mul_pos = i == deg_ && !is_one;

        if (i > horn.deg)
        {
            int ex = i - horn.deg;

            if (style == CHorner || ex <= 2)
            {
                ost << var;
                for (int k = 1; k < ex; ++k)
                    ost << "*" << var;
            }
            else if (style == PythonHorner
                     || style == GnuplotHorner)
            {
                ost << var;
                if (ex > 1)
                    ost << "**" << ex;
            }
            else if (style == PythonHornerPow)
            {
                if (ex == 1)
                    ost << var;
                else if (ex > 1)
                    ost << "pow(" << var << "," << ex << ")";
            }
            else
                fatal ("todo %d: %s", style, __PRETTY_FUNCTION__);

            horn.deg = i;
            if (is_one)
                continue;
            ost << " * ";
            paren = !mul_pos || mul_neg;
        }

        ost << (paren ? "(" : "");
        horn.n_paren += paren;

        if (mul_neg)
            ost << "-";
        ost << _(ai);
    }

    while (horn.n_paren--)
        ost << ")";

    return ost;
}

template<class A>
std::ostream& Poly<A>::print_VHDLtab (std::ostream& ost,
                                      const char* var, int prec2) const
{
    auto _var = [&ost, var]() -> void
    {
        for (const char *s = var; *s && *s != '#'; ++s)
            ost << *s;
    };

    auto _ = [&ost, prec2](const A& x) -> void
    {
        Float a0{x};
        Float a;
        a.copy (a0);
        if (prec2 > 0)
            a.set_precision_round (prec2, MPFR_RNDN);

        char *str = _best_print (a, a0, true);
        ost << str;
        std::free (str);
    };

    int n_used = 1 + deg();

    const char *hash = strchr (var, '#');
    if (hash)
    {
        sscanf (1 + hash, "%d", &n_used);
    }

    int n = _next_pow2 (n_used);

    for (int i = 0; i < n; ++i)
    {
        A ai = (i <= deg_) ? at(i) : Scalar::zero();

        ost << "\t";
        _print<A> (ost, ai, prec2, ListVHDLFloatX32);

        // Optional "," after VHDL hex X"beefbabe" and start VHDL comment.
        ost << ", "[i == n - 1] << " -- ";

        // Print VHDL comment
        if (i < n_used)
        {
            int sign;
            A bi = _split_sign<A> (ai, sign);
            _var();
            ost << i << "\t" << "+-"[sign < 0];
            _(bi);
        }
        else
            ost << "unused";

        ost << std::endl;
    }

    return ost;
}

template<class A>
void Poly<A>::print_VHDL_Table (std::ostream& ost, const char *var,
                                const char *typ, int prec2, int n_coeff) const
{
    char var0_n[10];
    sprintf (var0_n, "%c#%d", var[0], n_coeff);
    if (n_coeff < 1)
        var0_n[1] = '\0';

    ost << "CONSTANT " << var << " : " << typ << " :="
        << " -- Werte in FP (IEEE-754 single)" << std::endl;
    ost << "(" << std::endl;
    this->print (ost, Style::VHDLTab, var0_n, prec2);
    ost << ");" << std::endl;
}

/* Like used in AVR-LibC's math lookup tables like

.L_table:
    .byte   2
    .byte        0x89,0x88,0x08,0x3c    ; 1/120
    .byte   0xab,0xaa,0xaa,0x2a,0x3e    ; 1/6
    .byte   0x00,0x00,0x00,0x80,0x3f    ; 1

   Coefficient [0] comes last.  */
template<class A>
std::ostream& Poly<A>::print_FLT40tab (std::ostream& ost,
                                       const char* start) const
{
    if (! start)
        start = "\t.byte\t";

    ost << start << deg_ << std::endl;

    for (int i = deg_; i >= 0; --i)
    {
        A ai = at(i);

        ost << start << (i == deg_ ? "     " : "");
        _print<A> (ost, ai, i == deg_ ? 23 : 31, FLT40tab);
        char cmt[30];
        sprintf (cmt, "\t; % .10f", (double) ai);

        ost << cmt << std::endl;
    }

    return ost;
}

template<class A>
std::ostream& Poly<A>::print (std::ostream& ost, Style style, const char* txt,
                              int prec2) const
{
    bool list = style >= Style::List;
    bool isTeX = 0, isGpl = 0, isHTM = 0, isMinus = 0;
    bool horner = 0, nosplit = 0, down = 0, todo = 0;
    const A _0 = Scalar::zero();
    const A _1 = Scalar::one();

    using T = struct { bool& var; std::initializer_list<Style> styles; };
    const T ts[] =
    {
        { isTeX,  { TeX, TeXDown, TeXFrac, TeXFracDown, Desmos, DesmosDown } },
        { isHTM,    { HTML, HTMLDown, HTMLMinus, HTMLMinusDown } },
        { isMinus,  { HTMLMinus, HTMLMinusDown } },
        { isGpl,    { Gnuplot, GnuplotDown } },
        { horner,   { GnuplotHorner, PythonHorner, PythonHornerPow, CHorner } },
        { nosplit,  { None } },
        { down,     { TeXDown, TeXFracDown, GnuplotDown, HTMLDown, TextDown,
              HTMLMinusDown, DesmosDown, ListDoubleDown } },
        { todo,     { } },
    };
    for (auto& t : ts)
        for (auto& s : t.styles)
            t.var |= style == s;

    if (todo)
        fatal ("todo %d: %s", style, __PRETTY_FUNCTION__);

    bool skip0 = !list;
    bool split = !list && !nosplit && !horner;

    const char *var = txt ? txt : list ? ", " : "x";

    if (style == ListF7 || style == ListF7Normalized)
        if (!txt)
            txt = "F7_CONST_DEF (X, @)\n|";

    const char* l = "";
    const char* r = "";
    char t[1 + strlen (txt ? txt : "")];
    *t = '\0';

    if (txt)
    {
        std::strcpy (t, txt);
        char *at  = strchr (t, '@');
        char *bar = strchr (t, '|');
        if (at && bar && bar > at)
        {
            var = 1 + bar;
            l = t;
            r = 1 + at;
            *at = *bar = '\0';
        }
    }

    if (horner)
        return print_horner (ost, style, var, prec2, l, r);
    else if (style == VHDLTab)
        return print_VHDLtab (ost, txt, prec2);
    else if (style == FLT40tab)
        return print_FLT40tab (ost, txt);

    const char* const times = 0?0
        : isTeX ? " "
        : isGpl ? "*"
        : isHTM ? " "
        : " * ";

    const char* s_minus = 0?0
        : isMinus ? "&minus;"
        : "-";

    assert (!list || !split);

    auto _ = [&ost, prec2, style, l, r](const A& a) -> void
    {
        ost << l;
        _print<A> (ost, a, prec2, style);
        ost << r;
    };

    int n_out = 0;
    for (int j = 0; j <= deg_; ++j, ++n_out)
    {
        int i = down ? deg_ - j : j;
        int sign = 1;

        A ai = (split && n_out) ? _split_sign<A> (at(i), sign) : at(i);

        if (skip0 && _is_zero<A>(ai))
        {
            --n_out;
            continue;
        }

        if (split && !n_out && isMinus)
        {
            ai = _split_sign<A> (at(i), sign);
            ost << s_minus;
        }

        if (n_out)
        {
            if (list)
                // The separator between coefficients.
                ost << var;
            else
                // The separator between monomes.
                ost << " " << (sign >= 0 ? "+" : s_minus) << " ";
        }

        if (style == ListF7Normalized
            && j == deg_ - 1
            && at(deg_) == _1)
        {
            _f7_extra_flag = 8;
            _(ai);
            break;
        }

        if (list)
            _(ai);
        else
        {
            bool is_one = style != None && ai == _1;

            if (split && i && !n_out && ai < _0)
            {
                // If the 1st term is something like -x^2 but not -<number>.
                A bi = _split_sign<A> (at(i), sign);
                if ((is_one = bi == _1))
                    ost << s_minus;
            }

            if (i == 0)
                _(ai);
            else if (is_one)
                ost << var;
            else
            {
                _(ai);
                ost << times << var;
            }

            // The exponent.
            if (i >= 2)
            {
                const char* const l_pow = 0?0
                    : isTeX && i <= 9 ? "^"
                    : isTeX ? "^{"
                    : isGpl ? "**"
                    : isHTM ? "<sup>"
                    : "^";

                const char* const r_pow = 0?0
                    : isTeX && i <= 9 ? ""
                    : isTeX ? "}"
                    : isHTM ? "</sup>"
                    : "";

                ost << l_pow << i << r_pow;
            }
        }
    }

    if (!n_out)
        _(_0);

    return ost;
}
