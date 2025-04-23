#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include <stdio.h>
#include <stdarg.h>

/* Set to a function like gmp_vfprintf from <gmp.h> to support "%Zd"
   to print mpz_t.  Or set to mpfr_vfprintf from <mpfr.h> to support "%RNf"
   to print mpfr_t.  Notice that gmp.h and mpfr.h declare these functions
   if <stdarg.h> or <cstdarg> was included prior to <gmp.h> resp. <mpfr.h>.
   Note that mpfr_printf also supports the print modifiers from gmp_printf.
   If not set or nullptr, use vfprintf. */
extern int (*_diagnostic_vfprintf) (FILE*, const char*, va_list);

#ifdef DIAGNOSTIC_NO_FORMAT_CHECK
__attribute__((__noreturn__))
extern void fatal_at (const char *file, int line, const char *fmt, ...);
__attribute__((__noreturn__))
extern void error_at (const char *file, int line, const char *fmt, ...);
extern void warning_at (const char *file, int line, const char *fmt, ...);
extern void info_at (const char *file, int line, const char *func,
                     const char *fmt, ...);
extern void out_at (const char *file, int line, const char *func,
                    const char *fmt, ...);
#else
__attribute__((__noreturn__)) __attribute__((__format__(__printf__,3,4)))
extern void fatal_at (const char *file, int line, const char *fmt, ...);
__attribute__((__noreturn__)) __attribute__((__format__(__printf__,3,4)))
extern void error_at (const char *file, int line, const char *fmt, ...);
__attribute__((__format__(__printf__,3,4)))
extern void warning_at (const char *file, int line, const char *fmt, ...);
__attribute__((__format__(__printf__,4,5)))
extern void info_at (const char *file, int line, const char *func,
                     const char *fmt, ...);
__attribute__((__format__(__printf__,4,5)))
extern void out_at (const char *file, int line, const char *func,
                    const char *fmt, ...);
#endif /* DIAGNOSTIC_NO_FORMAT */

extern bool info_function;
extern bool out_context;

#define error(x...) \
    error_at (__FILE__, __LINE__, ##x)

#define fatal(x...) \
    fatal_at (__FILE__, __LINE__, ##x)

#define warning(x...) \
    warning_at (__FILE__, __LINE__, ##x)

#define info(active, x...)                                             \
    do {                                                               \
        if ((active))                                                  \
            info_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);    \
    } while (0)

#define out(x...)                                               \
    do {                                                        \
        out_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);  \
    } while (0)

// In many cases, it's too intrusive to use DIAGNOSTIC_NO_FORMAT_CHECK
// because that ditches format checks for an entire translation unit.
// Thus, provide alternatives that don't check the format.

#define errorX(x...)   errorX_at (__FILE__, __LINE__, ##x)
#define fatalX(x...)   fatalX_at (__FILE__, __LINE__, ##x)
#define warningX(x...) warningX_at (__FILE__, __LINE__, ##x)

#define infoX(active, x...)                                             \
    do {                                                               \
        if ((active))                                                  \
            infoX_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);    \
    } while (0)

#define outX(x...)                                               \
    do {                                                        \
        outX_at (__FILE__, __LINE__, __PRETTY_FUNCTION__, ##x);  \
    } while (0)

__attribute__((__noreturn__))
extern void fatalX_at (const char *file, int line, const char *fmt, ...);
__attribute__((__noreturn__))
extern void errorX_at (const char *file, int line, const char *fmt, ...);
extern void warningX_at (const char *file, int line, const char *fmt, ...);
extern void infoX_at (const char *file, int line, const char *func,
                      const char *fmt, ...);
extern void outX_at (const char *file, int line, const char *func,
                     const char *fmt, ...);

#if defined (__GNUC__)

#include <string.h>
#include <stdlib.h>

// In order to show trace addresses as source file locations:
//
// * Complile the project, or at least all modules in the call tree, with -g2.
// * Run the program a.exe and record the trace addresses in a file trace.txt.
// * Run from a shell or by means of "make trace":
//      addr2line -p -e a.exe -a -f -i -s @trace.txt | python trace.py
// Apart from that one can also dig into a map file which does not require
// re-compiling or re-running or debug info.

extern void _trace_write (void**, const char*);
extern char* _trace_str (void**);

static inline __attribute__((__always_inline__,__artificial__))
void** _trace_bra (int _n, void **_p0, int _np)
{
    // __builtin_return_address() must be called with a compile-time
    // constant *sigh*.

#define _tr_1p_(N) \
    if (_n >= N && N < _np) \
        *_p++ = __builtin_return_address (N)
    void **_p = _p0;
    _tr_1p_(0); _tr_1p_(1); _tr_1p_(2); _tr_1p_(3); _tr_1p_(4);
    _tr_1p_(5); _tr_1p_(6); _tr_1p_(7); _tr_1p_(8); _tr_1p_(9);
    _tr_1p_(10); _tr_1p_(11); _tr_1p_(12); _tr_1p_(13); _tr_1p_(14);
    _tr_1p_(15); _tr_1p_(16); _tr_1p_(17); _tr_1p_(18); _tr_1p_(19);
#undef _tr_1p_

    *_p = nullptr;
    return _p0;
}

static inline __attribute__((__always_inline__,__artificial__))
char* trace_str (int _n)
{
    int _np = 20;
    void* _p[1 + _np];
    return _trace_str (_trace_bra (_n, _p, _np));
}

static inline __attribute__((__always_inline__,__artificial__))
void write_trace (int _n, const char *_fname)
{
    int _np = 20;
    void* _p[1 + _np];
    _trace_write (_trace_bra (_n, _p, _np), _fname);
}

#endif /* gcc */

#endif /* DIAGNOSTIC_H */
