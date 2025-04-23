#include "diagnostic.h"

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>

int (*_diagnostic_vfprintf) (FILE*, const char*, va_list);

bool info_function;
bool out_context;

enum { x_out, x_info, x_warning, x_error, x_fatal, x_None };

static const int exit_code[] =
{
    [x_out] = EXIT_SUCCESS,
    [x_info]  = EXIT_SUCCESS,
    [x_warning] = EXIT_SUCCESS,
    [x_error] = -2,
    [x_fatal] = -42,
    [x_None] = 0,
};

static FILE* const stream[] =
{
    [x_out] = stdout,
    [x_info]  = stderr,
    [x_warning] = stderr,
    [x_error] = stderr,
    [x_fatal] = stderr,
    [x_None] = nullptr
};

int _cont = x_None;

static void
diagnostic (int what, const char *name, const char *file, int line,
            const char *fmt, va_list args)
{
    FILE *s = stream[what];
    fflush (stdout);
    fflush (stderr);

    auto vfp = _diagnostic_vfprintf ? _diagnostic_vfprintf : vfprintf;

    if (_cont != what)
        if (out_context || what != x_out)
            fprintf (s, "%s%s:%d: %s: ", _cont != x_None ? "\n" : "",
                     file, line, name);

    size_t len = std::strlen (fmt);
    _cont = what != x_out && len && fmt[len - 1] == ' ' ? what : x_None;

    vfp (s, fmt, args);
    if (what != x_out && _cont == x_None)
        fprintf (s, "\n");
    fflush (s);

#ifdef __SANITIZE_ADDRESS__
    if (what == x_fatal)
    {
        // Raise a signal, which is nice with, say -fsanitize=address to
        // show a call trace.
        fprintf (s, "%s:%d: RAISING SIGNAL...\n", __FILE__, __LINE__);fflush(s);
        int * volatile p = 0;
        *p = 0;
    }
#endif // -fsanitize=address

    if (exit_code[what] != EXIT_SUCCESS)
        std::exit (exit_code[what]);
}

#if !defined (__GNUC__) || defined (__llvm__)
# error
#elif defined (__x86_64__) || defined (__i386__)
// Print raw symbol name.  Defined in gcc/config/i386/i386.c[c].
#define RAW "p"
#else
#error
#endif

#define DEF_ALIAS(X, Y)                                             \
    do {                                                            \
        __asm volatile (".ifndef " #X ".defined"           "\n\t"   \
                        ".set    " #X ".defined, 1"        "\n\t"   \
                        ".global %" RAW "0"                "\n\t"   \
                        ".set    %" RAW "[x], %" RAW "[y]" "\n\t"   \
                        ".endif"                                    \
                        : /* no outputs */                          \
                        : [x] "X" (X), [y] "X" (Y));                \
        /* Trigger a diagnostic if the prototypes don't match.  */  \
        __typeof__(Y) *y = X;                                       \
        (void) y;                                                   \
    } while (0)

void out_at (const char *file, int line, const char *func,
             const char *fmt, ...)
{
    DEF_ALIAS (outX_at, out_at);
    va_list args;
    va_start (args, fmt);
    static char name[300];

    *name = '\0';
    if (out_context)
        sprintf (name, "<%s>", func);

    diagnostic (x_out, name, file, line, fmt, args);

    va_end (args);
}

void info_at (const char *file, int line, const char *func,
              const char *fmt, ...)
{
    DEF_ALIAS (infoX_at, info_at);
    va_list args;
    va_start (args, fmt);
    static char name[300];

    if (_cont != x_info)
    {
        if (info_function)
            sprintf (name, "info <%s>", func);
        else
            sprintf (name, "info");
    }

    diagnostic (x_info, name, file, line, fmt, args);

    va_end (args);
}

void warning_at (const char *file, int line, const char *fmt, ...)
{
    DEF_ALIAS (warningX_at, warning_at);
    va_list args;
    va_start (args, fmt);

    diagnostic (x_warning, "warning", file, line, fmt, args);

    va_end (args);
}

void error_at (const char *file, int line, const char *fmt, ...)
{
    DEF_ALIAS (errorX_at, error_at);
    va_list args;
    va_start (args, fmt);

    diagnostic (x_error, "error", file, line, fmt, args);

    va_end (args);

    std::exit (-2); // Get rid of "noreturn function does return".
}

void fatal_at (const char *file, int line, const char *fmt, ...)
{
    DEF_ALIAS (fatalX_at, fatal_at);
    va_list args;
    va_start (args, fmt);

    diagnostic (x_fatal, "fatal", file, line, fmt, args);

    va_end (args);

    std::exit (-42); // Get rid of "noreturn function does return".
}

#if 0
// FIXME: This requires that we know the assembly names of the aliased
// function, which might be different, e.g. _error_at on Windows.
#define ALIAS(X) __attribute__((__alias__(#X)))
void warningX_at (const char*, int, const char *, ...) ALIAS (warning_at);
void errorX_at   (const char*, int, const char *, ...) ALIAS (error_at);
void fatalX_at   (const char*, int, const char *, ...) ALIAS (fatal_at);
void outX_at  (const char*, int, const char*, const char*, ...) ALIAS (out_at);
void infoX_at (const char*, int, const char*, const char*, ...) ALIAS (info_at);
#endif

char* _trace_str (void **p)
{
    char *str = (char*) malloc (1000);
    str[0] = '\0';

    for (char *s = str; *p; ++p, s += strlen (s))
        sprintf (s, " %p", *p);

    return str;
}

void _trace_write (void **p, const char *fname)
{
    FILE *fout = fopen (fname, "w");
    if (fout)
    {
        for ( ; *p; ++p)
            fprintf (fout, "%p ", *p);
        fclose (fout);
    }
}
