#include "K.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Utility routines for math operations                          */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*
** Static arrays that hold precomputed values for the
** math functions.  These are all filled by one call to 
** init_math().
*/

static KScalar a_pow_half       [pow_half_VALS+1] = { -1.0 };
static int     done_pow_half = 0;  /* set to 1 in init_math */
static KScalar a_lnpow_half     [pow_half_VALS+1] = { -1.0 };
static int     done_lnpow_half = 0;  /* set to 1 in init_math */
static KScalar a_factorial      [factorial_VALS+1] = { -1.0 };
static int     done_factorial = 0;  /* set to 1 in init_math */
static KScalar a_lnfactorial    [factorial_VALS+1] = { -1.0 };
static int     done_lnfactorial = 0;  /* set to 1 in init_math */
static KScalar a_lnbinomial     [binomial_N_VALS+1][binomial_K_VALS+1] = { -1.0 };
static int     done_lnbinomial = 0;  /* set to 1 in init_math */

/*///////////////////////////////////////////////////////////////*/
KScalar     pow_half            (KInt n)
/*
** Returns (1/2)^n for integer powers of n
*/
{
    const char* thisfunction = "pow_half";
    if (n > pow_half_VALS) {
        char buf[200];
        sprintf(buf, "%s: arg large: n=%d", thisfunction, n);
        fatal(buf);
    }
    if (! done_pow_half)
        init_math();
    return a_pow_half[n];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     lnpow_half          (KInt n)
/*
** Returns n*ln(1/2) for integer powers of n
*/
{
    const char* thisfunction = "lnpow_half";
    if (n > pow_half_VALS) {
        char buf[200];
        sprintf(buf, "%s: arg large: n=%d", thisfunction, n);
        fatal(buf);
    }
    if (! done_lnpow_half)
        init_math();
    return a_lnpow_half[n];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     factorial           (KInt n)
{
    const char* thisfunction = "factorial";
    if (n < 0 || n > factorial_VALS) {
        char buf[200];
        sprintf(buf, "%s: arg < 0 or arg > max: n=%d", thisfunction, n);
        fatal(buf);
    }
    if (! done_factorial)
        init_math();
    return a_factorial[n];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     lnfactorial         (KInt n)
/*
** Returns ln(n!).
*/
{
    const char* thisfunction = "lnfactorial";
    if (n < 0 || n > factorial_VALS) {
        char buf[200];
        sprintf(buf, "%s: arg < 0 or > max: n=%d", thisfunction, n);
        fatal(buf);
    }
    if (! done_lnfactorial)
        init_math();
    return a_lnfactorial[n];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     binomial            (KInt n, KInt k)
/*
** Returns nCk = n!/(k!(n-k)!).  This will overflow in any
** realistic parameter space.  For general use, work in
** natural log space and use lnbinomial() instead.
*/
{
    const char* thisfunction = "binomial";
    KScalar t1, t2;
    if (k < 0 || n < 0 || k > n) {
        char buf[200];
        sprintf(buf, "%s: args incorrect: n=%d, k=%d", 
                thisfunction, n, k);
        fatal(buf);
    }
    t1 = lnfactorial(n) - lnfactorial(k) - lnfactorial(n-k);
    t2 = floor(0.5 + exp(t1));
    return t2;
}


/*///////////////////////////////////////////////////////////////*/
KScalar     lnbinomial          (KInt n, KInt k)
/*
** Returns ln(nCk) = ln(n!/(k!(n-k)!)).
*/
{
    const char* thisfunction = "lnbinomial";
    if (k < 0 || n < 0 || k > n) {
        char buf[200];
        sprintf(buf, "%s: args incorrect: n=%d, k=%d", 
                thisfunction, n, k);
        fatal(buf);
    }
    if (n > binomial_N_VALS) {
        char buf[200];
        sprintf(buf, "%s: arg large: n=%d", thisfunction, n);
        fatal(buf);
    }
    if (! done_lnbinomial)
        init_math();
    k = MIN(k, n-k);  /* binomial(n,k) == binomial(n,n-k) */
    return a_lnbinomial[n][k];
}


/*///////////////////////////////////////////////////////////////*/
void        init_math           (void)
/*
** Initialize the static arrays used by the math functions.
*/
{
    const char* thisfunction = "init_math";
    {
        /* pow_half */
        KInt i;
        a_pow_half[0] = 1.0;
        a_pow_half[1] = 0.5;
        for (i=2; i <= pow_half_VALS; i++) {
            a_pow_half[i] = 0.5 * a_pow_half[i-1];
        }
        done_pow_half = 1;
    }
    {
        /* lnpow_half */
        KInt i;
        KScalar loghalf = log(0.5);
        a_lnpow_half[0] = 0.0;
        for (i=1; i <= pow_half_VALS; i++) {
            a_lnpow_half[i] = i * loghalf;
        }
        done_lnpow_half = 1;
    }
    {
        /* factorial */
        KInt i;
        a_factorial[0] = 1.0;
        /* on a Windows PC with Visual C++, this will 
        ** quietly overflow around 171 */
        for (i=1; i <= factorial_VALS; i++) {
            a_factorial[i] = (a_factorial[i-1] * i);
        }
        done_factorial = 1;
    }
    {
        /* lnfactorial */
        KInt i;
        a_lnfactorial[0] = 0.0;
        a_lnfactorial[1] = 0.0;
        for (i=2; i <= factorial_VALS; i++) {
            a_lnfactorial[i] = a_lnfactorial[i-1] + log(i);
        }
        done_lnfactorial = 1;
    }
    {
        /* lnbinomial */
        KInt nn, kk, max_kk;
        if (! done_lnfactorial) {
            char buf[200];
            sprintf(buf, "%s: lnfactorial must precede lnbinomial",
                    thisfunction);
            fatal(buf);
        }
        for (nn=0; nn <= binomial_N_VALS; nn++) {
            a_lnbinomial[nn][0] = 0.0;
            a_lnbinomial[nn][1] = log((KScalar)nn);
            /* if nn is even, then max k=n/2; 
            ** else nn is odd, then max k=(n-1)/2 
            */
            max_kk = ((nn % 2) == 0) ? nn/2 : (nn-1)/2;
            for (kk=2; kk <= max_kk; kk++) {
                a_lnbinomial[nn][kk] = lnfactorial(nn) - 
                                       lnfactorial(kk) - 
                                       lnfactorial(nn-kk);
            }
        }
        done_lnbinomial = 1;
    }
}
