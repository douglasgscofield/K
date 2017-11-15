#include "K.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Utility routines for debugging                                */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
KScalar     sum_KArray          (KConfig K, KArray& A)
/*
** Computes sum of all array elements in A according to K
*/
{
    //const char* thisfunction = "sum_KArray";
    KInt i, j, g;
    KScalar ans = 0.0;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                ans += A[i][j][g];
            }
        }
    }
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     sum_KArray2         (KConfig K, KArray2& A2)
/*
** Computes sum of all array elements in A according to K
*/
{
    //const char* thisfunction = "sum_KArray2";
    KInt i, j;
    KScalar ans = 0.0;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            ans += A2[i][j];
        }
    }
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     sum_KVector1        (KConfig K, KVector1& v)
/*
** Computes sum of all vector elements in v according to K
*/
{
    //const char* thisfunction = "sum_KVector1";
    KInt i;
    KScalar ans = 0.0;
    for (i=0; i <= K->MI; i++) {
        ans += v[i];
    }
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
void        normalize_KArray    (KConfig K, KArray& a)
/*
** Normalizes A according to K so that it sums to 1.  
*/
{
    const char* thisfunction = "normalize_KArray";
    KInt i, j, g;
    KScalar sum = sum_KArray(K, a);
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_NORMALIZATION) 
        fprintf(stderr, "%s: prenormalization, sum_KArray = %lg\n", 
               thisfunction, sum);
    if (sum == 1.0)
        return;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                a[i][j][g] /= sum;
            }
        }
    }
    IF_DEBUG(DEBUG_NORMALIZATION) 
        fprintf(stderr, "%s: postnormalization, sum_KArray = %lg\n", 
               thisfunction, sum_KArray(K, a));
}

/*///////////////////////////////////////////////////////////////*/
void        truncate_KArray     (KConfig K, KArray& a, KScalar v)
/*
** Truncate a so that any load class-genotype that 
** is < v becomes 0.0
*/
{
    const char* thisfunction = "truncate_KArray";
    KInt i, j, g;
    KInt num = 0;
    KScalar sum = 0.0;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_TRUNCATE_DETAIL) {
        fprintf(stderr, "%s: truncating values lower than v=%lg \n", 
               thisfunction, v);
    }
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                if (a[i][j][g] > 0.0 && a[i][j][g] < v) {
                    sum += v - a[i][j][g];
                    num++;
                    a[i][j][g] = 0.0;
                    IF_DEBUG(DEBUG_TRUNCATE_DETAIL) {
                        fprintf(stderr, "%s: truncated a[%d][%d][%d]\n",
                               thisfunction, i, j, g);
                    }
                }
            }
        }
    }
    IF_DEBUG(DEBUG_TRUNCATE) {
        if (num > 0) {
            fprintf(stderr, "%s: truncated %d elements, sum %lg\n", thisfunction, num, sum);
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        fill_KArray         (KConfig K, KArray& a, KScalar val)
{
    //const char* thisfunction = "fill_KArray";
    KInt i, j, g;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                a[i][j][g] = val;
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        copy_KArray         (KConfig K, KArray& to, KArray& from)
{
    //const char* thisfunction = "copy_KArray";
    KInt i, j, g;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                to[i][j][g] = from[i][j][g];
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
int         isOK_KArray         (KConfig K, KArray& a)
/*
** checks to see if the classes in a[][] sum to 1
*/
{
    //const char* thisfunction = "isOK_KArray";
    KScalar t1;
    t1 = fabs(sum_KArray(K, a) - 1.0);
    return (t1 < NORMALIZATION_TOLERANCE) ? 1 : 0;
}

/*///////////////////////////////////////////////////////////////*/
int         isOK_KVector1       (KConfig K, KVector1& v)
/*
** checks to see if the classes in v[] sum to 1
*/
{
    //const char* thisfunction = "isOK_KVector1";
    KScalar t1;
    t1 = fabs(sum_KVector1(K, v) - 1.0);
    return (t1 < NORMALIZATION_TOLERANCE) ? 1 : 0;
}


/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Utility routines for error handling                           */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
void        fatal               (const char* msg)
{
    //const char* thisfunction = "fatal";
    fprintf(stderr, "FATAL ERROR: %s \n", msg);
    exit(1);
}

/*///////////////////////////////////////////////////////////////*/
void        warning             (const char* msg)
{
    //const char* thisfunction = "warning";
    fprintf(stderr, "WARNING: %s \n", msg);
}


