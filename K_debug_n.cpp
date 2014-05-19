#include "K_n.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Utility routines for debugging - for nested mutation classes  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
void        check_normalization_n   (KConfig_n KN, KArray_n& a,
                                     const char* caller, 
                                     char* array_name)
{
    const char* thisfunction = "check_normalization_n";
    IF_DEBUG(DEBUG_NORMALIZATION) {
        if (!isOK_KArray_n(KN, a)) {
            fprintf(stderr, "%s: %s not normalized: %f, tolerance=%lg\n", 
                   caller, array_name, sum_KArray_n(KN, a),
                   NORMALIZATION_TOLERANCE);
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        dump_KArray_n           (KConfig_n KN, KArray_n& a,
                                     KInt mi0, KInt mj0,
                                     KInt mi1, KInt mj1)
{
    const char* thisfunction = "dump_KArray_n";
    if (mi0 > KN->MI0 || mi0 < 0)
        mi0 = KN->MI0;
    if (mj0 > KN->MJ0 || mj0 < 0)
        mj0 = KN->MJ0;
    if (mi1 > KN->MI1 || mi1 < 0)
        mi1 = KN->MI1;
    if (mj1 > KN->MJ1 || mj1 < 0)
        mj1 = KN->MJ1;
    fprintf(stderr, "%s BEGIN KArray_n @ 0x%08x; max i0=%d, j0=%d, i1=%d, j1=%d\n", 
           thisfunction, a, mi0, mj0, mi1, mj1);
    fprintf(stderr, "i0\tj0\ti1\tj1\ta[i0][j0][i1][j1]\n");
    dump_values_KArray_n(KN, stdout, a, mi0, mj0, mi1, mj1);
    fprintf(stderr, "%s END for KArray_n a @ 0x%8x\n", thisfunction, a);
}

/*///////////////////////////////////////////////////////////////*/
void        dump_values_KArray_n    (KConfig_n KN,
                                     FILE* fp,
                                     KArray_n& a,
                                     KInt mi0, KInt mj0,
                                     KInt mi1, KInt mj1)
{
    KInt i0, j0, i1, j1;
    for (i0=0; i0 <= mi0; i0++) {
        for (j0=0; j0 <= mj0; j0++) {
            for (i1=0; i1 <= mi1; i1++) {
                for (j1=0; j1 <= mj1; j1++) {
                    fprintf(fp, "%d\t%d\t%d\t%d\t%lg\n", 
                            i0, j0, i1, j1, a[i0][j0][i1][j1]);
                }
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        dump_KArray_full_n      (KConfig_n KN, KArray_n& a)
{
    const char* thisfunction = "dump_KArray_full_n";
    fprintf(stderr, "%s BEGIN KArray_n @ 0x%08x; dump_KArray_n follows\n", 
           thisfunction, a);
    dump_KArray_n(KN, a, KN->MI0, KN->MJ0, KN->MI1, KN->MJ1);
    fprintf(stderr, "%s END for KArray_n a @ 0x%8x\n", thisfunction, a);
}

/*///////////////////////////////////////////////////////////////*/
void        dump_KVector_n          (KConfig_n KN, KVector_n& v,
                                     KInt mi0, KInt mi1)
{
    const char* thisfunction = "dump_KVector_n";
    KInt i0, i1;
    if (mi0 > KN->MI0 || mi0 < 0)
        mi0 = KN->MI0;
    if (mi1 > KN->MI1 || mi1 < 0)
        mi1 = KN->MI1;
    fprintf(stderr, "%s BEGIN KVector_n @ 0x%08x; max i0=%d, i1=%d\n", 
           thisfunction, v, mi0, mi1);
    fprintf(stderr, "i0\ti1\tv[i0][i1]\n");
    for (i0=0; i0 <= mi0; i0++) {
        for (i1=0; i1 < mi1; i1++) {
                fprintf(stderr, "%d\t%d\t%lg\n", i0, i1, v[i0][i1]);
        }
    }
    fprintf(stderr, "%s END for KVector_n v @ 0x%8x\n", thisfunction, v);
}

/*///////////////////////////////////////////////////////////////*/
void        dump_KVector_full_n     (KConfig_n KN, KVector_n& v)
{
    const char* thisfunction = "dump_KVector_full_n";
    fprintf(stderr, "%s BEGIN KVector_n @ 0x%08x; dump_KVector_n follows\n", 
           thisfunction, v);
    dump_KVector_n(KN, v, KN->MI0, KN->MI1);
    fprintf(stderr, "%s END for KVector_n v @ 0x%8x\n", thisfunction, v);
}

