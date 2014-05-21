#include "K_n.h"


/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Functions for determining equilibrium - nested mutations      */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/


/*///////////////////////////////////////////////////////////////*/
int         is_equilibrium_n    (KConfig_n KN)
/* 
** return 1 if KN satisfies equilibrium constraints as defined by user
*/
{
    const char* thisfunction = "is_equilibrium_n";
    int ans;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    ans = apply_is_equilibrium_n(KN, KN->x1, KN->x_prevgen);
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
int         apply_is_equilibrium_n  (KConfig_n KN,
                                     KArray_n& adults,
                                     KArray_n& prevgen)
/* 
** return 1 if KN satisfies equilibrium constraints as defined by user
*/
{
    const char* thisfunction = "apply_is_equilibrium_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (diff_epsilon_KArray_n(KN, adults, prevgen, KN->epsilon)) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s: no equilibrium\n", thisfunction);
        return 0;
    } else {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s: at equilibrium\n", thisfunction);
        return 1;
    }
}

/*///////////////////////////////////////////////////////////////*/
int         diff_epsilon_KArray_n   (KConfig_n KN,
                                     KArray_n& a1, KArray_n& a2,
                                     KScalar epsilon)
/*
** Checks to see if any of the proportions in a1 and a2 differ by
** more than epsilon.  If true, returns 1, else 0.  
*/
{
    const char* thisfunction = "diff_epsilon_KArray_n";
    KInt i0, j0, i1, j1;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_EQUILIBRIUM) {
        KInt n0, v0, n1, v1;
        KInt l0 = 0, lam0 = 0, l1 = 0, lam1 = 0;
        KScalar diff, maxdiff = 0.0;
        IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
            fprintf(stderr, "%s: BEGIN where (a1[][][][]-a2[][][][] != 0.0)\n",
                   thisfunction);
            fprintf(stderr, "i0\tj0\ti1\tj1\tdiff\n");
        }
        for (n0=0; n0 < KN->MI0; n0++) {
            for (v0=0; v0 < KN->MJ0; v0++) {
                for (n1=0; n1 < KN->MI1; n1++) {
                    for (v1=0; v1 < KN->MJ1; v1++) {
                        diff = fabs(a1[n0][v0][n1][v1] - a2[n0][v0][n1][v1]);
                        if (diff > maxdiff) {
                            maxdiff = diff; 
                            l0 = n0; lam0 = v0; l1 = n1; lam1 = v1;
                        }
                        IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
                            if (n0 < 10 && v0 < 10 &&
                                n1 < 10 && v1 < 10) {
                                diff = a1[n0][v0][n1][v1] - a2[n0][v0][n1][v1];
                                if (diff != 0.0) {
                                    fprintf(stderr, "%d\t%d\t%d\t%d\t%lg\n", 
                                           n0, v0, n1, v1, diff);
                                }
                            }
                        }
                    }
                }
            }
        }
        IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
            fprintf(stderr, "%s: END (a1[][][][]-a2[][][][] != 0.0)\n", thisfunction);
        }
        fprintf(stderr, "%s: epsilon = %lg, max diff @ [%d][%d][%d][%d] = %lg\n",
               thisfunction, epsilon, l0, lam0, l1, lam1, maxdiff);
    }
    /* TODO: sort out whether these should be "<=" or "<" */
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (j0=0; j0 <= KN->MJ0; j0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                for (j1=0; j1 <= KN->MJ1; j1++) {
                    if (fabs(a1[i0][j0][i1][j1] - a2[i0][j0][i1][j1]) > epsilon) {
                        /* only need one element that exceeds epsilon */
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}


