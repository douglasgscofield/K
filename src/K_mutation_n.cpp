#include "K_n.h"


/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for computing and applying mutation                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/* NOTE: The routines for Poisson terms are in K_math.c */

/*///////////////////////////////////////////////////////////////*/
void        compute_mutation_n  (KConfig_n KN)
/* 
** Apply mutations to each mutation class for adult load classes
**     held in KN->x1 with the new load classes in KN->x2.
** If you want to add mutations to an arbitrary KArray_n, then use 
**     the function 
**     apply_mutation_n(KConfig_n KN, KArray_n& to, KArray_n& from)
** KN->x2 is correct upon exit.
*/
{
    const char* thisfunction = "compute_mutation_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    apply_mutation_n(KN, KN->x2, KN->x1);
    KN->current_x = KN_CURRENT_X2;
}

/*///////////////////////////////////////////////////////////////*/
void        apply_mutation_n    (KConfig_n KN, KArray_n& to,
                                 KArray_n& from)
{
    const char* thisfunction = "apply_mutation_n";
    KInt i0, j0, i1, j1, n0, n1;
    KScalar sum;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->U[0] == 0.0 && KN->U[1] == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "no mutation\n");
        copy_KArray_n(KN, to, from);
        check_normalization_n(KN, to, thisfunction, "to");
        return;
    }
    for (j1=0; j1 <= KN->MJ1; j1++) {
        for (i1=0; i1 <= KN->MI1; i1++) {
            for (j0=0; j0 <= KN->MJ0; j0++) {
                for (i0=0; i0 <= KN->MI0; i0++) {
                    if ((KN->is_lethal[0] && j0 > 0 && KN->createlethal[0] == 0) ||
                        (KN->is_lethal[1] && j1 > 0 && KN->createlethal[1] == 0)) {
                        /* When the mutation class has homozygous lethals, 
                        ** it makes no sense to add heterozygous mutations to it.
                        */
                        to[i0][j0][i1][j1] = 0.0;
                        break;  /* this ends the for (i0=0;;) loop */
                    }
#if defined(TIMECRITICAL_INCLUDEDEBUG_n)
                    IF_DEBUG(DEBUG_TRACE2)
                        if (!(i0 % 10) && !(j0 % 10) &&
                            !(i1 % 10) && !(j1 % 10))
                            fprintf(stderr, "mut[%d,%d,%d,%d] ", 
                                   i0, j0, i1, j1);
#endif
                    sum = 0.0;
                    for (n1=0; n1 <= i1; n1++) {
                        for (n0=0; n0 <= i0; n0++) {
                            sum += from[n0][j0][n1][j1] * 
#if defined(MUTTERM_FUNCTION_n)
                                   mut_term_n(KN, i0-n0, i1-n1);
#else
                                   KN->mut_term[i0-n0][i1-n1];
#endif
                        }
                    }
                    to[i0][j0][i1][j1] = sum;
                }
            }
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
    check_normalization_n(KN, to, thisfunction, "to");
}

/*///////////////////////////////////////////////////////////////*/
/* TODO: IMPLEMENT THIS */
/*///////////////////////////////////////////////////////////////*/
void        apply_mutation_general_n    (KConfig K, KScalar U, 
                                         KArray& to, KArray& from)
{
    const char* thisfunction = "apply_mutation_general";
    KInt g, i, j, n1;
    KScalar sum;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->U == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "no mutation\n");
        copy_KArray(K, K->xp, K->x);
        return;
    }
    for (g=0; g < K->genotypes; g++) {
        for (j=0; j <= K->MJ; j++) {
            for (i=0; i <= K->MI; i++) {
                sum = 0.0;
                for (n1=0; n1 <= i; n1++) {
                    sum += from[n1][j][g] * 
                           mut_term_general(U, i - n1);
                }
                to[i][j][g] = sum;
            }
        }
    }
    check_normalization(K, to, thisfunction, "to");
}

/*///////////////////////////////////////////////////////////////*/
void        initiate_mut_term_n     (KConfig_n KN)
{
    const char* thisfunction = "initiate_mut_term_n";
    KInt i0, i1;
    KScalar t0, t1;
    KScalar checksum = 0.0;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    /*
    ** pois_term() behaves reasonably with U[0]==0.0 or U[1]==0.0
    **
    if (KN->U[0] == 0.0 || KN->U[1] == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "no mutation\n");
        not_implemented(thisfunction, "U[0]==0.0 or U[1]==0.0");
        for (i1=0; i1 <= KN->MI1; i1++) {
            for (i0=0; i0 <= KN->MI0; i0++) {
                KN->mut_term[i0][i1] = 0.0;
            }
        }
        return;
    }
    **/
    for (i0=0; i0 <= KN->MI0; i0++) {
        t0 = pois_term(KN->U[0], i0);
        for (i1=0; i1 <= KN->MI1; i1++) {
            t1 = pois_term(KN->U[1], i1) * t0;
            KN->mut_term[i0][i1] = t0;
            checksum += t1;
        }
    }
    /* verify that the contents of the mut_term[] array are OK */
    if (fabs(checksum - 1.0) > 0.000000001) {
        char buf[200];
        sprintf(buf, "%s: not valid pdf, 1 != %lf\n", 
                thisfunction, checksum);
        warning(buf);
    }
    IF_DEBUG(DEBUG_NORMALIZATION) {
        KInt i0, i1;
        KScalar termsum = 0.0;
        not_implemented(thisfunction, "debug normalization");
        for (i0=0; i0 <= KN->MI0; i0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                termsum += KN->mut_term[i0][i1];
            }
        }
        fprintf(stderr, "%s: U[0]=%lg, U[1]=%lg, sum mut_term=%lf\n", 
               thisfunction, KN->U[0], KN->U[1], termsum);
    }
}

/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term_n      (KConfig_n KN, KInt x0, KInt x1)
{
    const char* thisfunction = "mut_term_n";
    if (x0 > KN->MI0 || x1 > KN->MI1) {
        char buf[200];
        sprintf(buf, "%s: term too large", thisfunction);
        fatal(buf);
    }
    return KN->mut_term[x0][x1];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term_computed_n (KConfig_n KN, KInt x0, KInt x1)
{
    const char* thisfunction = "mut_term_computed_n";
    if (x0 > KN->MI0 || x1 > KN->MI1) {
        char buf[200];
        sprintf(buf, "%s: term too large", thisfunction);
        fatal(buf);
    }
    return pois_term(KN->U[0], x0) * pois_term(KN->U[1], x1);
}

/*///////////////////////////////////////////////////////////////*/
/* TODO: IMPLEMENT THIS */
/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term_general_n      (KScalar U, KInt x)
/*
** Returns the xth term of the Poisson distribution with mean U.
** Keeps internal arrays of terms generated for each U seen.  
**     The first call is slow, but all subsequent calls with the
**     same U simply require an array lookup.
*/
{
    const char* thisfunction = "mut_term_general";
    /* note that the term [][MAX_MI+1] of pois_general 
    ** holds the U value
    */
    static KScalar  pois_general    [pois_general_MAX][MAX_MI+2];
    static KInt     last_index      = -1;
    static KScalar  prev_U          = -1.0;
    static KInt     prev_index      = -1;
    not_implemented(thisfunction, "generalization of mutation");
    if (x > MAX_MI) {
        char buf[200];
        sprintf(buf, "%s: term too large", thisfunction);
        fatal(buf);
    }
    if (U == prev_U) {
        return pois_general[prev_index][x];
    } else {
        KInt i;
        KScalar t1;
        KScalar checksum = 0.0;
        /* find an existing entry in the pois_general array */
        for (i=0; i <= last_index; i++) {
            if (pois_general[i][MAX_MI+1] == U) {
                prev_U = U;
                prev_index = i;
                return pois_general[prev_index][x];
            }
        }
        /* entry not found; also, i == last_index+1 */
        if (i >= pois_general_MAX) {
            char buf[200];
            sprintf(buf, "%s: pois_general_MAX too low",
                    thisfunction);
            fatal(buf);
        }
        /* create new entry */
        last_index = i;
        for (i=0; i <= MAX_MI; i++) {
            t1 = pois_term(U, i);
            pois_general[last_index][i] = t1;
            checksum += t1;
        }
        /* verify that the contents of the mut_term[] array 
        ** are OK
        */
        if (fabs(checksum - 1.0) > 0.000000001) {
            char buf[200];
            sprintf(buf, "%s: not valid pdf, U=%lg: 1 != %f\n", 
                    thisfunction, U, checksum);
            warning(buf);
        }
        IF_DEBUG(DEBUG_NORMALIZATION) {
            KInt i;
            KScalar termsum = 0.0;
            fprintf(stderr, "%s: U = %lg\n", thisfunction, U);
            for (i=0; i <= MAX_MI; i++) {
                termsum += pois_general[last_index][i];
                fprintf(stderr, "pois_general[last_index][%d] = %lg, cumulative sum = %lg\n",
                       i, pois_general[last_index][i], termsum);
            }
        }
        pois_general[last_index][MAX_MI+1] = U;
        prev_U = U;
        prev_index = last_index;
        return pois_general[prev_index][x];
    }
}

