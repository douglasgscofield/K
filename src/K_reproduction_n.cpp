#include "K_n.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Functions used for all forms of reproduction                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
void        set_repro_n     (KConfig_n KN, KScalar s, KScalar a)
/*
** sets up the rates assuming one genotype
*/
{
    const char* thisfunction = "set_repro_n";
    KN->S = s;
    KN->A = a;
    KN->O = 1.0 - KN->S - KN->A;
    if (KN->O < 0.0 || KN->O > 1.0) {
        char buf[200];
        sprintf(buf, "%s: rate arguments incorrect", thisfunction);
        fatal(buf);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_self_progeny_n  (KConfig_n KN)
/*
** Compute proportions of selfed genotypes produced by population.
*/
{
    const char* thisfunction = "compute_self_progeny_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_PROGENY_TO_X1 &&
        KN->current_x != KN_CURRENT_X2) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    if (KN->current_x == KN_CURRENT_X2) {
        /* we're the first reproduction function to get
        ** called, so we have to zero the destination array,
        ** which is KN->x1
        */
        fill_KArray_n(KN, KN->x1, 0.0);
    }
    apply_self_progeny_n(KN, KN->x1, KN->x2);
    KN->current_x = KN_CURRENT_PROGENY_TO_X1;
}

/*///////////////////////////////////////////////////////////////*/
void        apply_self_progeny_n    (KConfig_n KN,
                                     KArray_n& to, KArray_n& from)
/*
** Compute proportions of selfed genotypes produced by population.
*/
{
    const char* thisfunction = "apply_self_progeny_n";
    KInt i0, j0, i1, j1, n0, v0, n1, v1;
    KScalar sum, loadclass;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->S == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s: no selfing\n", thisfunction);
        return;  /* nothing to add */
    }
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (j0=0; j0 <= KN->MJ0; j0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                for (j1=0; j1 <= KN->MJ1; j1++) {
#if defined(TIMECRITICAL_INCLUDEDEBUG_n)
                    IF_DEBUG(DEBUG_TRACE2)
                        if (!(i0 % 10) && !(j0 % 10) &&
                            !(i1 % 10) && !(j1 % 10)) {
                            fprintf(stderr, "s[%d,%d,%d,%d] ", 
                                   i0, j0, i1, j1);
                            fflush(stderr);
                        }
#endif
                    /* L(i0,j0,i1,j1) is the destination load class */
                    if ((KN->is_lethal[0] && j0 > 0 && KN->createlethal[0] == 0) ||
                        (KN->is_lethal[1] && j1 > 0 && KN->createlethal[1] == 0)) {
                        /* don't do any more to a class that will die */
                        continue;
                    }
                    /*
                    ** See the note in apply_self_progeny() regarding the
                    ** temptation to simplify the loop control here.
                    */
                    sum = 0.0;
                    for (n0=i0; n0 <= KN->MI0; n0++) {
                        for (v0=0; v0 <= j0; v0++) {
                            if (v0 + n0 < j0 || n0 + v0 < i0 + j0)
                                /* because s_self_n() is only valid with
                                ** j0 <= v0+n0, skip all v0+n0 < j0.
                                ** because s_self_n() is only valid with
                                ** i0+j0 <= n0+v0, skip all n0, v0 until
                                ** this is true.
                                */
                                continue;
                            for (n1=i1; n1 <= KN->MI1; n1++) {
                                for (v1=0; v1 <= j1; v1++) {
                                    if (v1 + n1 < j1 || n1 + v1 < i1 + j1)
                                        /* because s_self_n() is only valid with
                                        ** j1 <= v1+n1, skip all v1+n1 < j1.
                                        ** because s_self_n() is only valid with
                                        ** i1+j1 <= n1+v1, skip all n1, v1 until
                                        ** this is true.
                                        */
                                        continue;
                                    if ((loadclass=from[n0][v0][n1][v1]) == 0.0) 
                                        continue;
                                    /* 
                                    ** These checks all must be redundant.
                                    **
                                    ** if (i0 <= n0 && j0 >= v0 && 
                                    **     j0 <= n0+v0 && i0+j0 <= n0+v0 &&
                                    **     i1 <= n1 && j1 >= v1 && 
                                    **     j1 <= n1+v1 && i1+j1 <= n1+v1) {
                                    **     sum += s_self_n(i0,j0,i1,j1,n0,v0,n1,v1) * 
                                    **            KN->S *
                                    **            loadclass;
                                    ** }
                                    */
                                    sum += s_self_n(i0,j0,i1,j1,n0,v0,n1,v1) * 
                                           KN->S *
                                           loadclass;
                                }
                            }
                        }
                    }
                    /* final assignment to x"(s)i0,j0,i1,j1 */
                    to[i0][j0][i1][j1] += sum;
                }
            }
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_apomixis_progeny_n  (KConfig_n KN)
/*
** Compute proportions of apomixis genotypes produced by population.
*/
{
    const char* thisfunction = "compute_apomixis_progeny_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_PROGENY_TO_X1 &&
        KN->current_x != KN_CURRENT_X2) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    if (KN->current_x == KN_CURRENT_X2) {
        /* we're the first reproduction function to get
        ** called, so we have to zero the destination array,
        ** which is KN->x1
        */
        fill_KArray_n(KN, KN->x1, 0.0);
    }
    apply_apomixis_progeny_n(KN, KN->x1, KN->x2);
    KN->current_x = KN_CURRENT_PROGENY_TO_X1;
}
    
/*///////////////////////////////////////////////////////////////*/
void        apply_apomixis_progeny_n    (KConfig_n KN,
                                         KArray_n& to, KArray_n& from)
/*
** Compute proportions of apomixis genotypes produced by population.
*/
{
    const char* thisfunction = "apply_apomixis_progeny_n";
    KInt i0, j0, i1, j1, n0, v0, n1, v1;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->A == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s: no apomixis\n", thisfunction);
        return;  /* nothing to add */
    }
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (j0=0; j0 <= KN->MJ0; j0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                for (j1=0; j1 <= KN->MJ1; j1++) {
                    if ((KN->is_lethal[0] && j0 > 0
                         && KN->createlethal[0] == 0)
                        || (KN->is_lethal[1] && j1 > 0
                            && KN->createlethal[1] == 0)) {
                        // don't do any more to a class that will die **
                        continue;
                    }
                    IF_DEBUG(DEBUG_TRACE2)
                        if (!(i0 % 10) && !(j0 % 10) &&
                            !(i1 % 10) && !(j1 % 10)) {
                            fprintf(stderr, "a[%d,%d,%d,%d] ", 
                                   i0, j0, i1, j1);
                            fflush(stderr);
                        }
                    /* identical load class only */
                    n0 = i0;
                    v0 = j0;
                    n1 = i1;
                    v1 = j1;
                    /* final addition to to(a)i0,j0,i1,j1 */
                    to[i0][j0][i1][j1] += 
                        a_apomixis_n(i0,j0,i1,j1,n0,v0,n1,v1) *
                        KN->A *
                        from[n0][v0][n1][v1];
                }
            }
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_gametes_n   (KConfig_n KN)
/*
** Compute genotype proportions for male and female gametes.
*/
{
    const char* thisfunction = "compute_gametes_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_PROGENY_TO_X1 &&
        KN->current_x != KN_CURRENT_X2) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    if (KN->current_x == KN_CURRENT_X2) {
        /* we're the first reproduction function to get
        ** called, so we have to zero the destination array,
        ** which is KN->x1
        */
        fill_KArray_n(KN, KN->x1, 0.0);
    }
    apply_gametes_n(KN, KN->mgam, KN->fgam, KN->x2);
    KN->current_x = KN_CURRENT_GAMETES_X2_TO_X1;
}


/*///////////////////////////////////////////////////////////////*/
void        apply_gametes_n     (KConfig_n KN, 
                                 KVector_n& mgam, KVector_n& fgam,
                                 KArray_n& from)
/*
** Compute genotype proportions for male and female gametes.
*/
{
    const char* thisfunction = "apply_gametes_n";
    apply_gametes_full_n(KN, mgam, fgam, from);
    adjust_gametes_n(KN, mgam, fgam, from);
}


/*///////////////////////////////////////////////////////////////*/
void        apply_gametes_full_n    (KConfig_n KN, 
                                     KVector_n& mgam, KVector_n& fgam,
                                     KArray_n& from)
/*
** Compute genotype proportions for male and female gametes.
*/
{
    const char* thisfunction = "apply_gametes_full_n";
    KInt i0, i1, n0, n1, v0, v1;
    KScalar sum;
    KScalar t1, t2, t3, t4;
    static KScalar p[MAX_MI_n+1][MAX_MI_n+1][MAX_MI_n+1][MAX_MI_n+1];
    static int runonce = 0;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (!runonce) {
        runonce++;
        for (v1=0; v1 <= KN->MI1; v1++) {
            for (v0=0; v0 <= KN->MI0; v0++) {
                for (n1=v1; n1 <= KN->MI1; n1++) {
                    for (n0=v0; n0 <= KN->MI0; n0++) {
                        t1 = lnbinomial(n0, v0);
                        t2 = lnbinomial(n1, v1);
                        t3 = lnpow_half(n0 + n1);
                        t4 = exp(t1 + t2 + t3);
                        p[n0][v0][n1][v1] = t4;
                    }
                }
            }
        }
    }
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (i1=0; i1 <= KN->MI1; i1++) {
            IF_DEBUG(DEBUG_TRACE2)
                if (!(i0 % 10) && !(i1 % 10)) {
                    fprintf(stderr, "g[%d,%d] ", i0, i1);
                    fflush(stderr);
                }
            sum = 0.0;
            for (v0=0; v0 <= i0 && v0 <= KN->MJ0; v0++) {
                if (KN->is_lethal[0] && v0 > 0 && KN->createlethal[0] == 0)
                    /* When the mutation class is lethal, it's not
                    ** possible to have homozygous mutant adults
                    ** that contribute to the gamete pools
                    */
                    break;  /* this ends the for (v0=0;;) loop */
                for (v1=0; v1 <= i1 && v1 <= KN->MJ1; v1++) { 
                    if (KN->is_lethal[1] && v1 > 0 && KN->createlethal[1] == 0)
                        /* see above */
                        break;  /* this ends the for (v1=0;;) loop */
                    for (n0=(i0 - v0); n0 <= KN->MI0; n0++) {
                        for (n1=(i1 - v1); n1 <= KN->MI1; n1++) {
                            if (from[n0][v0][n1][v1] == 0.0)
                                continue;
                            sum += p[n0][i0-v0][n1][i1-v1] * 
                                   from[n0][v0][n1][v1];
                        }
                    }
                }
            }
            /* NO SPECIAL PROCESSING OF GAMETES SHOULD OCCUR HERE.
            ** THAT MUST ALL BE DONE IN adjust_gametes_n()
            */
            mgam[i0][i1] = sum;
            fgam[i0][i1] = sum;
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        adjust_gametes_n    (KConfig_n KN, 
                                 KVector_n& mgam, KVector_n& fgam,
                                 KArray_n& from)
/*
** Compute genotype proportions for male and female gametes.
*/
{
    const char* thisfunction = "adjust_gametes_n";
    KInt i0, i1;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (i1=0; i1 <= KN->MI1; i1++) {
            /* mgam[i0][i1] = mgam[i0][i1]; */
            fgam[i0][i1] *= KN->O;
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_zygotes_n   (KConfig_n KN)
/*
** Compute proportion of outcrossed zygotes produced by population.
*/
{
    const char* thisfunction = "compute_zygotes_n";
    /* Combines proportions of male and female gametes into
    ** expected genotypic classes.  Note that all classes 
    ** L(,j>0) = 0.
    */
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_GAMETES_X2_TO_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    apply_zygotes_n(KN, KN->x1, KN->mgam, KN->fgam);
    KN->current_x = KN_CURRENT_PROGENY_TO_X1;
}

/*///////////////////////////////////////////////////////////////*/
void        apply_zygotes_n   (KConfig_n KN, KArray_n& to,
                               KVector_n& mgam, KVector_n& fgam)
/*
** Compute proportion of outcrossed zygotes produced by population.
*/
{
    const char* thisfunction = "apply_zygotes_n";
    KInt i0, i1, k0, k1;
    KScalar sum;    
    /* Combines proportions of male and female gametes into
    ** expected genotypic classes.  Note that all classes 
    ** L(,j>0) = 0.
    */
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->O == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s: no outcrossing\n",
                                      thisfunction);
        return;
    }
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (i1=0; i1 <= KN->MI1; i1++) {
            IF_DEBUG(DEBUG_TRACE2)
                if (!(i0 % 10) && !(i1 % 10)) {
                    fprintf(stderr, "z[%d,0,%d,0] ", i0, i1);
                    fflush(stderr);
                }
            sum = 0.0;
            for (k0=0; k0 <= i0; k0++) {
                for (k1=0; k1 <= i1; k1++) {
                    sum += ((fgam[k0][k1] * mgam[i0-k0][i1-k1]) +
                            (fgam[k0][i1-k1] * mgam[i0-k0][k1]) +
                            (fgam[i0-k0][k1] * mgam[k0][i1-k1]) +
                            (fgam[i0-k0][i1-k1] * mgam[k0][k1]));
                }
            }
            /* Gamete pools can each sum to one (selfing and apomixis
            ** decrease the sum of K->fgam) so when we combine them
            ** into the zygote pool, we have to divide by 4 to get the
            ** proportions of zygotes.
            */
            /* note this adds the proportions to the 'to' array */
            to[i0][0][i1][0] += sum * 0.25;
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_summed_progeny_n    (KConfig_n KN)
{
    const char* thisfunction = "compute_summed_progeny_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_PROGENY_TO_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    apply_summed_progeny_n(KN, KN->x1, KN->x2);
    KN->current_x = KN_CURRENT_X1;
}

/*///////////////////////////////////////////////////////////////*/
void        apply_summed_progeny_n  (KConfig_n KN,
                                     KArray_n& to, KArray_n& from)
{
    const char* thisfunction = "apply_summed_progeny_n";
    /* nothing to do, all this has been done already by each
    ** reproduction operation, which added new progeny 
    ** proportions to KN->x1
    */
    check_normalization_n(KN, to, thisfunction, "to");
}

/*///////////////////////////////////////////////////////////////*/
KScalar     s_self_n        (KInt i0, KInt j0, KInt i1, KInt j1, 
                             KInt n0, KInt v0, KInt n1, KInt v1)
/*
** Compute selfing proportion produced by given class.
** Proportion of L(i,j) produced when L(n,v) is selfed.
** The proportion of L(i,j) is returned.
** No need to know population data, works solely from i, j, n, v.
*/
{
    const char* thisfunction = "s_self_n";
    KScalar ans;
    KScalar t1, t2, t3, t4, t5;
    /* check for classes that cannot be produced through selfing */
    if (i0 > n0 || j0 < v0 || j0 > n0+v0 || i0+j0 > n0+v0 ||
        i1 > n1 || j1 < v1 || j1 > n1+v1 || i1+j1 > n1+v1)
        return 0.0;
    /* if the class can be produced, then compute its expected proportion */
    /* log scale */
    t1 = lnbinomial(n0, i0);
    t2 = lnbinomial((n0 - i0), (j0 - v0));
    t3 = lnbinomial(n1, i1);
    t4 = lnbinomial((n1 - i1), (j1 - v1));
    t5 = lnpow_half(2*(n0 + n1) - i0 - i1);
    ans = exp(t1 + t2 + t3 + t4 + t5);
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     a_apomixis_n    (KInt i0, KInt j0, KInt i1, KInt j1, 
                             KInt n0, KInt v0, KInt n1, KInt v1)
/*
** Compute apomictic proportion produced by given class.
** Proportion of L(i0,j0, i1, j1) produced when 
** L(n0,v0,n1,v1) undergoes apomixis.
*/
{
    const char* thisfunction = "a_apomixis_n";
    if (i0 == n0 && j0 == v0 && i1 == n1 && j1 == v1)
        return 1.0;
    else
        return 0.0;
}


