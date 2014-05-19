#include "K.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Functions used for all forms of reproduction                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
void        set_repro               (KConfig K, 
                                     KScalar s, KScalar ds,
                                     KScalar a, KScalar da)
/*
** sets up the rates assuming one genotype
*/
{
    const char* thisfunction = "set_repro";
    KInt g = 0;
    if (K->genotypes != 1) {
        char buf[200];
        sprintf(buf, "%s: more than one genotype seen, >1 unimplemented",
                thisfunction);
        fatal(buf);
    }
    K->S[g] = s;
    K->D_S[g] = ds;
    K->A[g] = a;
    K->D_A[g] = da;
    K->O[g] = 1.0 - K->S[g] - K->A[g];
    if (K->O[g] < 0.0 || K->O[g] > 1.0) {
        char buf[200];
        sprintf(buf, "%s: rate arguments incorrect", thisfunction);
        fatal(buf);
    }
    set_repro_resources(K);
    set_transform_one_genotype(K);
}

/*///////////////////////////////////////////////////////////////*/
void        set_repro_resources     (KConfig K)
/*
** sets up the arrays that drive the reproductive resource algorithms
*/
{
    //const char* thisfunction = "set_repro_resources_kondrashov";
    KInt g;
    for (g=0; g < K->genotypes; g++) {
        K->rsrc_SO[g] = 0.5 * (K->S[g] + K->D_S[g]*K->S[g]);
        K->rsrc_AO[g] = 0.5 * (K->A[g] + K->D_A[g]*K->A[g]);
        K->rsrc_OO[g] = 0.5 * K->O[g];
        K->rsrc_OP[g] = 0.5 * (1 - K->D_S[g]*K->S[g] - K->D_A[g]*K->A[g]);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        set_repro_resources_kondrashov  (KConfig K)
/*
** sets up the arrays that drive the reproductive resource algorithms
*/
{
    //const char* thisfunction = "set_repro_resources_kondrashov";
    KInt g;
    for (g=0; g < K->genotypes; g++) {
        K->rsrc_SO[g] = 0.5 * (K->S[g] + K->D_S[g]*K->S[g]);
        K->rsrc_AO[g] = 0.5 * (K->A[g] + K->D_A[g]*K->A[g]);
        K->rsrc_OO[g] = 0.5 * K->O[g];
        K->rsrc_OP[g] = 0.5 * (1 - K->D_S[g]*K->S[g] - K->D_A[g]*K->A[g]);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_self_progeny    (KConfig K)
{
    const char* thisfunction = "compute_self_progeny";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_self_progeny(K, K->xpps, K->xp);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_self_progeny  (KConfig K, KArray& to, KArray& from)
/*
** Compute proportions of selfed genotypes produced by population.
** This follows Charlesworth et al. 1990.
** K->xpps is correct on exit
*/
{
    const char* thisfunction = "apply_self_progeny";
    KInt i, j, n, v, g, xi;
    KScalar sum, loadclass;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    for (g=0; g < K->genotypes; g++) {
        /* g is the destination genotype */
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                if (K->is_lethal && j > 0 && K->createlethal == 0) {
                    // don't do any more to a class that will die **
                    to[i][j][g] = 0.0;
                    continue;
                }
                /* L(i,j) is the destination load class */
                sum = 0.0;
                xi = 0;
                if (K->S[xi] == 0.0) {
                    to[i][j][g] = 0.0;
                    continue;
                }
                /*
                ** Note: it is tempting to optimize the for(v=0;;) loop
                ** below, by starting at v=i+j-n, and then removing the
                ** if(... || ...) within the loop body.  Don't do it!  And,
                ** the if() conditionals are redundant, in that the left
                ** expression is always true if the right expression is true,
                ** but leave them both!
                **
                ** You'll get non-zero mean self fitness during periods
                ** of selective interference, when you should be getting
                ** zero or very-near-zero self fitness.
                **
                ** I didn't try to figure out exactly *why* it didn't
                ** work, because I ran out of time, but I do know there was
                ** no apparent interaction with the is_lethal code (because
                ** I ran with the -nolethal option).
                */
                for (n=i; n <= K->MI; n++) {
                    for (v=0; v <= j; v++) {
                        if (K->is_lethal && v > 0 && K->createlethal == 0)
                            break;  /* don't look at classes that died */
                        if (n+v < j || n+v < i+j)
                            continue;
                        if ((loadclass = from[n][v][xi]) == 0.0) 
                            continue;
                        IF_DEBUG(DEBUG_TRACE3) 
                            fprintf(stderr, "%s: [%d][%d] <- [%d][%d]\n",
                                   thisfunction, i, j, n,v);
                        sum += s_self(i, j, n, v) * 
                               K->S[g] * loadclass;
                    }
                }
                /* final assignment to x"(s)i,j(g) */
                to[i][j][g] = sum;
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        apply_self_progeny_stats    (KConfig K, 
                                         KArray& to, KScalar fromval,
                                         KInt fi, KInt fj, KInt fg)
/*
** Compute proportions of selfed genotypes produced by population.
** This follows Charlesworth et al. 1990.
** K->xpps is correct on exit
*/
{
    const char* thisfunction = "apply_self_progeny_stats";
    KInt i, j, g;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    for (g=0; g < K->genotypes; g++) {
        /* g is the destination genotype */
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                to[i][j][g] = s_self(i, j, fi, fj) * fromval;
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_self_progeny_kondrashov (KConfig K)
/*
** Compute proportions of selfed genotypes produced by population.
** This is related to Kondrashov 1985, p.640, eq. (2), x"i(k) part,
**     term on right side of "+" involving s().
** K->xpps is correct on exit
*/
{
    const char* thisfunction = "compute_self_progeny_kondrashov";
    KScalar sum, innersum, t1;
    KInt i, j, n, v, g, xi;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "i,j =");
    for (g=0; g < K->genotypes; g++) {
        /* g is the destination genotype */
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                /* L(i,j) is the destination load class */
                IF_DEBUG(DEBUG_TRACE2) {
                    fprintf(stderr, " %d,%d ", i, j);
                }
                IF_DEBUG(DEBUG_TRACE3) {
                    if (i%10 == 0) 
                        fprintf(stderr, "%s: [%d][%d] <-\n", 
                               thisfunction, i, j);
                }
                if (K->S[g] == 0.0) {
                    K->xpps[i][j][g] = 0.0;
                    continue;
                }
                sum = 0.0;
                for (n=0; n <= K->MI; n++) {
                    for (v=0; v <= K->MJ; v++) {
                        IF_DEBUG(DEBUG_TRACE3) 
                            fprintf(stderr, "%s: [%d][%d] <- [%d][%d]\n",
                                   thisfunction, i, j, n,v);
                        if (i <= n && j >= v && 
                            j <= n+v && i+j <= n+v) {
                            t1 = K->gamma[n][v] * 
                                 s_self(i, j, n, v);
                            innersum = 0.0;
                            for (xi=0; xi < K->genotypes; xi++) {
                                innersum += K->SO[n][v][xi] * 
                                            K->trfm_S[xi][g];
                            }
                            sum += t1 * innersum;
                        }
                    }
                }
                /* final assignment to x"(s)i,j(g) */
                K->xpps[i][j][g] = sum;
            }
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_apomixis_progeny    (KConfig K)
/*
** Compute proportions of apomixis genotypes produced by population.
** This follows Charlesworth et al. 1990.
** K->xpps is correct on exit
*/
{
    const char* thisfunction = "compute_apomixis_progeny";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_apomixis_progeny(K, K->xppa, K->xp);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_apomixis_progeny  (KConfig K, KArray& to,
                                                KArray& from)
/*
** Compute proportions of apomixis genotypes produced by population.
** This follows Charlesworth et al. 1990.
** K->xpps is correct on exit
*/
{
    const char* thisfunction = "apply_apomixis_progeny";
    KInt i, j, n, v, g, xi;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    for (g=0; g < K->genotypes; g++) {
        /* g is the destination genotype */
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                if (K->is_lethal && j > 0 && K->createlethal == 0) {
                    /* don't look at/create classes that died/will die */
                    to[i][j][g] = 0.0;
                    continue;
                }
                /* L(i,j) is the destination load class */
                if (K->A[g] == 0.0) {
                    to[i][j][g] = 0.0;
                    continue;
                }
                /* note that for this definition of apomixis, 
                ** all we care about is i==n, j==v, g==xi
                */
                n = i;
                v = j;
                xi = g;
                /* final assignment to x"(a)i,j(g) */
                to[i][j][g] = a_apomixis(i, j, n, v) *
                              K->A[xi] *
                              from[n][v][xi];
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_apomixis_progeny_kondrashov (KConfig K)
/*
** Compute proportions of selfed genotypes produced by population.
** This is related to Kondrashov 1985, p.640, eq. (2), x"i(k) part,
**     term on right side of "+" involving s().
** K->xppo is correct on exit
*/
{
    const char* thisfunction = "compute_apomixis_progeny_kondrashov";
    KInt i, j, n, v, g, xi;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "i,j =");
    for (g=0; g < K->genotypes; g++) {
        /* g is the destination genotype */
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                /* L(i,j) is the destination load class */
                IF_DEBUG(DEBUG_TRACE2) {
                    fprintf(stderr, " %d,%d ", i, j);
                }
                IF_DEBUG(DEBUG_TRACE3) {
                    if (i%10 == 0) 
                        fprintf(stderr, "%s: [%d][%d] <-\n", 
                               thisfunction, i, j);
                }
                if (K->A[g] == 0.0) {
                    K->xppa[i][j][g] = 0.0;
                    continue;
                }
                /* note that for this definition of apomixis, 
                ** all we care about is i==n, j==v, g==xi
                */
                n = i;
                v = j;
                xi = g;
                /* final assignment to x"(a)i,j(g) */
                K->xppa[i][j][g] = K->alpha[n][v] * 
                                   a_apomixis(i, j, n, v) *
                                   K->AO[n][v][xi] *
                                   K->trfm_A[xi][g];
            }
        }
    }
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_outcross_progeny_kondrashov (KConfig K)
/*
** Compute proportions of outcrossed genotypes
** This is related to Kondrashov 1985, p.640, eq. (2), 
**     x"i(k) part, term involving b()
** K->xppo is correct on exit
*/
{
    const char* thisfunction = "compute_outcross_progeny_kondrashov";
    KInt g, i, j, n, v, l, lam, xi, zeta;
    KScalar sum, innersum, t1, t2;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "i =");
    for (g=0; g < K->genotypes; g++) {
        for (i=0; i <= K->MI; i++) {
            /* note for outcrossing we only consider j==0 */
            for (j=0; j <= 0; j++) {
                IF_DEBUG(DEBUG_TRACE2) {
                    fprintf(stderr, " %d", i);
                }
                IF_DEBUG(DEBUG_TRACE3) {
                    fprintf(stderr, "%s: [%d][%d] <-\n",
                           thisfunction, i, j);
                }
                if (K->O[g] == 0.0) {
                    /* skip if this genotype does not outcross */
                    K->xppo[i][j][g] = 0.0;
                    continue;
                }
                sum = 0.0;
                for (n=0; n <= K->MI; n++) {
                    for (v=0; v <= K->MJ; v++) {
                        /* L(n,v) represents 'ovule' load class */
                        IF_DEBUG(DEBUG_TRACE3) {
                            fprintf(stderr, "%s: [%d][%d] <- [%d][%d] x ...\n",
                                   thisfunction, i, j, n, v);
                        }
                        if (K->beta[n][v] == 0.0) {
                            /* skip if none of the population is
                            ** allocating outcross ovules to this 
                            ** load class
                            */
                            continue;
                        }
                        for (l=0; l <= K->MI; l++) {
                            for (lam=0; lam <= K->MJ; lam++) {
                                /* L(l,lam) represents 'pollen' 
                                ** load class 
                                */
                                IF_DEBUG(DEBUG_TRACE3) {
                                    fprintf(stderr, "%s: [%d][%d] <- [%d][%d] x [%d][%d]\n",
                                           thisfunction, i, j, 
                                           n, v, l, lam);
                                }
                                if (K->rho[l][lam] == 0.0) {
                                    /* skip if none of the population
                                    ** is allocating outcross pollen to
                                    ** this load class
                                    */
                                    continue;
                                }
                                if ((t1 = o_outcross(i,j,n,v,l,lam)) == 0.0) {
                                    /* skip if no load class will 
                                    ** allocate here or produce progeny 
                                    ** here; what happens with each 
                                    ** genotype doesn't matter if t1==0.0
                                    */
                                    continue;
                                }
                                t1 *= K->beta[n][v] * K->rho[l][lam];
                                innersum = 0.0;
                                for (xi=0; xi < K->genotypes; xi++) {
                                    for (zeta=0; zeta < K->genotypes; zeta++) {
                                        t2 = K->OO[n][v][xi] * 
                                             K->OP[l][lam][zeta];
                                        innersum += t2 * 
                                                    K->trfm_O[xi][zeta][g];
                                    }
                                }
                                sum += t1 * innersum;
                            }
                        }
                    }
                }
                /* final assignment to x"(o)i,j(g) */
                K->xppo[i][j][g] = sum;
            }
            IF_DEBUG(DEBUG_TRACE3) fprintf(stderr, "%s: to [%d][>0]\n",
                                          thisfunction, i);
            for (j=1; j <= K->MJ; j++) {
                /* not possible to produce homozygotes 
                ** via outcrossing
                */
                K->xppo[i][j][g] = 0.0;
            }
        }
        IF_DEBUG(DEBUG_TRACE2) fprintf(stderr, "\n");
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_gametes     (KConfig K)
{
    const char* thisfunction = "compute_gametes";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_gametes(K, K->mgam, K->fgam, K->xp);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_gametes       (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from)
/*
** This is the function as used during normal model operations.
** The gamete proportions are computed from all source classes,
**     then any adjustments are made in a separate function,
**     which receives all the same information as the original
**     function.
*/
{
    apply_gametes_full(K, mgam, fgam, from);
    adjust_gametes(K, mgam, fgam, from);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_gametes_full  (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from)
/*
** Compute genotype proportions for male and female gametes.
** Compute proportion of gametes carrying i mutations, where 
**     i <= K->MI.
** The algorithm is O(K->MI^3).
** This implements Charlesworth et al. 1990, p.1472, eq. (3). 
** mgam and fgam vectors are correct upon exit.
*/
{
    const char* thisfunction = "apply_gametes_full";
    KInt i, n, v, g, xi;
    KScalar sum;
    KScalar t1, t2;
    /* the process of segregation is the same for both gamete 
    ** pools
    */
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    g = 0;
    xi = 0;
    for (i=0; i <= K->MI; i++) {
        sum = 0.0;
        for (v=0; v <= i; v++) {
            if (K->is_lethal && v > 0 && K->createlethal == 0)
                /* When the mutation class is lethal, it's not
                ** possible to have homozygous mutant adults
                ** that contribute to the gamete pools
                */
                break;  /* this ends the for (v=0;;) loop */
            if (v > K->MJ)
                /* It's not possible to have a contribution from
                ** homozogyous loci greater than the number of 
                ** such loci.
                */
                break;  /* this ends the for (v=0;;) loop */
            for (n=(i - v); n <= K->MI; n++) {
                if (from[n][v][xi] == 0.0)
                    continue;
                /* log space */
                /* t1 = lnpow_half(i - v);  << this is a typo in the Charlesworth et al paper */
                t1 = lnpow_half(n);
                t2 = lnbinomial(n, i - v);
                /*
                t2 = lnfactorial(n) - lnfactorial(i - v) - 
                     lnfactorial(n + v - i);
                */
                sum += exp(t1 + t2) * from[n][v][xi];
                /*
                ** linear space -- overflows easily
                ** t1 = pow_half(i - v);
                ** t2 = factorial(n) /
                **      (factorial(i - v) * 
                **      factorial(n + v - i));
                ** sum += (t1 * t2 * K->qp[n][v]);
                */
            }
        }
        /* NO SPECIAL PROCESSING OF GAMETES SHOULD OCCUR HERE.
        ** THAT MUST ALL BE DONE IN adjust_gametes()
        */
        if ((void*)mgam != NULL)
            mgam[i] = sum;
        if ((void*)fgam != NULL)
            fgam[i] = sum;
    }
}

/*///////////////////////////////////////////////////////////////*/
void        adjust_gametes      (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from)
/*
** Adjusts the gamete ratios according to whatever adjustments
**     are in use.
*/
{
    const char* thisfunction = "adjust_gametes";
    KInt i, g;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    g = 0;
    for (i=0; i <= K->MI; i++) {
        /* any special processing for proportions of 
        ** male & female gametes can occur here.
        */
        if ((void*)mgam != NULL) {
            /* mgam[i] = mgam[i]; */
        }
        if ((void*)fgam != NULL) {
            fgam[i] *= K->O[g];
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        apply_gametes_stats (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KScalar fromval,
                                 KInt fi, KInt fj, KInt fg)
/*
** Compute genotype proportions for male and female gametes.
** Compute proportion of gametes carrying i mutations, where 
**     i <= K->MI.
** The algorithm is O(K->MI^3).
** This implements Charlesworth et al. 1990, p.1472, eq. (3). 
** mgam and fgam vectors are correct upon exit.
** [fi][fj][fg] represents the class from which to generate
**     the gamete distribution.  If any are <0, this is
**     identical to apply_gametes() except fgam[] is not
**     scaled by the outcrossing rate.
*/
{
    const char* thisfunction = "apply_gametes_stats";
    KInt i;
    KScalar t1, ans;
    /* the process of segregation is the same for both gamete 
    ** pools
    */
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    /* we have specified the source class explicitly */
    for (i=0; i <= K->MI; i++) {
        if (fi < i - fj || fj > i || fj > K->MJ) {
            ans = 0.0;
        } else {
            t1 = lnpow_half(fi) + lnbinomial(fi, i - fj);
            ans = exp(t1) * fromval;
        }
        if ((void*)mgam != NULL)
            mgam[i] = ans;
        if ((void*)fgam != NULL)
            fgam[i] = ans;
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_zygotes     (KConfig K)
{
    const char* thisfunction = "compute_zygotes";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_zygotes(K, K->xppo, K->mgam, K->fgam);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_zygotes       (KConfig K, KArray& to,
                                 KVector1& mgam, KVector1& fgam)
/*
** TODO: make into genotypes-based with O_transform()
** Compute proportion of outcrossed zygotes produced by population.
** The algorithm is O(K->MI * K->MJ)
** K->qppx is correct on exit.
*/
{
    const char* thisfunction = "apply_zygotes";
    KInt i, j, k, g;
    KScalar sum;    
    /* Combines proportions of male and female gametes into
    ** expected genotypic classes.  Note that all classes 
    ** L(,j>0) = 0.
    */
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    g = 0;
    for (i=0; i <= K->MI; i++) {
        sum = 0.0;
        for (k=0; k <= i; k++) {
            sum += ((mgam[k] * fgam[i - k]) +
                    (mgam[i - k] * fgam[k]));
        }
        /* Gamete pools can each sum to one (selfing and apomixis
        ** decrease the sum of K->fgam) so when we combine them
        ** into the zygote pool, we have to divide by 2 to get the
        ** proportions of zygotes.
        */
        sum *= 0.5;
        to[i][0][g] = sum;
        for (j=1; j <= K->MJ; j++) {
            to[i][j][g] = 0.0;
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_summed_progeny  (KConfig K)
{
    const char* thisfunction = "compute_summed_progeny";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_summed_progeny(K, K->xpp, K->xpps, K->xppa, K->xppo);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_summed_progeny    (KConfig K, KArray& to,
                                     KArray& froms,
                                     KArray& froma,
                                     KArray& fromo)
{
    const char* thisfunction = "apply_summed_progeny";
    KInt i, j, g;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                to[i][j][g] = 0.0;
                if (froms != NULL)
                    to[i][j][g] += froms[i][j][g];
                if (froma != NULL)
                    to[i][j][g] += froma[i][j][g];
                if (fromo != NULL)
                    to[i][j][g] += fromo[i][j][g];
            }
        }
    }
    IF_DEBUG(DEBUG_NORMALIZATION)
        if (K->is_lethal && K->S[0] > 0.0 && K->createlethal == 0)
            /* so that we can produce homozygous lethal progeny */
            fprintf(stderr, "%s: lethal mut class, normalized array not expected here\n", thisfunction);
    check_normalization(K, to, thisfunction, "to");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_auxiliary_values    (KConfig K)
{
    const char* thisfunction = "compute_auxiliary_values";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    not_implemented(thisfunction, "not running in kondrashov mode");
}

/*///////////////////////////////////////////////////////////////*/
void        compute_auxiliary_values_kondrashov (KConfig K)
{
    const char* thisfunction = "compute_auxiliary_values_kondrashov";
    KInt g, i, j;
    KScalar t1, sum, SOsum, AOsum, OOsum, OPsum;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    sum = 0.0;
    for (g=0; g < K->genotypes; g++) {
        t1 = K->rsrc_SO[g] + K->rsrc_AO[g] + K->rsrc_OO[g];
        for (i=0; i <= K->MI; i++) {
            for (j=0; j <= K->MJ; j++) {
                sum += K->x[i][j][g] * t1;
            }
        }
    }
    K->F_female = sum;
    K->F_male = 1 - K->F_female;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            SOsum = AOsum = OOsum = OPsum = 0.0;
            for (g=0; g < K->genotypes; g++) {
                t1 = K->xp[i][j][g];
                SOsum += t1 * K->rsrc_SO[g];
                AOsum += t1 * K->rsrc_AO[g];
                OOsum += t1 * K->rsrc_OO[g];
                OPsum += t1 * K->rsrc_OP[g];
            }
            K->gamma[i][j] = SOsum / K->F_female;
            K->alpha[i][j] = AOsum / K->F_female;
            K->beta[i][j] = OOsum / K->F_female;
            K->rho[i][j] = OPsum / K->F_male;
            for (g=0; g < K->genotypes; g++) {
                if (SOsum != 0.0)
                    K->SO[i][j][g] = K->xp[i][j][g] * 
                                     K->rsrc_SO[g] / SOsum;
                else
                    K->SO[i][j][g] = 0.0;
                if (AOsum != 0.0)
                    K->AO[i][j][g] = K->xp[i][j][g] * 
                                     K->rsrc_AO[g] / AOsum;
                else
                    K->AO[i][j][g] = 0.0;
                if (OOsum != 0.0)
                    K->OO[i][j][g] = K->xp[i][j][g] * 
                                     K->rsrc_OO[g] / OOsum;
                else
                    K->OO[i][j][g] = 0.0;
                if (OPsum != 0.0)
                    K->OP[i][j][g] = K->xp[i][j][g] * 
                                     K->rsrc_OP[g] / OPsum;
                else
                    K->OP[i][j][g] = 0.0;
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
KScalar     s_self              (KInt i, KInt j, KInt n, KInt v)
/*
** Compute selfing proportion produced by given class.
** Proportion of L(i,j) produced when L(n,v) is selfed.
** The proportion of L(i,j) is returned.
** No need to know population data, works solely from i, j, n, v.
*/
{
    //const char* thisfunction = "s_self";
    KScalar ans;
    KScalar t1, t2, t3;
    /* check for classes that cannot be produced through selfing */
    /* the final impossible condition (i+j > n+v) was not made explicit
    **    by Kondrashov, but it is sensible -- you can't have more mutant
    **    loci (hets plus homs) after selfing than existed before selfing */
    if (i > n || j < v || j > n+v || i+j > n+v)
        return 0.0;
    /* if the class can be produced, then compute its expected proportion */
    /* 
    ** linear scale
    ** t1 = binomial(n, i);
    ** t2 = binomial((n - i), (j - v));
    ** t3 = pow_half(2*n - i);
    ** ans = t1 * t2 * t3;
    */
    /* log scale */
    t1 = lnbinomial(n, i);
    t2 = lnbinomial((n - i), (j - v));
    t3 = lnpow_half(2*n - i);
    ans = exp(t1 + t2 + t3);
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     a_apomixis          (KInt i, KInt j, KInt n, KInt v)
/*
** Compute apomictic proportion produced by given class.
** Proportion of L(i,j) produced when L(n,v) undergoes apomixis.
** The proportion of L(i,j) is returned.
** No need to know population data, works solely from i, j, n, v.
** Apomixis duplicates the parent's genotype exactly.
*/
{
    //const char* thisfunction = "a_apomixis";
    if (i == n && j == v)
        return 1.0;
    else
        return 0.0;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     o_outcross          (KInt i, KInt j, KInt n, KInt v, 
                                 KInt l, KInt lam)
/*
** Compute outcross proportion produced by given mating.
** Proportion of L(i,j) produced when L(n,v) mates with L(l,lam).
** Defined and used by Kondrashov, not Charlesworth et al.
** This algorithm is O(1)
** No need to know population data, works solely from i, j, n, v, l, lam.
** The proportion of L(i,j) is returned.
*/
{
    const char* thisfunction = "o_outcross";
    KScalar ans;
    KScalar t1, t2;
    /* check for classes that cannot be produced through outcrossing */
    if (i < (v + lam) ||
        i > (n + l + v + lam) ||
        j != 0) {
        ans = 0.0;
        /*
        IF_DEBUG(DEBUG_OUTCROSS) {
            fprintf(stderr, "o_outcross: L(%d,%d) <- L(%d,%d) X L(%d,%d) = %lg\n",
                   i, j, n, v, l, lam, ans);
        }
        */
        return ans;
    }
    /* if the class can be produced, then compute its expected 
    ** proportion
    */
    /*
    ** linear scale
    ** t1 = binomial((n + l), (i - v - lam));
    ** t2 = pow_half(n + l);
    ** ans = t1 * t2;
    */
    /* log scale */
    t1 = lnbinomial((n + l), (i - v - lam));
    t2 = lnpow_half(n + l);
    ans = exp(t1 + t2);
    IF_DEBUG(DEBUG_OUTCROSS) {
        fprintf(stderr, "%s: L(%d,%d) <- L(%d,%d) X L(%d,%d) = %lg\n",
               thisfunction, i, j, n, v, l, lam, ans);
    }
    return ans;    
}

