#include "K.h"


/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for computing and applying mutation                  */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/* NOTE: The routines for Poisson terms are in K_math.c */

/*///////////////////////////////////////////////////////////////*/
void        compute_mutation    (KConfig K)
/* 
** Add mutation to adults for each genotype class held in K->x
**     and place new frequencies in K->xp.
** If you want to add mutations to an arbitrary KArray, then use the
**     subfunction apply_mutation(KConfig K, KArray& to, KArray& from)
** K->xp is correct upon exit.
*/
{
    const char* thisfunction = "compute_mutation";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->U != 0.0) {
        apply_mutation(K, K->xp, K->x);
    } else {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "no mutation\n");
        copy_KArray(K, K->xp, K->x);
    }
    check_normalization(K, K->xp, thisfunction, "K->xp");
}

/*///////////////////////////////////////////////////////////////*/
void        apply_mutation      (KConfig K, KArray& to, KArray& from)
{
    const char* thisfunction = "apply_mutation";
    KInt g, i, j, n1;
    KScalar sum;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    for (g=0; g < K->genotypes; g++) {
        for (j=0; j <= K->MJ; j++) {
            for (i=0; i <= K->MI; i++) {
                sum = 0.0;
                for (n1=0; n1 <= i; n1++) {
                    sum += from[n1][j][g] * 
#if defined(MUTTERM_FUNCTION)
                           mut_term(K, i-n1);
#else
                           K->mut_term[i-n1];
#endif
                }
                to[i][j][g] = sum;
            }
        }
    }
    check_normalization(K, to, thisfunction, "to");
}

/*///////////////////////////////////////////////////////////////*/
void        apply_mutation_general  (KConfig K, KScalar U, 
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
void        initiate_mut_term       (KConfig K)
{
    const char* thisfunction = "initiate_mut_term";
    KInt i;
    KScalar t1;
    KScalar checksum = 0.0;
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    if (K->U == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "no mutation\n");
        for (i=0; i <= K->MI; i++) {
            K->mut_term[i] = pois_term(0.0, i);
        }
        return;
    }
    /* with K->U = 1, K->mut_term[150] = 0.1750276207e-262 */
    /* with K->U = 30, K->mut_term[150] = 0.6475820422e-41 */
    /* with K->U = 50, K->mut_term[150] = 0.1226329680e-7 */
    /* with K->U = 70, K->mut_term[150] = 0.1018151936e15 */
    for (i=0; i <= K->MI; i++) {
        t1 = pois_term(K->U, i);
        K->mut_term[i] = t1;
        checksum += t1;
    }
    /* verify that the contents of the mut_term[] array are OK */
    if (fabs(checksum - 1.0) > 0.000000001) {
        char buf[200];
        sprintf(buf, "%s: not valid pdf, U=%lg: 1 != %lf\n", 
                thisfunction, K->U, checksum);
        warning(buf);
    }
    IF_DEBUG(DEBUG_NORMALIZATION) {
        KInt i;
        KScalar termsum = 0.0;
        for (i=0; i <= K->MI; i++) {
            termsum += K->mut_term[i];
        }
        fprintf(stderr, "%s: U=%lg, sum mut_term=%lf\n", 
               thisfunction, K->U, termsum);
    }
}

/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term                (KConfig K, KInt t)
{
    const char* thisfunction = "mut_term";
    if (t > K->MI) {
        char buf[200];
        sprintf(buf, "%s: term t=%d too large", thisfunction, t);
        fatal(buf);
    }
    return K->mut_term[t];
}

/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term_computed       (KConfig K, KInt t)
{
    const char* thisfunction = "mut_term_computed";
    if (t > K->MI) {
        char buf[200];
        sprintf(buf, "%s: term t=%d too large", thisfunction, t);
        fatal(buf);
    }
    return pois_term(K->U, t);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     mut_term_general        (KScalar U, KInt t)
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
    if (t > MAX_MI) {
        char buf[200];
        sprintf(buf, "%s: term too large", thisfunction);
        fatal(buf);
    }
    if (U == prev_U) {
        return pois_general[prev_index][t];
    } else {
        KInt i;
        KScalar t1;
        KScalar checksum = 0.0;
        /* find an existing entry in the pois_general array */
        for (i=0; i <= last_index; i++) {
            if (pois_general[i][MAX_MI+1] == U) {
                prev_U = U;
                prev_index = i;
                return pois_general[prev_index][t];
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
        return pois_general[prev_index][t];
    }
}

/*///////////////////////////////////////////////////////////////*/
KScalar     pois_term           (KScalar U, KInt term)
/*
** Computes and returns the term-th term of the Poisson 
** distribution with mean U.  Sets term to 0.0 if its
** value is below POISSON_TRUNCATE.  Also modified so
** that if U == 0.0, then the zeroth term is 1.0, and all
** the rest of the terms are 0.0.
*/
{
    const char* thisfunction = "pois_term";
    KScalar ans;
    /* linear space -- overflows ca. i==150
    ** t1 = exp_neg_U * pow(K->U, i) / factorial(i);
    ** alternate log space
    ** t1 = exp(-K->U) * exp(i*log(K->U) - lnfactorial(i));
	*/
    /* log space */
    if (U < 0.0) {
        char buf[200];
        sprintf(buf, "%s: U < 0", thisfunction);
        fatal(buf);
    } else if (U == 0.0) {
        return (term == 0) ? 1.0 : 0.0;
    }
    ans = exp(-U + term*log(U) - lnfactorial(term));
    if (ans < POISSON_TRUNCATE)
        ans = 0.0;
    return ans;
}

