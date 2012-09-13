#include "K.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for computing/displaying statistics of model results */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

static      struct struct_KStats    KStats;

/*///////////////////////////////////////////////////////////////*/
void        stats_print             (KConfig K)
/*
** Print all fields in the KStats structure
*/
{
    if (K->option_table) {
        stats_print_table(K);
    } else {
        stats_print_verbose(K);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_table       (KConfig K)
/*
** Print all fields in the KStats structure, as a table
*/
{
    /* stats_print_table_heading(K); */
    printf("%s %d\t%lg\t%lg\t%lg\t%lg\t%s\t\
%lg\t%lg\t\
%lg\t%lg\t\
%lg\t%lg\t\
%lg\t%lg\t%lg\t\
%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
           (K->generation > GENERATION_CUTOFF) ? "*" : " ",
           K->generation,
           K->U,
           K->S[0],
           K->A[0],
           K->O[0],
           get_fitness_function_name(K->fitness_function),
           K->fit_s,
           K->fit_h,
           KStats.mean_hetloci,
           KStats.var_hetloci,
           KStats.mean_homloci,
           KStats.var_homloci,
           KStats.mean_totmuts,
           KStats.var_totmuts,
           KStats.var_to_mean_totmuts_ratio,
           KStats.mean_fitness_self_progeny,
           KStats.mean_fitness_apomixis_progeny,
           KStats.mean_fitness_outcross_progeny,
           KStats.population_mean_fitness,
           KStats.inbreeding_depression,
           KStats.secondary_selfing_rate
           );
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_table_heading   (KConfig K)
/*
** Print heading for stats produced
*/
{
    printf("generations\tU\tS_0\tA_0\tO_0\tfitness_function\t\
fit_s\tfit_h\tmean_hetloci\tvar_hetloci\tmean_homloci\tvar_homloci\tmean_totmuts\tvar_totmuts\tvar_mean_totmuts\t\
w_self\tw_apomixis\tw_outcross\tw_popmean\tIBD\tS_secondary\n");
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_verbose     (KConfig K)
/*
** Print all fields in the KStats structure, verbosely
*/
{
    printf("stats begin ========================================\n");
    printf("generations = %d (out of %d max)\n", K->generation, 
           GENERATION_CUTOFF);
    printf("U = %lg\n", K->U);
    printf("S[0] = %lg\n", K->S[0]);
    printf("A[0] = %lg\n", K->A[0]);
    printf("O[0] = %lg\n", K->O[0]);
    printf("fitness function = %s\n",
           get_fitness_function_name(K->fitness_function));
    printf("fit_s = %lg\n", K->fit_s);
    printf("fit_h = %lg\n", K->fit_h);
    printf("option_truncate = %d\n", K->option_truncate);
    printf("option_nolethal = %d\n", K->option_nolethal);

    printf("mean_hetloci = %lg\n", KStats.mean_hetloci);
    printf("var_hetloci = %lg\n", KStats.var_hetloci);
    printf("mean_homloci = %lg\n", KStats.mean_homloci);
    printf("var_homloci = %lg\n", KStats.var_homloci);
    printf("mean_totmuts = %lg\n", KStats.mean_totmuts);
    printf("var_totmuts = %lg\n", KStats.var_totmuts);
    printf("var / mean totmuts = %lg\n",
           KStats.var_to_mean_totmuts_ratio);
    printf("mean_fitness_self_progeny = %lg\n", 
           KStats.mean_fitness_self_progeny);
    printf("mean_fitness_apomixis_progeny = %lg\n", 
           KStats.mean_fitness_apomixis_progeny);
    printf("mean_fitness_outcross_progeny = %lg\n", 
           KStats.mean_fitness_outcross_progeny);
    printf("population_mean_fitness = %lg\n", 
           KStats.population_mean_fitness);
    printf("inbreeding_depression = %lg\n", 
           KStats.inbreeding_depression);
    printf("secondary_selfing_rate = %lg\n", 
           KStats.secondary_selfing_rate);
    printf("stats end ========================================\n");
}

/*///////////////////////////////////////////////////////////////*/
void        stats_all               (KConfig K)
/*
** Compute inbreeding depression (1 - self/outcross) from
** the load class date for selfed and outcrossed progeny.
*/
{
    KStats.mean_hetloci = -999.9;
    KStats.var_hetloci = -999.9;
    KStats.mean_homloci = -999.9;
    KStats.var_homloci = -999.9;
    KStats.mean_totmuts = -999.9;
    KStats.var_totmuts = -999.9;
    KStats.var_to_mean_totmuts_ratio = -999.9;
    KStats.mean_fitness_self_progeny = -999.9;
    KStats.mean_fitness_apomixis_progeny = -999.9;
    KStats.mean_fitness_outcross_progeny = -999.9;
    KStats.population_mean_fitness = -999.9;
    KStats.inbreeding_depression = -999.9;
    KStats.secondary_selfing_rate = -999.9;

    stats_muts(K);
    if (KStats.mean_totmuts == 0.0) {
        KStats.var_to_mean_totmuts_ratio = 0.0;
    } else {
        KStats.var_to_mean_totmuts_ratio = 
            KStats.var_totmuts / KStats.mean_totmuts;
    }
    KStats.mean_fitness_self_progeny = 
        stats_mean_fitness_self_progeny(K);
    KStats.mean_fitness_apomixis_progeny = 
        stats_mean_fitness_apomixis_progeny(K);
    KStats.mean_fitness_outcross_progeny = 
        stats_mean_fitness_outcross_progeny(K);
    KStats.population_mean_fitness = 
        stats_population_mean_fitness(K);
    KStats.inbreeding_depression = stats_inbreeding_depression_old(K);
    // KStats.inbreeding_depression = stats_inbreeding_depression(K);
    KStats.secondary_selfing_rate = stats_secondary_selfing_rate(K);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_self_progeny (KConfig K)
/*
** Compute mean fitness of self progeny
*/
{
    const char* thisfunction = "stats_mean_fitness_self_progeny";
    KScalar wmean, sum;
    if ((sum = sum_KArray(K, K->xpps)) == 0.0) {
        /* no self progeny were produced, so we have to
        ** create some and examine their fitness.
        */
    }
    wmean = mean_fitness(K, K->xpps);

    return mean_fitness(K, K->xpps);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_apomixis_progeny (KConfig K)
/*
** Compute mean fitness of apomixis progeny
*/
{
    const char* thisfunction = "stats_mean_fitness_apomixis_progeny";
    return mean_fitness(K, K->xppa);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_outcross_progeny (KConfig K)
/*
** Compute mean fitness of outcross progeny
*/
{
    const char* thisfunction = "stats_mean_fitness_outcross_progeny";
    return mean_fitness(K, K->xppo);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_inbreeding_depression (KConfig K)
/*
** Compute inbreeding depression (1 - self/outcross) from
** the load class data for selfed and outcrossed progeny.
*/
{
    const char* thisfunction = "stats_inbreeding_depression";
    KInt i, j, g;
    KArray a, via_fgam, via_mgam;
    KVector1 thismgam, thisfgam, allmgam, allfgam;
    KScalar w_self, w_out, ibd, thisibd;
    IF_DEBUG(DEBUG_TRACE1) {
        printf("%s\n", thisfunction);
    }
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    ibd = 0.0;
    apply_gametes_full(K, allmgam, allfgam, K->x);
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                /* compute fitness of selfed and outcrossed
                ** progeny of this load class by creating dummy
                ** arrays representing the progeny resulting
                ** from matings involving just these load classes.
                **
                ** Selfing is straightforward.
                **
                ** For outcrossing, matings involving both female
                ** and male gametes must be considered.  This is
                ** done by mating female gametes produced by this
                ** load class with the population male gametes,
                ** and the male gametes produced by this load
                ** class with  the population female gametes.
                */
                if (K->x[i][j][g] == 0.0)
                    continue;
                /* fill_KArray(K, a1, 0.0);
                ** a1[i][j][g] = K->x[i][j][g];
                ** apply_self_progeny(K, a2, a1);
                */
                apply_self_progeny_stats(K, a, K->x[i][j][g], i, j, g);
                w_self = cumulative_fitness(K, a);
                apply_gametes_stats(K, thismgam, thisfgam, K->x[i][j][g], i, j, g);
                /* apply_gametes(K, allmgam, allfgam, K->x); */
                apply_zygotes(K, via_fgam, allmgam, thisfgam);
                apply_zygotes(K, via_mgam, thismgam, allfgam);
                w_out = cumulative_fitness(K, via_fgam) +
                        cumulative_fitness(K, via_mgam);
                w_out *= 0.5;
                thisibd = 1.0 - (w_self / w_out);
                ibd += thisibd * K->x[i][j][g];
            }
        }
    }
    return ibd;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_inbreeding_depression_old (KConfig K)
/*
** Compute inbreeding depression (1 - self/outcross) from
** the load class date for selfed and outcrossed progeny.
*/
{
    const char* thisfunction = "stats_inbreeding_depression_old";
    KScalar ibd;
    IF_DEBUG(DEBUG_TRACE1) {
        printf("%s\n", thisfunction);
    }
    if (KStats.mean_fitness_self_progeny < 0.0 ||
        KStats.mean_fitness_outcross_progeny < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    if (KStats.mean_fitness_self_progeny == 0.0 ||
        KStats.mean_fitness_outcross_progeny == 0.0) {
        /*
        ** char buf[200];
        ** sprintf(buf, "%s: no self or outcross progeny", thisfunction);
        ** warning(buf);
        */
        ibd = 0.0;
    } else {
        ibd = 1.0 - (KStats.mean_fitness_self_progeny / 
                     KStats.mean_fitness_outcross_progeny);
    }
    return ibd;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_secondary_selfing_rate    (KConfig K)
/*
** Compute secondary selfing rate
**     primary selfing rate * (self_w_mean / outcross_w_mean)
*/
{
    const char* thisfunction = "stats_secondary_selfing_rate";
    KInt g;
    KScalar ans;
    IF_DEBUG(DEBUG_TRACE1) {
        printf("%s\n", thisfunction);
    }
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    if (KStats.mean_fitness_self_progeny < 0.0 ||
        KStats.population_mean_fitness < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    g = 0;
    if (KStats.population_mean_fitness == 0.0) {
        ans = 0.0;
    } else {
        ans = K->S[g] * (KStats.mean_fitness_self_progeny /
                         KStats.population_mean_fitness);
    }
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_population_mean_fitness   (KConfig K)
/*
** Compute mean fitness of population
*/
{
    const char* thisfunction = "stats_mean_fitness";
    KScalar ans;
    KInt g;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
    }
    if (KStats.mean_fitness_self_progeny < 0.0 ||
        KStats.mean_fitness_apomixis_progeny < 0.0 ||
        KStats.mean_fitness_outcross_progeny < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    g = 0;
    ans = (K->S[g] * KStats.mean_fitness_self_progeny) +
          (K->A[g] * KStats.mean_fitness_apomixis_progeny) +
          (K->O[g] * KStats.mean_fitness_outcross_progeny);
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
void        stats_muts      (KConfig K)
/*
** Compute statistics for the number of mutant
** alleles carried by each adult plant.
*/
{
    const char* thisfunction = "stats_muts";
    KScalar mean_hetloci, mean_homloci, mean_totmuts;
    KScalar var_hetloci, var_homloci, var_totmuts, t1;
    KInt hetloci, homloci, totmuts;
    KInt i, j, g;
    IF_DEBUG(DEBUG_TRACE1) {
        printf("%s\n", thisfunction);
    }
    mean_hetloci = 0.0;
    mean_homloci = 0.0;
    mean_totmuts = 0.0;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                hetloci = i;
                homloci = j;
                totmuts = hetloci + 2*homloci;
                mean_hetloci += hetloci * K->x[i][j][g];
                mean_homloci += homloci * K->x[i][j][g];
                mean_totmuts += totmuts * K->x[i][j][g];
            }
        }
    }
    var_hetloci = 0.0;
    var_homloci = 0.0;
    var_totmuts = 0.0;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                hetloci = i;
                homloci = j;
                totmuts = hetloci + 2*homloci;
                t1 = hetloci - mean_hetloci;
                var_hetloci += (t1 * t1) * K->x[i][j][g];
                t1 = homloci - mean_homloci;
                var_homloci += (t1 * t1) * K->x[i][j][g];
                t1 = totmuts - mean_totmuts;
                var_totmuts += (t1 * t1) * K->x[i][j][g];
            }
        }
    }
    KStats.mean_hetloci = mean_hetloci;
    KStats.mean_homloci = mean_homloci;
    KStats.mean_totmuts = mean_totmuts;
    KStats.var_hetloci = var_hetloci;
    KStats.var_homloci = var_homloci;
    KStats.var_totmuts = var_totmuts;
}


/* //////////////////////////////////////////////////////////////
KScalar     stats_variance_lethals  (KConfig K)

// Compute the variance of the number of lethal
// alleles carried by each adult plant.

{
    const char* thisfunction = "stats_variance_lethals";
    KScalar variance_lethals, t1, t2, t3;
    KInt i, j, g;
    IF_DEBUG(DEBUG_TRACE1) {
        printf("%s\n", thisfunction);
    }
    if (KStats.mean_totmuts < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    //printf("i\tj\tg\tx(ijg)\ti+j*2\tx\t(x-meanx)\tcumsum\t()^2\tcumsum^2\n");
    t3 = 0.0;
    variance_lethals = 0.0;
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                // keep a running sum-of-squares
                t1 = (i + j*2) - KStats.mean_lethals;
                t2 = t1 * t1;
                t3 = K->x[i][j][g] * t2;
                variance_lethals += t3;
                //printf("%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\n",
                //    i, j, g, K->x[i][j][g], i+j*2, t1, t2, t3, t4, variance_lethals);
            }
        }
    }
    return variance_lethals;
}
*/

