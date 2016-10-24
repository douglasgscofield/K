#include "K.h"

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// Routines for computing/displaying statistics of model results
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

static      struct struct_KStats    KStats;

//////////////////////////////////////////////////////////////////
void        stats_print             (KConfig K)

// Print all fields in the KStats structure

{
    if (K->option_table) {
        stats_print_table(K);
    } else {
        stats_print_verbose(K);
    }
}

//////////////////////////////////////////////////////////////////
void        stats_print_table       (KConfig K)

// Print all fields in the KStats structure, as a table

{
    /* stats_print_table_heading(K); */
    cout << (K->generation > GENERATION_CUTOFF) ? "x" : "+";
    cout << sep << K->generation;
    cout << sep << K->U;
    cout << sep << K->S[0];
    cout << sep << K->A[0];
    cout << sep << K->O[0];
    cout << sep << get_fitness_function_name(K->fitness_function);
    cout << sep << K->fit_s;
    cout << sep << K->fit_h;
    cout << sep << KStats.mean_hetloci;
    cout << sep << KStats.var_hetloci;
    cout << sep << KStats.mean_homloci;
    cout << sep << KStats.var_homloci;
    cout << sep << KStats.mean_totmuts;
    cout << sep << KStats.var_totmuts;
    cout << sep << KStats.var_to_mean_totmuts_ratio;
    cout << sep << KStats.mean_fitness_self_progeny;
    cout << sep << KStats.mean_fitness_apomixis_progeny;
    cout << sep << KStats.mean_fitness_outcross_progeny;
    cout << sep << KStats.population_mean_fitness;
    cout << sep << KStats.inbreeding_depression;
    cout << sep << KStats.secondary_selfing_rate;
    cout << endl;
}

//////////////////////////////////////////////////////////////////
void        stats_print_table_heading   (KConfig K)

// Print heading for stats produced

{
    cout << "converge";
    cout << sep << "gen";
    cout << sep << "U";
    cout << sep << "S_0";
    cout << sep << "A_0";
    cout << sep << "O_0";
    cout << sep << "fitfunc";
    cout << sep << "fit_s";
    cout << sep << "fit_h";
    cout << sep << "m_het";
    cout << sep << "v_het";
    cout << sep << "m_hom";
    cout << sep << "v_hom";
    cout << sep << "m_totmuts";
    cout << sep << "v_totmuts";
    cout << sep << "v_m_totmuts";
    cout << sep << "w_self";
    cout << sep << "w_apom";
    cout << sep << "w_outc";
    cout << sep << "w_m_pop";
    cout << sep << "IBD";
    cout << sep << "S_sec;";
    cout << endl;
}

//////////////////////////////////////////////////////////////////
void        stats_print_verbose     (KConfig K)

// Print all fields in the KStats structure, verbosely

{
    cout << "stats begin ========================================" << endl;
    cout << "generations = " << K->generation << " (out of " << GENERATION_CUTOFF << " max)" << endl;
    cout << "fitness function = " << get_fitness_function_name(K->fitness_function) 
        << "  option_truncate = " << K->option_truncate << "  option_nolethal = " << K->option_nolethal
        << endl;
    cout << "U = " << K->U << "  S[0] = " << K->S[0] << "  A[0] = " << K->A[0] << "  O[0] = " << K->O[0] << endl;
    cout << "fit_s = " << K->fit_s << "  fit_h = " << K->fit_h << endl;
    cout << "hetloci mean = " << KStats.mean_hetloci << "  var = " << KStats.var_hetloci << endl;
    cout << "homloci mean = " << KStats.mean_homloci << "  var = " << KStats.var_homloci << endl;
    cout << "totmuts mean = " << KStats.mean_totmuts << "  var = " << KStats.var_totmuts << "  var / mean totmuts = " << KStats.var_to_mean_totmuts_ratio << endl;
    cout << "mean_fitness self_progeny = " << KStats.mean_fitness_self_progeny << "  apomixis_progeny = " << KStats.mean_fitness_apomixis_progeny << "  outcross_progeny = " << KStats.mean_fitness_outcross_progeny << endl;
    cout << "population_mean_fitness = " << KStats.population_mean_fitness << "  inbreeding_depression = " << KStats.inbreeding_depression << "  secondary_selfing_rate = " << KStats.secondary_selfing_rate << endl;
    cout << "stats end ========================================" << endl;
}

//////////////////////////////////////////////////////////////////
void        stats_all               (KConfig K)

// Compute inbreeding depression (1 - self/outcross) from
// the load class date for selfed and outcrossed progeny.

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

//////////////////////////////////////////////////////////////////
KScalar     stats_mean_fitness_self_progeny (KConfig K)

// Compute mean fitness of self progeny

{
    KScalar wmean, sum;
    if ((sum = sum_KArray(K, K->xpps)) == 0.0) {
        /* no self progeny were produced, so we have to
        ** create some and examine their fitness.
        */
    }
    wmean = mean_fitness(K, K->xpps);

    return mean_fitness(K, K->xpps);
}

//////////////////////////////////////////////////////////////////
KScalar     stats_mean_fitness_apomixis_progeny (KConfig K)

// Compute mean fitness of apomixis progeny

{
    return mean_fitness(K, K->xppa);
}

//////////////////////////////////////////////////////////////////
KScalar     stats_mean_fitness_outcross_progeny (KConfig K)

// Compute mean fitness of outcross progeny

{
    return mean_fitness(K, K->xppo);
}

//////////////////////////////////////////////////////////////////
KScalar     stats_inbreeding_depression (KConfig K)

// Compute inbreeding depression (1 - self/outcross) from
// the load class data for selfed and outcrossed progeny.

{
    const char * thisfunction = "stats_inbreeding_depression";
    KInt i, j, g;
    KArray a, via_fgam, via_mgam;
    KVector1 thismgam, thisfgam, allmgam, allfgam;
    KScalar w_self, w_out, ibd, thisibd;
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

//////////////////////////////////////////////////////////////////
KScalar     stats_inbreeding_depression_old (KConfig K)

// Compute inbreeding depression (1 - self/outcross) from
// the load class date for selfed and outcrossed progeny.

{
    const char * thisfunction = "stats_inbreeding_depression_old";
    KScalar ibd;
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

//////////////////////////////////////////////////////////////////
KScalar     stats_secondary_selfing_rate    (KConfig K)

// Compute secondary selfing rate
//     primary selfing rate * (self_w_mean / outcross_w_mean)

{
    const char * thisfunction = "stats_secondary_selfing_rate";
    KInt g;
    KScalar ans;
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

///////////////////////////////////////////////////////////////////
KScalar     stats_population_mean_fitness   (KConfig K)
/*
** Compute mean fitness of population
*/
{
    const char * thisfunction = "stats_population_mean_fitness";
    KScalar ans;
    KInt g;
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

///////////////////////////////////////////////////////////////////
void        stats_muts      (KConfig K)
/*
** Compute statistics for the number of mutant
** alleles carried by each adult plant.
*/
{
    //const char * thisfunction = "stats_muts";
    KScalar mean_hetloci, mean_homloci, mean_totmuts;
    KScalar var_hetloci, var_homloci, var_totmuts, t1;
    KInt hetloci, homloci, totmuts;
    KInt i, j, g;
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


/*
///////////////////////////////////////////////////////////////
KScalar     stats_variance_lethals  (KConfig K)

// Compute the variance of the number of lethal
// alleles carried by each adult plant.

{
    const char * thisfunction = "stats_variance_lethals";
    KScalar variance_lethals, t1, t2, t3;
    KInt i, j, g;
    if (KStats.mean_totmuts < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    //fprintf(stderr, "i\tj\tg\tx(ijg)\ti+j*2\tx\t(x-meanx)\tcumsum\t()^2\tcumsum^2\n");
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
                //fprintf(stderr, "%d\t%d\t%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\n",
                //    i, j, g, K->x[i][j][g], i+j*2, t1, t2, t3, t4, variance_lethals);
            }
        }
    }
    return variance_lethals;
}
*/

