#include "K_n.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Routines for computing/displaying statistics of model results */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

static      struct struct_KStats_n  KStats_n;

/*///////////////////////////////////////////////////////////////*/
void        stats_print_n           (KConfig_n KN)
/*
** Print all fields in the KStats_n structure
*/
{
    const char* thisfunction = "stats_print_n";
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->option_table) {
        stats_print_table_n(KN);
    } else {
        stats_print_verbose_n(KN);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_table_n     (KConfig_n KN)
/*
** Print all fields in the KStats_n structure, as a table
*/
{
    const char* thisfunction = "stats_print_table_n";
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    /* stats_print_table_heading_n(KN); */
    printf("%s %d\t%lg\t%lg\t%lg\t\
%lg\t%s\t%lg\t%lg\t\
%lg\t%s\t%lg\t%lg\t\
%lg\t%lg\t%lg\t\
%lg\t%lg\t%lg\t\
%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
           (KN->generation > GENERATION_CUTOFF) ? "*" : " ",
           KN->generation,
           KN->S,
           KN->A,
           KN->O,  
           /* */
           KN->U[0],
           get_fitness_function_name(KN->fitness_function[0]),
           KN->fit_s[0],
           KN->fit_h[0],
           /* */
           KN->U[1],
           get_fitness_function_name(KN->fitness_function[1]),
           KN->fit_s[1],
           KN->fit_h[1],
           /* */
           KStats_n.mean_hetloci[0],
           KStats_n.var_hetloci[0],
           KStats_n.mean_homloci[0],
           KStats_n.var_homloci[0],
           KStats_n.mean_totmuts[0],
           KStats_n.var_totmuts[0],
           KStats_n.var_to_mean_totmuts_ratio[0],
           /* */
           KStats_n.mean_hetloci[1],
           KStats_n.var_hetloci[1],
           KStats_n.mean_homloci[1],
           KStats_n.var_homloci[1],
           KStats_n.mean_totmuts[1],
           KStats_n.var_totmuts[1],
           KStats_n.var_to_mean_totmuts_ratio[1],
           /* */
           KStats_n.mean_fitness_self_progeny,
           KStats_n.mean_fitness_apomixis_progeny,
           KStats_n.mean_fitness_outcross_progeny,
           KStats_n.population_mean_fitness,
           KStats_n.inbreeding_depression,
           KStats_n.secondary_selfing_rate
           );
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_table_heading_n (KConfig_n KN)
/*
** Print heading for stats produced
*/
{
    const char* thisfunction = "stats_print_table_heading_n";
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    printf("generations\tS\tA\tO\t\
U_0\tfitness_function_0\tfit_s_0\tfit_h_0\t\
U_1\tfitness_function_1\tfit_s_1\tfit_h_1\t\
mean_lethals_0\tvar_lethals_0\tvar_mean_lethals_0\t\
mean_lethals_1\tvar_lethals_1\tvar_mean_lethals_1\t\
w_self\tw_apomixis\tw_outcross\tw_popmean\tIBD\tS_secondary\n");
}

/*///////////////////////////////////////////////////////////////*/
void        stats_print_verbose_n   (KConfig_n KN)
/*
** Print all fields in the KStats_n structure, verbosely
*/
{
    const char* thisfunction = "stats_print_verbose_n";
    KMutClass m;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    printf("stats begin ========================================\n");
    printf("generations = %d (out of %d max)\n", KN->generation, 
           GENERATION_CUTOFF);
    printf("S = %lg\n", KN->S);
    printf("A = %lg\n", KN->A);
    printf("O = %lg\n", KN->O);
    printf("option_truncate = %d\n", KN->option_truncate);
    printf("option_nolethal = %d\n", KN->option_nolethal);
    printf("mutation classes = %d\n", KN->mutclasses);
    for (m=0; m < KN->mutclasses; m++) {
        printf("U[%d] = %lg\n", m, KN->U[m]);
        printf("fitness function[%d] = %s\n", m,
               get_fitness_function_name(KN->fitness_function[m]));
        printf("fit_s[%d] = %lg\n", m, KN->fit_s[m]);
        printf("fit_h[%d] = %lg\n", m, KN->fit_h[m]);
    }
    for (m=0; m < KN->mutclasses; m++) {
        printf("mean_hetloci[%d] = %lg\n", m,
               KStats_n.mean_hetloci[m]);
        printf("var_hetloci[%d] = %lg\n", m,
               KStats_n.var_hetloci[m]);
        printf("mean_homloci[%d] = %lg\n", m,
               KStats_n.mean_homloci[m]);
        printf("var_homloci[%d] = %lg\n", m,
               KStats_n.var_homloci[m]);
        printf("mean_totmuts[%d] = %lg\n", m,
               KStats_n.mean_totmuts[m]);
        printf("var_totmuts[%d] = %lg\n", m,
               KStats_n.var_totmuts[m]);
        printf("var/mean totmuts[%d] = %lg\n", m,
               KStats_n.var_to_mean_totmuts_ratio[m]);
    }
    printf("mean_fitness_self_progeny = %lg\n", 
           KStats_n.mean_fitness_self_progeny);
    printf("mean_fitness_apomixis_progeny = %lg\n", 
           KStats_n.mean_fitness_apomixis_progeny);
    printf("mean_fitness_outcross_progeny = %lg\n", 
           KStats_n.mean_fitness_outcross_progeny);
    printf("population_mean_fitness = %lg\n", 
           KStats_n.population_mean_fitness);
    printf("inbreeding_depression = %lg\n", 
           KStats_n.inbreeding_depression);
    printf("secondary_selfing_rate = %lg\n", 
           KStats_n.secondary_selfing_rate);
    printf("stats end ========================================\n");
}

/*///////////////////////////////////////////////////////////////*/
void        stats_all_n         (KConfig_n KN)
{
    const char* thisfunction = "stats_all_n";
    KMutClass m;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    KStats_n.mean_fitness_self_progeny = -999.9;
    KStats_n.mean_fitness_apomixis_progeny = -999.9;
    KStats_n.mean_fitness_outcross_progeny = -999.9;
    KStats_n.population_mean_fitness = -999.9;
    KStats_n.inbreeding_depression = -999.9;
    KStats_n.secondary_selfing_rate = -999.9;
    for (m=0; m < KN->mutclasses; m++) {
        stats_muts_n(KN, m);
        if (KStats_n.mean_totmuts[m] == 0.0) {
            KStats_n.var_to_mean_totmuts_ratio[m] = 0.0;
        } else {
            KStats_n.var_to_mean_totmuts_ratio[m] = 
                    KStats_n.var_totmuts[m] / 
                    KStats_n.mean_totmuts[m];
        }
    }
    KStats_n.mean_fitness_self_progeny = 
        stats_mean_fitness_self_progeny_n(KN);
    KStats_n.mean_fitness_apomixis_progeny = 
        stats_mean_fitness_apomixis_progeny_n(KN);
    KStats_n.mean_fitness_outcross_progeny = 
        stats_mean_fitness_outcross_progeny_n(KN);
    KStats_n.population_mean_fitness = 
        stats_population_mean_fitness_n(KN);
    KStats_n.inbreeding_depression = 
        stats_inbreeding_depression_old_n(KN);
    KStats_n.secondary_selfing_rate = 
        stats_secondary_selfing_rate_n(KN);
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_self_progeny_n (KConfig_n KN)
/*
** Compute mean fitness of self progeny.  To do this with the
** KN->x1,x2 type of data structures, we have to generate
** self progeny to a temporary KArray_n, then examine the
** genotypes in that array for fitness.
**
** If there were no selfed progeny produced, then we of course
** have zero mean fitness due to selfed progeny.
*/
{
    const char* thisfunction = "stats_mean_fitness_self_progeny_n";
    KScalar wmean;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->S == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) printf("%s: S==0.0, self wmean=0.0\n", 
                                      thisfunction);
        wmean = 0.0;
    } else {
        //void* a = alloc_KArray_n();
        KArray_n a;
        apply_self_progeny_n(KN, a, KN->x1);
        wmean = mean_fitness_n(KN, a);
        //free_KArray_n(a);
    }
    return wmean;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_apomixis_progeny_n (KConfig_n KN)
/*
** Compute mean fitness of apomixis progeny
*/
{
    const char* thisfunction = "stats_mean_fitness_apomixis_progeny_n";
    KScalar wmean;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->A == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) printf("%s: A==0.0, apomict wmean=0.0\n", 
                                      thisfunction);
        wmean = 0.0;
    } else {
        //void* a = alloc_KArray_n();
        KArray_n a;
        apply_apomixis_progeny_n(KN, a, KN->x1);
        wmean = mean_fitness_n(KN, a);
        //free_KArray_n(a);
    }
    return wmean;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_mean_fitness_outcross_progeny_n (KConfig_n KN)
/*
** Compute mean fitness of outcross progeny
*/
{
    const char* thisfunction = "stats_mean_fitness_outcross_progeny_n";
    KScalar wmean;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->O == 0.0) {
        IF_DEBUG(DEBUG_TRACE1) printf("%s: O==0.0, outcross wmean=0.0\n", 
                                      thisfunction);
        wmean = 0.0;
    } else {
        //void* a = alloc_KArray_n();
        //void* vm = alloc_KVector_n();
        //void* vf = alloc_KVector_n();
        KVector_n vm;
        KVector_n vf;
        apply_gametes_n(KN, vm, vf, KN->x1);
        {
            KArray_n a;
            apply_zygotes_n(KN, a, vm, vf);
            wmean = mean_fitness_n(KN, a);
        }
        //free_KVector_n(vf);
        //free_KVector_n(vm);
        //free_KArray_n(a);
    }
    return wmean;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_inbreeding_depression_n   (KConfig_n KN)
/*
** Compute inbreeding depression (1 - self/outcross) from
** the load class data for selfed and outcrossed progeny.
**
** Because inbreeding depression can have a value even if
** selfed progeny aren't being produced, we have to dummy
** up a value of inbreeding depression if that is the case.
** To do this, we compute inbreeding depression on a class
** by class basis by creating outcrossed and selfed progeny
** for each class and examining their fitnesses.
*/
{
    const char* thisfunction = "stats_inbreeding_depression_n";
    KScalar ibd;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    /* for each load class */
    /* compute fitnesses of selfed progeny */
    /* compute fitnesses of outcrossed progeny */
    /* compute inbreeding depression for each load class,
    ** weighted by the load class proportion.
    */
    if (KStats_n.mean_fitness_self_progeny < 0.0 ||
        KStats_n.mean_fitness_outcross_progeny < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    if (KStats_n.mean_fitness_self_progeny == 0.0 ||
        KStats_n.mean_fitness_outcross_progeny == 0.0) {
        /*
        ** char buf[200];
        ** sprintf(buf, "%s: no self or outcross progeny", thisfunction);
        ** warning(buf);
        */
        ibd = 0.0;
    } else {
        ibd = 1.0 - (KStats_n.mean_fitness_self_progeny / 
                     KStats_n.mean_fitness_outcross_progeny);
    }
    return ibd;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_inbreeding_depression_old_n   (KConfig_n KN)
/*
** Compute inbreeding depression (1 - self/outcross) from
** the load class data for selfed and outcrossed progeny.
**
** Because inbreeding depression can have a value even if
** selfed progeny aren't being produced, we have to dummy
** up a value of inbreeding depression if that is the case.
** To do this, we compute inbreeding depression on a class
** by class basis by creating outcrossed and selfed progeny
** for each class and examining their fitnesses.
*/
{
    const char* thisfunction = "stats_inbreeding_depression_n";
    KScalar ibd;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    /* compute fitnesses of selfed progeny */
    /* compute fitnesses of outcrossed progeny */
    /* compute inbreeding depression for each load class,
    ** weighted by the load class proportion.
    */
    if (KStats_n.mean_fitness_self_progeny < 0.0 ||
        KStats_n.mean_fitness_outcross_progeny < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    if (KStats_n.mean_fitness_self_progeny == 0.0 ||
        KStats_n.mean_fitness_outcross_progeny == 0.0) {
        /*
        ** char buf[200];
        ** sprintf(buf, "%s: no self or outcross progeny", thisfunction);
        ** warning(buf);
        */
        ibd = 0.0;
    } else {
        ibd = 1.0 - (KStats_n.mean_fitness_self_progeny / 
                     KStats_n.mean_fitness_outcross_progeny);
    }
    return ibd;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_secondary_selfing_rate_n    (KConfig_n KN)
/*
** Compute secondary selfing rate
**     primary selfing rate * (self_w_mean / outcross_w_mean)
*/
{
    const char* thisfunction = "stats_secondary_selfing_rate_n";
    KScalar ans;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KStats_n.mean_fitness_self_progeny < 0.0 ||
        KStats_n.population_mean_fitness < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    if (KStats_n.population_mean_fitness == 0.0) {
        ans = 0.0;
    } else {
        ans = KN->S * (KStats_n.mean_fitness_self_progeny /
                       KStats_n.population_mean_fitness);
    }
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
KScalar     stats_population_mean_fitness_n   (KConfig_n KN)
/*
** Compute mean fitness of population
*/
{
    const char* thisfunction = "stats_mean_fitness";
    KScalar ans;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KStats_n.mean_fitness_self_progeny < 0.0 ||
        KStats_n.mean_fitness_apomixis_progeny < 0.0 ||
        KStats_n.mean_fitness_outcross_progeny < 0.0) {
        char buf[200];
        sprintf(buf, "%s: wrong stats order", thisfunction);
        fatal(buf);
    }
    ans = (KN->S * KStats_n.mean_fitness_self_progeny) +
          (KN->A * KStats_n.mean_fitness_apomixis_progeny) +
          (KN->O * KStats_n.mean_fitness_outcross_progeny);
    return ans;
}

/*///////////////////////////////////////////////////////////////*/
void        stats_muts_n    (KConfig_n KN, KMutClass m)
/*
** Compute a variety of statistics for mutant alleles carried by 
** each adult plant in the given mutation class
*/
{
    const char* thisfunction = "stats_muts_n";
    KScalar t1;
    KScalar mean_hetloci, mean_homloci, mean_totmuts;
    KScalar var_hetloci, var_homloci, var_totmuts;
    KInt i0, j0, i1, j1, hetloci, homloci, totmuts;
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    mean_hetloci = 0.0;
    mean_homloci = 0.0;
    mean_totmuts = 0.0;
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (j0=0; j0 <= KN->MJ0; j0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                for (j1=0; j1 <= KN->MJ1; j1++) {
                    switch(m) {
                    case 0:
                        hetloci = i0;
                        homloci = j0;
                        break;
                    case 1:
                        hetloci = i1;
                        homloci = j1;
                        break;
                    default:
                        not_implemented(thisfunction, "m > 1");
                        break;
                    }
                    totmuts = hetloci + 2*homloci;
                    mean_hetloci += hetloci * KN->x1[i0][j0][i1][j1];
                    mean_homloci += homloci * KN->x1[i0][j0][i1][j1];
                    mean_totmuts += totmuts * KN->x1[i0][j1][i1][j1];
                }
            }
        }
    }
    var_hetloci = 0.0;
    var_homloci = 0.0;
    var_totmuts = 0.0;
    for (i0=0; i0 <= KN->MI0; i0++) {
        for (j0=0; j0 <= KN->MJ0; j0++) {
            for (i1=0; i1 <= KN->MI1; i1++) {
                for (j1=0; j1 <= KN->MJ1; j1++) {
                    switch(m) {
                    case 0:
                        hetloci = i0;
                        homloci = j0;
                        break;
                    case 1:
                        hetloci = i1;
                        homloci = j1;
                        break;
                    default:
                        not_implemented(thisfunction, "m > 1");
                        break;
                    }
                    totmuts = hetloci + 2*homloci;
                    t1 = hetloci - mean_hetloci;
                    var_hetloci += (t1 * t1) * KN->x1[i0][j0][i1][j1];
                    t1 = homloci - mean_homloci;
                    var_homloci += (t1 * t1) * KN->x1[i0][j0][i1][j1];
                    t1 = totmuts - mean_totmuts;
                    var_totmuts += (t1 * t1) * KN->x1[i0][j0][i1][j1];
                }
            }
        }
    }
    KStats_n.mean_hetloci[m] = mean_hetloci;
    KStats_n.mean_homloci[m] = mean_homloci;
    KStats_n.mean_totmuts[m] = mean_totmuts;
    KStats_n.var_hetloci[m] = var_hetloci;
    KStats_n.var_homloci[m] = var_homloci;
    KStats_n.var_totmuts[m] = var_totmuts;
}

