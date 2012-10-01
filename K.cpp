#include "K.h"
#include "Trajectory.h"


/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Global variables (keep to an absolute minimum)                */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
int         main_unnested       (int argc, char* argv[])
{
    const char* thisfunction = "main_unnested";
    KConfig K;
    KInt    MI =   500;
    KInt    MJ =   40;
    KInt    g  =   1;
    KScalar U  =   1.0;
    KScalar s  =   1.0;
    KScalar h  =   0.0;
    KScalar S  =   0.01;  /* selfing rate */
    Trajectory trajectory;
    
    if (0) {
        K = initiate_quick(MI, MJ, g, U, s, h, S);
        K->fitness_function = FITNESS_MULTIPLICATIVE;
    } else {
        K = initiate_KConfig();

        initiate_load_classes(K, MI, MJ);
        initiate_genotypes(K, g);
        K->U = U;
        K->fitness_function = FITNESS_MULTIPLICATIVE;
        K->fit_s = s;
        K->fit_h = h;
        K->S[0] = S;
        K->epsilon = 0.00000001;

        if (cmdline_args(K, argc, argv)) {
            fatal("Usage: K --help");
        }

        set_repro(K, K->S[0], 0.0, 0.0, 0.0);
    }

    // set_debug(DEBUG_LETHALS);
    // set_debug(DEBUG_TRACE1);
    // set_debug(DEBUG_TRACE2);
    // set_debug(DEBUG_GENERATIONS);
    // set_debug(DEBUG_FOLLOW_EQUILIBRIUM);
    // set_debug(DEBUG_EQUILIBRIUM);
    // set_debug(DEBUG_NORMALIZATION);
    // set_debug(DEBUG_TRUNCATE);
    // set_debug(DEBUG_TRUNCATE_DETAIL);

    initiate_model_state(K);

    compute_adults_initial(K);

    // K->option_table = 1;
	// K->load_savefile = 1;
    if (K->load_savefile) {
        load_savefile(K, K->x);
    } else {
        fill_KArray(K, K->x, 0.0);
        K->x[10][0][0] = 1.0;
    }

    // fill_KArray(K, K->x, 0.0);
    // K->x[20][0][0] = 1.0;
    /* the above method leaves all adults in the population
    **     in one load class, L(0,0).  Now, we'll shift these
    **     around a little to examine the results of outcrossing */
    /* add a few mutations */
    /*
    for (i=0; i < 10; i++) {
        apply_mutation(K, temp, K->x);
        copy_KArray(K, K->x, temp);
    }
    */

    // printf("K --------------------------------------------------\n");
    // printf("U\ts\th\tS\n%lg\t%lg\t%lg\t%lg\n", K->U, K->fit_s, K->fit_h, K->S[0]);

    while (! is_equilibrium(K)) {

        if (K->generation > GENERATION_CUTOFF) {
			IF_DEBUG(DEBUG_GENERATIONS)
				printf("exceeded GENERATION_CUTOFF=%d, stopping\n", GENERATION_CUTOFF);
			IF_DEBUG(DEBUG_TRACE1)
				printf("exceeded GENERATION_CUTOFF=%d, stopping\n", GENERATION_CUTOFF);
            break;
        }

        IF_DEBUG(DEBUG_GENERATIONS)
			printf("generation %d\n", K->generation);
        IF_DEBUG(DEBUG_TRACE1)
			printf("generation %d\n", K->generation);

        if (K->option_truncate)
            truncate_KArray(K, K->x, LOADCLASS_TRUNCATE);

        compute_mutation(K);
        compute_self_progeny(K);
        compute_apomixis_progeny(K);
        // compute_outcross_progeny(K);  // inefficient Kondrashov method
        compute_gametes(K);
        compute_zygotes(K);
        compute_summed_progeny(K);
        compute_selection(K);
        compute_adults_nextgen(K);

        IF_DEBUG(DEBUG_EQUILIBRIUM) {
            printf("checking equilibrium at end of generation %d------------\n", 
				   K->generation - 1);
            printf("dump of K->x[..][0][0]\n");
            dump_KArray(K, K->x, K->MI, 0, 1);
            printf("calling is_equilibrium(K) just to check...\n");
            (void) is_equilibrium(K);
            printf("end of equilibrium check ---------------\n");
        }

        normalize_KArray(K, K->x);

        /*
        if (K->generation - 1 == 0) {
            trajectory.start(K, "trajectory.txt", "end_of_gen", 5);
            trajectory.debug(true);
        } else if (trajectory.active()) {
            trajectory.check(K);
        }
        */

        // dump_KArray(K, K->x, 100, 0, 0);

    }

    // printf("\n");


    if (!K->option_nolethal && K->is_lethal) {
        /*
        ** We need one more run through everything to get the
        ** progeny arrays built correctly.  We don't do mutation
        ** here.
        */
        K->createlethal = 1;
        IF_DEBUG(DEBUG_LETHALS)
            printf("%s: one more run, to get progeny arrays\n", thisfunction);
        IF_DEBUG(DEBUG_LETHALS)
            printf("%s: K->createlethal=%d\n", thisfunction, K->createlethal);
        compute_self_progeny(K);
        compute_apomixis_progeny(K);
        compute_gametes(K);
        compute_zygotes(K);
        compute_summed_progeny(K);
        compute_selection(K);
        compute_adults_nextgen(K);
        K->generation--;
        normalize_KArray(K, K->x);
    }

    stats_all(K);
    stats_print(K);

    if (trajectory.active() && trajectory.lastgen() != K->generation) {
        trajectory.write(K);
        trajectory.stop();
    }

    // dump_KArray(K, K->x, 10, 0, 0);

    //K->save_savefile = 1;
	if (K->save_savefile) {
        save_savefile(K, K->x);
    }

    return 0;
}

