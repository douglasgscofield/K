#include "K.h"

/*
** TODO: test without mutation, without selection, just looking
**       for the equilibrium distribution of load classes
*/

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Global variables (keep to an absolute minimum)                */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////// */
int
main_nested(
	int argc,
	char *argv[])
{
	const char *thisfunction = "main_nested";
	KConfig_n KN;
	KInt MI0 = 50;
	KInt MJ0 = 10;
	KInt MI1 = 50;
	KInt MJ1 = 10;
	KScalar_n U = { 1.0, 0.0 };
	KScalar_n s = { 1.0, 0.0 };
	KScalar_n h = { 0.0, 0.0 };
	KScalar S = 0.01;			/* selfing rate */
	KScalar A = 0.0;			/* apomixis rate */

	KN = initiate_KConfig_n();

	initiate_load_classes_n(KN, MI0, MJ0, MI1, MJ1);
	KN->U[0] = U[0];
	KN->U[1] = U[1];
	KN->fitness_function[0] = FITNESS_MULTIPLICATIVE_n;
	KN->fitness_function[1] = FITNESS_MULTIPLICATIVE_n;
	KN->fit_s[0] = s[0];
	KN->fit_s[1] = s[1];
	KN->fit_h[0] = h[0];
	KN->fit_h[1] = h[1];
	KN->S = S;
	KN->A = A;
	KN->epsilon = 0.00000001;

	if (cmdline_args_n(KN, argc, argv)) {
		fatal("Usage: K --help");
	}

	set_repro_n(KN, KN->S, KN->A);

	/* set_debug(DEBUG_LETHALS); /* */
	/* set_debug(DEBUG_GENERATIONS); /* */
	/* set_debug(DEBUG_TRACE1); /* */
	/* set_debug(DEBUG_TRACE2); /* */
	/* set_debug(DEBUG_EQUILIBRIUM); /* */
	/* set_debug(DEBUG_NORMALIZATION); /* */
	/* set_debug(DEBUG_TRUNCATE); /* */

	initiate_model_state_n(KN);

	compute_adults_initial_n(KN);

	if (KN->load_savefile) {
		if (KN->current_x != KN_CURRENT_X1)
			fatal("load_savefile code not cognizant of current x");
		load_savefile_n(KN, KN->x1);
	} else {
		fill_KArray_n(KN, KN->x1, 0.0);
		KN->x1[6][0][2][0] = 1.0;
	}

	/*
	 ** fill_KArray_n(KN, KN->x1, 0.0);
	 ** KN->x1[6][0][2][0] = 1.0;
	 */

	/* the above method leaves all adults in the population
	 **     in one load class, L(0,0).  Now, we'll shift these
	 **     around a little to examine the results of outcrossing */
	/* add a few mutations */
	/*
	   for (i=0; i < 10; i++) {
	   apply_mutation(K, temp, K->x);
	   copy_KArray(K, K->x, temp);
	   }
	   /* */
	/*
	   printf("K --------------------------------------------------\n");
	   printf("U\ts\th\tS\n%lg\t%lg\t%lg\t%lg\n", K->U, K->fit_s, K->fit_h, K->S[0]);
	   /* */
	while (!is_equilibrium_n(KN)) {
		if (KN->generation > GENERATION_CUTOFF_n) {
			IF_DEBUG(DEBUG_TRACE1)
				printf("exceeded GENERATION_CUTOFF_n=%d, stopping\n",
					   GENERATION_CUTOFF_n);
			IF_DEBUG(DEBUG_GENERATIONS)
				printf("exceeded GENERATION_CUTOFF_n=%d, stopping\n",
					   GENERATION_CUTOFF_n);
			break;
		}
		IF_DEBUG(DEBUG_TRACE1)
			printf("generation %d\n", KN->generation);
		IF_DEBUG(DEBUG_GENERATIONS)
			printf("generation %d\n", KN->generation);
		compute_adults_prevgen_n(KN);
		if (KN->option_truncate)
			truncate_KArray_n(KN, KN->x1, LOADCLASS_TRUNCATE);
		compute_mutation_n(KN);
		compute_self_progeny_n(KN);
		compute_apomixis_progeny_n(KN);
		compute_gametes_n(KN);
		compute_zygotes_n(KN);

		compute_summed_progeny_n(KN);
		compute_selection_n(KN);
		compute_adults_nextgen_n(KN);

		if (KN->current_x != KN_CURRENT_X1) {
			char buf[200];
			sprintf(buf, "%s: wrong current x array = %d", thisfunction,
					KN->current_x);
			fatal(buf);
		}

		IF_DEBUG(DEBUG_EQUILIBRIUM) {
			printf
				("checking equilibrium at end of generation %d------------\n",
				 KN->generation - 1);
			printf("dump of KN->x[..][0][..][0]\n");
			dump_KArray_n(KN, KN->x1, KN->MI0, 0, KN->MI1, 0);
			printf("calling is_equilibrium(K) just to check...\n");
			(void) is_equilibrium_n(KN);
			printf("end of equilibrium check ---------------\n");
		}
		 /**/ normalize_KArray_n(KN, KN->x1);
		 /**/
			/*
			   dump_KArray_n(KN, KN->x, 50, 0, 50, 0);
			   /* */
	 /**/}
	/*
	   printf("\n");
	   /* */

	if (KN->is_lethal[0]) {
		/*
		 ** Because of the way stats are gathered in the two-class
		 ** model, it's sufficient just to allow the creation of
		 ** homozygous lethal mutant progeny, by setting
		 ** KN->createlethal[].  The stats routines
		 ** will reconstitute all the progeny classes for us.
		 */
		KN->createlethal[0] = 1;
		IF_DEBUG(DEBUG_LETHALS)
			printf("%s: KN->createlethal[0]=%d\n", thisfunction,
				   KN->createlethal[0]);
	}
	if (KN->is_lethal[1]) {
		KN->createlethal[1] = 1;
		IF_DEBUG(DEBUG_LETHALS)
			printf("%s: KN->createlethal[1]=%d\n", thisfunction,
				   KN->createlethal[1]);
	}

	/**/ stats_all_n(KN);
	stats_print_n(KN);
	 /**/
		/*
		   dump_KArray_n(KN, KN->x1, 50, 0, 50, 0);
		   /* */
		if (KN->save_savefile) {
		save_savefile_n(KN, KN->x1);
	}

	return 0;
}
