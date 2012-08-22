#include "K.h"

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Fitness model operations - for nested mutation classes        */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* File-scope types and variables */

struct struct_KFitness_n fitness_functions_n[fitness_function_VALS_n + 1] = {
	/* The order of the entries in this array should exactly match
	 ** the values of the corresponding FITNESS_ value #define'd in
	 ** the file K_selection_n.h.
	 */
/* 0 */ {
		 FITNESS_INVALID_n,		/* thisnum */
		 (KPtr_ff_n) 0,			/* func */
		 "-invalid nested-",	/* name */
		 0						/* mustcompute */
		 },

/* 1 */ {
		 FITNESS_MULTIPLICATIVE_n,
		 fitness_multiplicative_n,
		 "multiplicative nested",
		 0},

/* 2 */ {
		 FITNESS_KONDRASHOV_n,
		 fitness_kondrashov_n,
		 "kondrashov nested",
		 0}
};

/* MAX_fitness_function should be set to the index of the last
** entry in the fitness_functions[] table.
*/
static KInt MAX_fitness_function_n = FITNESS_KONDRASHOV_n;


/* ///////////////////////////////////////////////////////////// */
void
compute_selection_n(
	KConfig_n KN)
/*
** K->X is correct upon exit.
*/
{
	const char *thisfunction = "compute_selection_n";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	if (KN->current_x != KN_CURRENT_X1) {
		char buf[200];
		sprintf(buf, "%s: wrong current x array = %d", thisfunction,
				KN->current_x);
		fatal(buf);
	}
	apply_selection_n(KN, KN->x2, KN->x1);
	KN->current_x = KN_CURRENT_X2;
}

/* ///////////////////////////////////////////////////////////// */
void
apply_selection_n(
	KConfig_n KN,
	KArray_n & to,
	KArray_n & from)
/*
** 'to' is correct upon exit.
*/
{
	const char *thisfunction = "apply_selection_n";
	KInt i0, j0, i1, j1;
	KScalar recip_w_mean;
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	recip_w_mean = 1.0 / mean_fitness_allprogeny_n(KN, from);
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					IF_DEBUG(DEBUG_TRACE2)
						if (!(i0 % 10) && !(j0 % 10) && !(i1 % 10) &&
							!(j1 % 10))
						printf("sel[%d,%d,%d,%d] ", i0, j0, i1, j1);
					to[i0][j0][i1][j1] =
						from[i0][j0][i1][j1] * fitness_n(KN, i0, j0, i1,
														 j1) * recip_w_mean;
				}
			}
		}
	}
	IF_DEBUG(DEBUG_TRACE2) printf("\n");
	check_normalization_n(KN, to, thisfunction, "to");
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mfitness_computed_n(
	KConfig_n KN,
	KMutClass m,
	KInt i,
	KInt j)
/*
** Returns the fitness of the mut class m load class L(i,j) as
** determined by KN->fitness_function[m].
*/
{
	const char *thisfunction = "mfitness_computed_n";
	KScalar w;
	if (KN->fitness_function[m] == FITNESS_INVALID ||
		KN->fitness_function[m] < 0 ||
		KN->fitness_function[m] > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	/* Note that the genotype passed to the fitness function is
	 ** always 0 for nested mutation classes.
	 */
	w = (*(get_fitness_function_n(KN->fitness_function[m]))) (KN, m, i, j);
	return w;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mfitness_n(
	KConfig_n KN,
	KMutClass m,
	KInt i,
	KInt j)
{
	const char *thisfunction = "mfitness_n";
	if (fitness_functions[KN->fitness_function[m]].mustcompute) {
		return mfitness_computed_n(KN, m, i, j);
	} else {
		return KN->fitness_precomputed[m][i][j];
	}
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_from_mfitness(
	KConfig_n KN,
	KScalar_n mw)
/*
** This is mostly a place-holder function until I can implement
** a real user-defined mechanism for combining the mut class 
** fitnesses
*/
{
	KMutClass m;
	KScalar w = 1.0;
	for (m = 0; m < KN->mutclasses; m++) {
		w *= mw[m];
	}
	return w;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_computed_n(
	KConfig_n KN,
	KInt i0,
	KInt j0,
	KInt i1,
	KInt j1)
{
	const char *thisfunction = "fitness_computed_n";
	KMutClass m;
	KScalar w;
	KScalar_n mw;
	m = 0;
	mw[m] = mfitness_computed_n(KN, m, i0, j0);
	m = 1;
	mw[m] = mfitness_computed_n(KN, m, i1, j1);
	w = fitness_from_mfitness(KN, mw);
	return w;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_n(
	KConfig_n KN,
	KInt i0,
	KInt j0,
	KInt i1,
	KInt j1)
{
	const char *thisfunction = "fitness_n";
	KScalar w;
	KScalar_n mw;
	KMutClass m;
	m = 0;
	if (fitness_functions[KN->fitness_function[m]].mustcompute) {
		mw[m] = mfitness_computed_n(KN, m, i0, j0);
	} else {
		mw[m] = KN->fitness_precomputed[m][i0][j0];
	}
	m = 1;
	if (fitness_functions[KN->fitness_function[m]].mustcompute) {
		mw[m] = mfitness_computed_n(KN, m, i1, j1);
	} else {
		mw[m] = KN->fitness_precomputed[m][i1][j1];
	}
	w = fitness_from_mfitness(KN, mw);
	return w;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mean_fitness_n(
	KConfig_n KN,
	KArray_n & a)
/*
** The mean fitness of the load class-genotype array
*/
{
	const char *thisfunction = "mean_fitness_n";
	KScalar sum;
	sum = sum_KArray_n(KN, a);
	if (sum == 0.0) {
		return 0.0;
	}
	return cumulative_fitness_n(KN, a) / sum;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mean_fitness_allprogeny_n(
	KConfig_n KN,
	KArray_n & a)
/*
** The mean fitness of the load class-genotype array
*/
{
	const char *thisfunction = "mean_fitness_n";
	KScalar sum;
	sum = sum_KArray_n(KN, a);
	if (sum == 0.0) {
		return 0.0;
	}
	IF_DEBUG(DEBUG_LETHALS) printf("%s: sum_KArray_n was=%g\n", thisfunction,
								   sum);
	IF_DEBUG(DEBUG_LETHALS) printf("%s: KN->createlethal[0]=%d\n",
								   thisfunction, KN->createlethal[0]);
	IF_DEBUG(DEBUG_LETHALS) printf("%s: KN->createlethal[1]=%d\n",
								   thisfunction, KN->createlethal[1]);
	if ((KN->is_lethal[0] || KN->is_lethal[1]) && sum < 1.0)
		/*
		 ** see mean_fitness_allprogeny()
		 */
		sum = 1.0;
	IF_DEBUG(DEBUG_LETHALS) printf("%s: sum_KArray_n is=%g\n", thisfunction,
								   sum);
	return cumulative_fitness_n(KN, a) / sum;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
cumulative_fitness_n(
	KConfig_n KN,
	KArray_n & a)
/*
** The cumulative fitness of the load class-genotype array
*/
{
	const char *thisfunction = "cumulative_fitness_n";
	KInt i0, j0, i1, j1;
	KScalar ans = 0.0;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					ans += a[i0][j0][i1][j1] * fitness_n(KN, i0, j0, i1, j1);
				}
			}
		}
	}
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mean_fitness_outcross_n(
	KConfig_n KN,
	KVector_n & out)
/*
** The mean fitness of outcross progeny, which are represented
** by a KVector_n.
*/
{
	const char *thisfunction = "mean_fitness_outcross_n";
	KScalar sum;
	sum = sum_KVector_n(KN, out);
	if (sum == 0.0) {
		return 0.0;
	}
	return cumulative_fitness_outcross_n(KN, out) / sum;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
cumulative_fitness_outcross_n(
	KConfig_n KN,
	KVector_n & out)
/*
** The cumulative fitness of outcross progeny, which are represented
** by a KVector_n.
*/
{
	const char *thisfunction = "cumulative_fitness_outcross_n";
	KInt i0, i1;
	KScalar ans = 0.0;
	for (i0 = 0; i0 <= KN->MI1; i0++) {
		for (i1 = 0; i1 <= KN->MI1; i1++) {
			/* note that for fitness of outcrossed progeny,
			 ** the number of heterozygous mutations is
			 ** always 0
			 */
			ans += out[i0][i1] * fitness_n(KN, i0, 0, i1, 0);
		}
	}
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
void
initiate_fitness_precomputed_n(
	KConfig_n KN)
/*
** Compute precomputed fitness array.  Note that this computes
** m self-contained arrays, one for each mutation class.
*/
{
	const char *thisfunction = "initiate_fitness_precomputed_n";
	KMutClass m;
	KInt i, j;
	m = 0;
	for (i = 0; i <= KN->MI0; i++) {
		for (j = 0; j <= KN->MJ0; j++) {
			KN->fitness_precomputed[m][i][j] =
				mfitness_computed_n(KN, m, i, j);
		}
	}
	m = 1;
	for (i = 0; i <= KN->MI1; i++) {
		for (j = 0; j <= KN->MJ1; j++) {
			KN->fitness_precomputed[m][i][j] =
				mfitness_computed_n(KN, m, i, j);
		}
	}
}

/* ///////////////////////////////////////////////////////////// */
KInt
register_fitness_function_n(
	KPtr_ff_n pff_n,
	char *ff_name)
/*
** Registers the fitness function specified by pff_n and returns 
** the number of the fitness model's registration slot.
*/
{
	const char *thisfunction = "register_fitness_function_n";
	++MAX_fitness_function_n;
	if (MAX_fitness_function_n > fitness_function_VALS_n) {
		char buf[200];
		sprintf(buf, "%s: exceeded fitness_function_VALS_n", thisfunction);
		fatal(buf);
	}
	fitness_functions_n[MAX_fitness_function_n].thisnum =
		MAX_fitness_function_n;
	fitness_functions_n[MAX_fitness_function_n].func = pff_n;
	fitness_functions_n[MAX_fitness_function_n].name = ff_name;
	fitness_functions_n[MAX_fitness_function_n].mustcompute = 0;
	return MAX_fitness_function_n;
}

/* ///////////////////////////////////////////////////////////// */
void
must_compute_fitness_function_n(
	KInt ff)
/*
** Specifies that the fitness model must be computed each time
** that fitness is required; we cannot use values stored in
** K->fitness_precomputed[][][]
*/
{
	const char *thisfunction = "must_compute_fitness_function_n";
	if (ff <= 0 || ff > fitness_function_VALS_n) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	fitness_functions_n[ff].mustcompute = 1;
}

/* ///////////////////////////////////////////////////////////// */
KPtr_ff_n
get_fitness_function_n(
	KInt ff)
/*
** Returns a pointer to the function corresponding to the
** fitness model specified by ff
*/
{
	const char *thisfunction = "get_fitness_function_n";
	if (ff <= 0 || ff > fitness_function_VALS_n) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	return fitness_functions_n[ff].func;
}

/* ///////////////////////////////////////////////////////////// */
char *
get_fitness_function_name_n(
	KInt ff)
/*
** Returns a pointer to the function corresponding to the
** fitness model specified by ff
*/
{
	const char *thisfunction = "get_fitness_function_name_n";
	if (ff <= 0 || ff > fitness_function_VALS_n) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	return fitness_functions_n[ff].name;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_multiplicative_n(
	KConfig_n KN,
	KMutClass m,
	KInt i,
	KInt j)
/*
** Compute fitness based on load class using multiplicative 
** fitness model a la Charlesworth et al (1990)
*/
{
	const char *thisfunction = "fitness_multiplicative_n";
	KScalar t1, t2, ans;
	t1 = pow((1 - KN->fit_h[m] * KN->fit_s[m]), (KScalar) i);
	t2 = pow((1 - KN->fit_s[m]), (KScalar) j);
	ans = t1 * t2;
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_kondrashov_n(
	KConfig_n KN,
	KMutClass m,
	KInt i,
	KInt j)
/* 
** Compute fitness based on load class using Kondrashov (1985)
*/
{
	const char *thisfunction = "fitness_kondrashov_n";
	KScalar t1, t2, ans;
	t1 = ((KScalar) i + (KN->fit_d[m] * (KScalar) j)) / KN->fit_k[m];
	t2 = pow(t1, (KScalar) KN->fit_alpha[m]);
	ans = 1 - t2;
	return ans;
}
