#include "K.h"

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Fitness model operations                                      */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* File-scope types and variables */

struct struct_KFitness fitness_functions[fitness_function_VALS + 1] = {
	/* The order of the entries in this array should exactly match
	 ** the values of the corresponding FITNESS_ value #define'd in
	 ** the file K_fitness.h.
	 */
/* 0 */ {
		 FITNESS_INVALID,		/* thisnum */
		 (KPtr_ff) 0,			/* func */
		 (char *) "-invalid-",	/* name */
		 0						/* mustcompute */
		 },

/* 1 */ {
		 FITNESS_MULTIPLICATIVE,
		 fitness_multiplicative,
		 (char *) "multiplicative",
		 0},

/* 2 */ {
		 FITNESS_KONDRASHOV,
		 fitness_kondrashov,
		 (char *) "kondrashov",
		 0}
};

/* MAX_fitness_function should be set to the index of the last
** entry in the fitness_functions[] table.
*/
static KInt MAX_fitness_function = FITNESS_KONDRASHOV;


/* ///////////////////////////////////////////////////////////// */
void
compute_selection(
	KConfig K)
/*
** K->X is correct upon exit.
*/
{
	const char *thisfunction = "compute_selection";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	apply_selection(K, K->X, K->xpp);
	check_normalization(K, K->X, thisfunction, "K->X");
}

/* ///////////////////////////////////////////////////////////// */
void
apply_selection(
	KConfig K,
	KArray & to,
	KArray & from)
/*
** to is correct upon exit.
*/
{
	const char *thisfunction = "apply_selection";
	KInt g, i, j;
	KScalar recip_w_mean;
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	recip_w_mean = 1.0 / mean_fitness_allprogeny(K, from);
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			for (g = 0; g < K->genotypes; g++) {
				to[i][j][g] =
					from[i][j][g] * fitness(K, i, j, g) * recip_w_mean;
			}
		}
	}
	check_normalization(K, to, thisfunction, "to");
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_computed(
	KConfig K,
	KInt i,
	KInt j,
	KInt g)
/*
** Returns the fitness as a function of i, j, g using the 
** model specified in K->fitness_function.
*/
{
	const char *thisfunction = "fitness_computed";
	KScalar w;
	if (K->fitness_function == FITNESS_INVALID || K->fitness_function < 0 ||
		K->fitness_function > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	w = (*(get_fitness_function(K->fitness_function))) (K, i, j, g);
	return w;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness(
	KConfig K,
	KInt i,
	KInt j,
	KInt g)
/*
** Returns the fitness as a function of i, j, g using the 
** fitness function specified in K->fitness_function.  If 
** fitness is precomputed, then all the corresponding 
** fitnesses are already in K->fitness_precomputed[][][]
*/
{
	const char *thisfunction = "fitness";
	if (fitness_functions[K->fitness_function].mustcompute) {
		return fitness_computed(K, i, j, g);
	} else {
		return K->fitness_precomputed[i][j][g];
	}
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mean_fitness(
	KConfig K,
	KArray & a)
/*
** The mean fitness of the load class-genotype array
*/
{
	const char *thisfunction = "mean_fitness";
	KScalar sum;
	sum = sum_KArray(K, a);
	if (sum == 0.0) {
		return 0.0;
	}
	return cumulative_fitness(K, a) / sum;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
mean_fitness_allprogeny(
	KConfig K,
	KArray & a)
/*
** The mean fitness of the load class-genotype array
*/
{
	const char *thisfunction = "mean_fitness_allprogeny";
	KScalar sum;
	sum = sum_KArray(K, a);
	if (sum == 0.0) {
		return 0.0;
	}
	IF_DEBUG(DEBUG_LETHALS) printf("%s: sum_KArray was=%g\n", thisfunction,
								   sum);
	IF_DEBUG(DEBUG_LETHALS) printf("%s: K->createlethal=%d\n", thisfunction,
								   K->createlethal);
	if (K->is_lethal && sum < 1.0)
		/*
		 ** sum should be 1.0 in any case, its computation here is
		 ** done to deal with normalization issues that may arise
		 ** during computation.  Hopefully, taking this approach
		 ** with K->is_lethal won't affect the approach to equilibrium.
		 ** It is necessary with K->is_lethal because we're not producing
		 ** the selfed progeny in classes j>0, so that screws up the
		 ** proportion of total progeny so that with S>0, sum<1.0.
		 */
		sum = 1.0;
	IF_DEBUG(DEBUG_LETHALS) printf("%s: sum_KArray is=%g\n", thisfunction,
								   sum);
	return cumulative_fitness(K, a) / sum;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
cumulative_fitness(
	KConfig K,
	KArray & a)
/*
** The cumulative fitness of the load class-genotype array
*/
{
	const char *thisfunction = "cumulative_fitness";
	KInt i, j, g;
	KScalar ans = 0.0;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			for (g = 0; g < K->genotypes; g++) {
				ans += (a[i][j][g] * fitness(K, i, j, g));
			}
		}
	}
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
void
initiate_fitness_precomputed(
	KConfig K)
/*
** Compute fitness array.
** The fitness function is fixed into K->fitness at model
** initiation.  K->fitness_precomputed is correct upon exit.
*/
{
	const char *thisfunction = "initiate_fitness_precomputed";
	KInt i, j, g;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			for (g = 0; g < K->genotypes; g++) {
				K->fitness_precomputed[i][j][g] = fitness_computed(K, i, j, g);
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////// */
KInt
register_fitness_function(
	KPtr_ff pff,
	char *ff_name)
/*
** Registers the fitness function specified by pfm and returns 
** the number of the fitness model's registration slot.
*/
{
	const char *thisfunction = "register_fitness_function";
	++MAX_fitness_function;
	if (MAX_fitness_function > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: exceeded fitness_function_VALS", thisfunction);
		fatal(buf);
	}
	fitness_functions[MAX_fitness_function].thisnum = MAX_fitness_function;
	fitness_functions[MAX_fitness_function].func = pff;
	fitness_functions[MAX_fitness_function].name = ff_name;
	fitness_functions[MAX_fitness_function].mustcompute = 0;
	return MAX_fitness_function;
}

/* ///////////////////////////////////////////////////////////// */
void
must_compute_fitness_function(
	KInt ff)
/*
** Specifies that the fitness model must be computed each time
** that fitness is required; we cannot use values stored in
** K->fitness_precomputed[][][]
*/
{
	const char *thisfunction = "must_compute_fitness_function";
	if (ff <= 0 || ff > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	fitness_functions[ff].mustcompute = 1;
}

/* ///////////////////////////////////////////////////////////// */
KPtr_ff
get_fitness_function(
	KInt ff)
/*
** Returns a pointer to the function corresponding to the
** fitness model specified by ff
*/
{
	const char *thisfunction = "get_fitness_function";
	if (ff <= 0 || ff > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	return fitness_functions[ff].func;
}

/* ///////////////////////////////////////////////////////////// */
char *
get_fitness_function_name(
	KInt ff)
/*
** Returns a pointer to the function corresponding to the
** fitness model specified by ff
*/
{
	const char *thisfunction = "get_fitness_function_name";
	if (ff <= 0 || ff > fitness_function_VALS) {
		char buf[200];
		sprintf(buf, "%s: model not set or out of range", thisfunction);
		fatal(buf);
	}
	return fitness_functions[ff].name;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_multiplicative(
	KConfig K,
	KInt i,
	KInt j,
	KInt g)
/*
** Compute fitness based on load class using multiplicative 
** fitness model a la Charlesworth et al (1991)
*/
{
	const char *thisfunction = "fitness_multiplicative";
	KScalar t1, t2, ans;
	t1 = pow((1 - K->fit_h * K->fit_s), (KScalar) i);
	t2 = pow((1 - K->fit_s), (KScalar) j);
	ans = t1 * t2;
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
fitness_kondrashov(
	KConfig K,
	KInt i,
	KInt j,
	KInt g)
/* 
** Compute fitness based on load class using Kondrashov (1985)
*/
{
	const char *thisfunction = "fitness_kondrashov";
	KScalar t1, t2, ans;
	t1 = ((KScalar) i + (K->fit_d * (KScalar) j)) / K->fit_k;
	t2 = pow(t1, (KScalar) K->fit_alpha);
	ans = 1 - t2;
	return ans;
}
