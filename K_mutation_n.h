/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Methods for handling mutation in nested mutation models       */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

void compute_mutation_n(
	KConfig_n KN);
void apply_mutation_n(
	KConfig_n K,
	KArray_n & to,
	KArray_n & from);
void apply_mutation_general_n(
	KConfig K,
	KScalar U,
	KArray & to,
	KArray & from);
void initiate_mut_term_n(
	KConfig_n KN);
KScalar mut_term_n(
	KConfig_n KN,
	KInt x0,
	KInt x1);
KScalar mut_term_computed_n(
	KConfig_n KN,
	KInt x0,
	KInt x1);
KScalar mut_term_general_n(
	KScalar U,
	KInt x);
