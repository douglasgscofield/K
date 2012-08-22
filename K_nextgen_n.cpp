#include "K.h"


/* ///////////////////////////////////////////////////////////// */
void
compute_adults_nextgen_n(
	KConfig_n KN)
/*
** Create reproductive adults that will begin the next generation
** K->x1 is correct and ready for the next iteration upon exit.
*/
{
	const char *thisfunction = "compute_adults_nextgen_n";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	if (KN->current_x != KN_CURRENT_X2) {
		char buf[200];
		sprintf(buf, "%s: wrong current x array = %d", thisfunction,
				KN->current_x);
		fatal(buf);
	}
	apply_adults_nextgen_n(KN, KN->x1, KN->x2);
	KN->current_x = KN_CURRENT_X1;
}

/* ///////////////////////////////////////////////////////////// */
void
apply_adults_nextgen_n(
	KConfig_n KN,
	KArray_n & nextgen,
	KArray_n & from)
/*
** Create reproductive adults that will begin the next generation
** and store adults from this generation in prevgen
*/
{
	const char *thisfunction = "apply_adults_nextgen_n";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	copy_KArray_n(KN, nextgen, from);
	KN->generation++;
	check_normalization_n(KN, nextgen, thisfunction, "nextgen");
}
