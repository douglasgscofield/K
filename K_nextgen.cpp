#include "K.h"


/* ///////////////////////////////////////////////////////////// */
void
compute_adults_nextgen(
	KConfig K)
/*
** Create reproductive adults that will begin the next generation
** K->x is correct and ready for the next iteration upon exit.
*/
{
	const char *thisfunction = "compute_adults_nextgen";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	copy_KArray(K, K->x_prevgen, K->x);
	copy_KArray(K, K->x, K->X);
	K->generation++;
	check_normalization(K, K->x_prevgen, thisfunction, "K->x_prevgen");
	check_normalization(K, K->x, thisfunction, "K->x");
}
