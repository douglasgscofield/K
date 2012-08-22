#include "K_n.h"


/*///////////////////////////////////////////////////////////////*/
void        compute_adults_prevgen_n(KConfig_n KN)
/*
** Store away adult values from previous generation
*/
{
    const char* thisfunction = "compute_adults_prevgen_n";
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (KN->current_x != KN_CURRENT_X1) {
        char buf[200];
        sprintf(buf, "%s: wrong current x array = %d",
                thisfunction, KN->current_x);
        fatal(buf);
    }
    apply_adults_prevgen_n(KN, KN->x_prevgen, KN->x1);
}

/*///////////////////////////////////////////////////////////////*/
void        apply_adults_prevgen_n  (KConfig_n KN,
                                     KArray_n& prevgen,
                                     KArray_n& from)
/*
** Create reproductive adults that will begin the next generation
** and store adults from this generation in prevgen
*/
{
    const char* thisfunction = "apply_adults_prevgen_n";
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    copy_KArray_n(KN, prevgen, from);
    check_normalization_n(KN, prevgen, thisfunction, "prevgen");
}

