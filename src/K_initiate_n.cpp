#include "K_n.h"

/*///////////////////////////////////////////////////////////////*/
KConfig_n   initiate_KConfig_n  (void)
/*
** allocate zeroed space for KConfig_n
*/
{
    const char* thisfunction = "initiate_KConfig_n";
    KConfig_n KN;
    KN = (KConfig_n)calloc((size_t)1, sizeof(struct struct_KConfig_n));
    return KN;
}

/*///////////////////////////////////////////////////////////////*/
void        initiate_load_classes_n (KConfig_n KN, 
                                     KInt MI0, KInt MJ0,
                                     KInt MI1, KInt MJ1)
{
    const char* thisfunction = "initiate_load_classes_n";
    int err = 0;
    char buf[200];
    KN->mutclasses = 2;
    sprintf(buf, "%s: ", thisfunction);
    /* check array extent, that they don't go too far */
    if (MI0 > MAX_MI_n) {
        strcat(buf, "MI0 too large, adjust MAX_MI_n; "); err++;
    } else if (MI0 < 0) {
        strcat(buf, "MI0 too small, <0; "); err++;
    }
    if (MJ0 > MAX_MJ_n) {
        strcat(buf, "MJ0 too large, adjust MAX_MJ_n; "); err++;
    } else if (MJ0 < 0) {
        strcat(buf, "MJ0 too small, <0; "); err++;
    }
    if (MI0 < MJ0) {
        strcat(buf, "MI0 < MJ0; "); err++;
    }
    if (MI1 > MAX_MI_n) {
        strcat(buf, "MI1 too large, adjust MAX_MI_n; "); err++;
    } else if (MI1 < 0) {
        strcat(buf, "MI1 too small, <0; "); err++;
    }
    if (MJ1 > MAX_MJ_n) {
        strcat(buf, "MJ1 too large, adjust MAX_MJ_n; "); err++;
    } else if (MJ1 < 0) {
        strcat(buf, "MJ1 too small, <0; "); err++;
    }
    if (MI1 < MJ1) {
        strcat(buf, "MI1 < MJ1; "); err++;
    }
    if (err)
        fatal(buf);
    KN->MI0 = MI0;
    KN->MJ0 = MJ0;
    KN->MI1 = MI1;
    KN->MJ1 = MJ1;
}

/*///////////////////////////////////////////////////////////////*/
void        initiate_model_state_n  (KConfig_n KN)
{
    const char* thisfunction = "initiate_model_state_n";
    initiate_mut_term_n(KN);
    initiate_fitness_precomputed_n(KN);
    KN->generation = 0;
    if (!KN->option_nolethal) {
        if (KN->fit_s[0] == 1.0) {
            IF_DEBUG(DEBUG_LETHALS)
                fprintf(stderr, "%s: mutation class 0 is lethal, so KN->is_lethal[0]=1\n",
                       thisfunction);
            KN->is_lethal[0] = 1;
            KN->createlethal[0] = 0;
        }
        if (KN->fit_s[1] == 1.0) {
            IF_DEBUG(DEBUG_LETHALS)
                fprintf(stderr, "%s: mutation class 1 is lethal, so KN->is_lethal[1]=1\n",
                       thisfunction);
            KN->is_lethal[1] = 1;
            KN->createlethal[1] = 0;
        }
        IF_DEBUG(DEBUG_LETHALS)
            fprintf(stderr, "%s: KN->is_lethal[0]=%d\n", thisfunction, KN->is_lethal[0]);
        IF_DEBUG(DEBUG_LETHALS)
            fprintf(stderr, "%s: KN->is_lethal[1]=%d\n", thisfunction, KN->is_lethal[1]);
        IF_DEBUG(DEBUG_LETHALS)
            fprintf(stderr, "%s: if either of these are non-zero, expect normalization problems\n",
                   thisfunction);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_adults_initial_n    (KConfig_n KN)
/*
** Initiate adult frequencies to KN->x1
*/
{
    const char* thisfunction = "compute_adults_initial_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    apply_adults_initial_n(KN, KN->x1);
    KN->current_x = KN_CURRENT_X1;
}

/*///////////////////////////////////////////////////////////////*/
void        apply_adults_initial_n  (KConfig_n KN, KArray_n& a)
/*
** Initiate adult frequencies
*/
{
    const char* thisfunction = "apply_adults_initial_n";
    IF_DEBUG(DEBUG_TRACE1) fprintf(stderr, "%s\n", thisfunction);
    fill_KArray_n(KN, a, 0.0);
    a[0][0][0][0] = 1.0;
    check_normalization_n(KN, a, thisfunction, "a");
}

