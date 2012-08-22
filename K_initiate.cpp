#include "K.h"

/*///////////////////////////////////////////////////////////////*/
KConfig     initiate_KConfig    (void)
/*
** allocate zeroed space for KConfig
*/
{
    const char* thisfunction = "initiate_KConfig";
    KConfig K;
    K = (KConfig)calloc((size_t)1, sizeof(struct struct_KConfig));
    return K;
}

/*///////////////////////////////////////////////////////////////*/
void		initiate_load_classes  (KConfig K, KInt MI, KInt MJ)
{
	const char* thisfunction = "initiate_load_classes";
    int err = 0;
    char buf[200];
    sprintf(buf, "%s: ", thisfunction);
    /* check array extent, that they don't go too far */
    if (MI > MAX_MI) {
        strcat(buf, "MI too large, adjust MAX_MI; ");
        err++;
    } else if (MI < 0) {
        strcat(buf, "MI too small, <0; ");
        err++;
	}
    if (MJ > MAX_MJ) {
        strcat(buf, "MJ too large, adjust MAX_MJ; ");
        err++;
    } else if (MJ < 0) {
        strcat(buf, "MJ too small, <0; ");
        err++;
    }
    if (MI < MJ) {
        strcat(buf, "MI < MJ; ");
        err++;
    }
    if (err)
        fatal(buf);
	K->MI = MI;
	K->MJ = MJ;
}

/*///////////////////////////////////////////////////////////////*/
void        initiate_genotypes  (KConfig K, KInt g)
{
    const char* thisfunction = "initiate_genotypes";
    if (g != 1) {
        not_implemented(thisfunction, "genotypes != 1");
	}
    K->genotypes = g;
    set_transform_one_genotype(K);
}

/*///////////////////////////////////////////////////////////////*/
KConfig     initiate_quick      (KInt MI, KInt MJ, KInt g,
                                 KScalar U,
                                 KScalar s, KScalar h,
								 KScalar S)
{
    const char* thisfunction = "initiate_quick";
	/* for now, only support one genotype in initiate_quick */
    KConfig K;
    K = initiate_KConfig();
	initiate_load_classes(K, MI, MJ);
	initiate_genotypes(K, g);
    K->U = U;
    K->fit_s = s;
    K->fit_h = h;
	set_repro(K, S, 0.0, 0.0, 0.0);
    return K;
}

/*///////////////////////////////////////////////////////////////*/
void        initiate_model_state	(KConfig K)
{
	const char* thisfunction = "initiate_model_state";
    initiate_mut_term(K);
    initiate_fitness_precomputed(K);
    K->generation = 0;
    if (!K->option_nolethal) {
        if (K->fit_s == 1.0) {
            IF_DEBUG(DEBUG_LETHALS)
                printf("%s: mutation class is lethal, so K->is_lethal=1\n",
                       thisfunction);
            K->is_lethal = 1;
            K->createlethal = 0;
        }
        IF_DEBUG(DEBUG_LETHALS)
            printf("%s: K->is_lethal=%d\n", thisfunction, K->is_lethal);
        IF_DEBUG(DEBUG_LETHALS)
            printf("%s: if this is non-zero, expect normalization problems\n",
                   thisfunction);
    }
}

/*///////////////////////////////////////////////////////////////*/
void        compute_adults_initial (KConfig K)
/*
** Initiate adult frequencies
*/
{
    const char* thisfunction = "compute_adults_initial";
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    if (K->genotypes != 1) {
        not_implemented(thisfunction, "genotypes != 1");
	}
    fill_KArray(K, K->x, 0.0);
    K->x[0][0][0] = 1.0;
}

/*///////////////////////////////////////////////////////////////*/
/* void        compute_adults_initial  (KConfig K) */
/*
** Initiate adult frequencies according to requested genotype
** frequencies
*/
/*
{
    const char* thisfunction = "compute_adults_initial";
    KInt i, j, g;
    not_implemented(thisfunction);
    IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
    for (i=0; i <= K->MI; i++) {
        for (j=0; j <= K->MJ; j++) {
            for (g=0; g < K->genotypes; g++) {
                IF_DEBUG(DEBUG_TRACE2) {
                    printf("%s: [%d][%d][%d]\n",
						   thisfunction, i, j, g);
                }
                if (i == 0 && j == 0) {
                    K->x[i][j][g] = K->initial_genotype_freqs[g];
                } else {
                    K->x[i][j][g] = 0;
                }
            }
        }
    }
    K->generation = 0;
    check_normalization(K, K->x, thisfunction, "K->x");
}
*/

