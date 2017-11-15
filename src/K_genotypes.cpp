#include "K.h"

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Functions for handling the non-mutable genotypes              */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
void        set_transform_one_genotype  (KConfig K)
/*
** These arrays assume one genotype that does not change
*/
{
    const char* thisfunction = "set_transform_one_genotype";
    KInt xi, zeta, g;
    /* O_transform holds proportion of g-th genotype appearing 
    ** in progeny of outcrossed mating of xi-th and zeta-th 
    ** genotypes.
    */
    /* Indices for O_transform are [xi][zeta][g] */
    KScalar trfm_O[1][1][1] = 
                {{{1.0}}};
    /* A_transform holds proportion of k-th genotype appearing
    ** in progeny of selfed or apomixis mating of xi-th 
    ** genotype.
    */
    /* Under selfing, trfm_S is the diagonal of trfm_O. */
    /* Under apomixis, trfm_A is identity where xi == g */
    /* Indices for trfm_S and trfm_A are [xi][g]; */
    KScalar trfm_A[1][1] = 
                {{1.0}};
    if (K->genotypes != 1) {
        char buf[200];
        sprintf(buf, "%s: K->genotypes=%d, changing to 1\n",
                thisfunction, K->genotypes);
        warning(buf);
        K->genotypes = 1;
    }
    for (xi=0; xi < K->genotypes; xi++) {
        for (g=0; g < K->genotypes; g++) {
            K->trfm_S[xi][g] = trfm_O[xi][xi][g];
            K->trfm_A[xi][g] = trfm_A[xi][g];
        }
    }
    for (xi=0; xi < K->genotypes; xi++) {
        for (zeta=0; zeta < K->genotypes; zeta++) {
            for (g=0; g < K->genotypes; g++) {
                K->trfm_O[xi][zeta][g] = trfm_O[xi][zeta][g];
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
void        set_transform_default       (KConfig K)
/*
** These arrays assume three genotypes controlled by one 
** biallelic locus
*/
{
    const char* thisfunction = "set_transform_default";
    KInt xi, zeta, g;
    /* O_transform holds proportion of g-th genotype appearing 
    ** in progeny of outcrossed mating of xi-th and zeta-th 
    ** genotypes.
    */
    /* Indices for O_transform are [xi][zeta][g] */
    KScalar trfm_O[3][3][3] = 
                {{{1.0,0.0,0.0}, {0.5,0.5,0.0}, {0.0,1.0,0.0}},
                 {{0.5,0.5,0.0}, {0.25,0.5,0.25}, {0.0,0.5,0.5}},
                 {{0.0,1.0,0.0}, {0.0,0.5,0.5}, {0.0,0.0,1.0}}};
    /* A_transform holds proportion of k-th genotype appearing
    ** in progeny of selfed or apomixis mating of xi-th 
    ** genotype.
    */
    /* Under selfing, trfm_S is the diagonal of trfm_O. */
    /* Under apomixis, trfm_A is identity where xi == g */
    /* Indices for trfm_S and trfm_A are [xi][g]; */
    KScalar trfm_A[3][3] = 
                {{1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    if (K->genotypes != 3) {
        char buf[200];
        sprintf(buf, "%s: K->genotypes != 3, changing to 3\n", thisfunction);
        warning(buf);
        K->genotypes = 3;
    }
    for (xi=0; xi < K->genotypes; xi++) {
        for (g=0; g < K->genotypes; g++) {
            K->trfm_S[xi][g] = trfm_O[xi][xi][g];
            K->trfm_A[xi][g] = trfm_A[xi][g];
        }
    }
    for (xi=0; xi < K->genotypes; xi++) {
        for (zeta=0; zeta < K->genotypes; zeta++) {
            for (g=0; g < K->genotypes; g++) {
                K->trfm_O[xi][zeta][g] = trfm_O[xi][zeta][g];
            }
        }
    }
}

/*///////////////////////////////////////////////////////////////*/
int         isOK_repro_transformations  (KConfig K)
/*
** checks to see if the transformation arrays in K are OK
*/
{
    const char* thisfunction = "isOK_repro_transformations";
    {
        char buf[200];
        sprintf(buf, "%s: not yet implemented", thisfunction);
        fatal(buf);
    }
    return (1) ? 1 : 0;
}

