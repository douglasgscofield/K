/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Types etc. to implement nested mutation classes               */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Setup options                                                 */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

#define NESTED

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Datatype extremes                                             */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

#define N_MUTCLASSES        2

/* Note that each mutation class does not get its own extremes,
** these are common to all mutation classes.  This may get further
** refined as the implementation progresses.
*/
#define MAX_MI_N            200
#define MAX_MJ_N            40



/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Types                                                         */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

typedef     int         KMutClass;

typedef     KScalar     KVector_n   [MAX_MI_N+1][MAX_MI_N+1];
typedef     KScalar     KVector_n   [MAX_MI_N+1][MAX_MI_N+1];
typedef     KScalar     KArray_n    [MAX_MI_N+1][MAX_MJ_N+1][MAX_MI_N+1][MAX_MJ_N+1];
typedef     KScalar     KScalar_n   [N_MUTCLASSES];
typedef     KInt        KInt_n      [N_MUTCLASSES];
typedef     KScalar     KArray2_n   [N_MUTCLASSES][MAX_MI_N+1][MAX_MJ_N+1];

typedef     struct struct_KConfig_n*   KConfig_n;

struct      struct_KConfig_n   {
    /* generation */
        int         generation;
    /* parameters of the model */
        KScalar_n   U;
        KScalar     S;
        KScalar     A;
        KScalar     O;
    /* fitness */
        KInt_n      fitness_function;
        KArray2_n   fitness_precomputed;
        KScalar_n   fit_s;
        KScalar_n   fit_h;
        KInt_n      fit_k;
        KScalar_n   fit_d;
        KScalar_n   fit_alpha;
    /* numbers of the classes */
        KInt        MI0;
        KInt        MJ0;
        KInt        MI1;
        KInt        MJ1;
    /* arrays used during model iterations */
        KArray_n    x;
        KArray_n    xp;
        KVector_n   fgam;
        KVector_n   mgam;
        KVector_n   xppo;
        KArray_n    xpps;
        KArray_n    xppa;
        KArray_n    xpp;
        KArray_n    X;
    /* values used to determine if equilibrium is reached */
        KArray_n    x_prevgen;  /* " " genotype frequencies */
        KScalar_n   epsilon;    /* fabs(q[i,j]-q_prevgen[i,j]) must be <= */
    /* arrays used to speed up model operation */
        KVector_n   mut_term;  /* product of mut terms for classes 0 and 1 */
    /* arrays and values used to keep track of statistics */
        KGenotype1  self_fitness;
        KGenotype1  apomixis_fitness;
        KGenotype1  outcross_fitness;
    /* options and extras */
        /* if != 0, then set all elements of K->x that
        ** are < LOADCLASS_TRUNCATE to 0.0 at the
        ** beginning of each generation */
        int         option_truncate;   
        /* if != 0, then print stats in table format */
        int         option_table;
        /* bits are set according to the debug flags
        ** in effect */
        int         debug_flags;
};

#include "K_mutation_n.h"
#include "K_reproduction_n.h"
