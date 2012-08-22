/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Types etc. to implement nested mutation classes               */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Setup options                                                 */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

#define NESTED

#define GENERATION_CUTOFF_n 500
/*
**  Define this to make a bounds-checking function call, rather than
**  direct KN->mut_term[][] array access, in apply_mutation_n()
*/
/* #define MUTTERM_FUNCTION_n 1 */
/*
**  Define this to include within-loop debug statements in 
**  time-critical areas.
*/
/* #define TIMECRITICAL_INCLUDEDEBUG_n 1 */


/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Datatype extremes                                             */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

#define N_MUTCLASSES        2

/* Note that each mutation class does not get its own extremes,
** these are common to all mutation classes.  This may get furthe
** refined as the implementation progresses.
*/
#define MAX_MI_n            50
#define MAX_MJ_n            10



/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Types                                                         */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

typedef int KMutClass;

typedef KScalar KVector_n[MAX_MI_n + 1][MAX_MI_n + 1];
typedef KScalar KArray_n[MAX_MI_n + 1][MAX_MJ_n + 1][MAX_MI_n + 1][MAX_MJ_n +
																   1];
typedef KScalar KScalar_n[N_MUTCLASSES];
typedef KInt KInt_n[N_MUTCLASSES];
typedef KScalar KArray2_n[N_MUTCLASSES][MAX_MI_n + 1][MAX_MJ_n + 1];

typedef struct struct_KConfig_n *KConfig_n;

struct struct_KConfig_n
{
	KMutClass mutclasses;
	KInt_n is_lethal;			/* 1 if mutation class is lethal */
	KInt_n createlethal;		/* 1 if lethal and should create progeny with j>0 */
	/* generation */
	KInt generation;
	/* parameters of the model */
	KScalar_n U;
	KScalar S;
	KScalar A;
	KScalar O;
	/* fitness */
	KInt_n fitness_function;
	KArray2_n fitness_precomputed;
	KScalar_n fit_s;
	KScalar_n fit_h;
	KInt_n fit_k;
	KScalar_n fit_d;
	KScalar_n fit_alpha;
	/* numbers of the classes */
	KInt MI0;
	KInt MJ0;
	KInt MI1;
	KInt MJ1;
	/* arrays used during model iterations */
	int current_x;				/* set to the 'number' of the array just set */
#define KN_CURRENT_X1               1
#define KN_CURRENT_X2               2
#define KN_CURRENT_GAMETES_X2_TO_X1 3
#define KN_CURRENT_PROGENY_TO_X1    4
	KArray_n x1;
	KArray_n x2;
	KVector_n fgam;
	KVector_n mgam;
	/* values used to determine if equilibrium is reached */
	KArray_n x_prevgen;			/* " " genotype frequencies */
	KScalar epsilon;
	/* arrays used to speed up model operation */
	KVector_n mut_term;			/* product of mut terms for classes 0 and 1 */
	/* arrays and values used to keep track of statistics */
	KGenotype1 self_fitness;
	KGenotype1 apomixis_fitness;
	KGenotype1 outcross_fitness;
	/* options and extras */
	/* if != 0, then set all elements of K->x that
	 ** are < LOADCLASS_TRUNCATE to 0.0 at the
	 ** beginning of each generation */
	int option_truncate;
	/* if != 0, then print stats in table format */
	int option_table;
	/* if != 0, then do not shortcut computations if lethals involved */
	int option_nolethal;
	/* bits are set according to the debug flags
	 ** in effect */
	int debug_flags;
	/* if != 0, create a savefile */
	int save_savefile;
	char *save_savefile_name;
	/* if != 0, load a savefile */
	int load_savefile;
	char *load_savefile_name;
};

int main_nested(
	int argc,
	char *argv[]);

#include "K_cmdline_n.h"
#include "K_debug_n.h"
#include "K_equilibrium_n.h"
#include "K_initiate_n.h"
#include "K_mutation_n.h"
#include "K_nextgen_n.h"
#include "K_prevgen_n.h"
#include "K_reproduction_n.h"
#include "K_savefile_n.h"
#include "K_selection_n.h"
#include "K_stats_n.h"
#include "K_util_n.h"
