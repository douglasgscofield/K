#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <valarray>
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Global variable external declarations -- all defined in K.c  
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Setup options                                                 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

#define DEBUG

//#define POISSON_TRUNCATE        1.0e-100
//#define LOADCLASS_TRUNCATE      1.0e-100
#define POISSON_TRUNCATE        1.0e-50
#define LOADCLASS_TRUNCATE      1.0e-50
#define NORMALIZATION_TOLERANCE 1.0e-9
#define GENERATION_CUTOFF       2000
/*
**  Define this to make a bounds-checking function call, rather than
**  direct K->mut_term[] array access, in apply_mutation()
*/
/* #define MUTTERM_FUNCTION 1 */
/*
**  Define this to include within-loop debug statements in 
**  time-critical areas.
*/
/* #define TIMECRITICAL_INCLUDEDEBUG 1 */



/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Datatype extremes                                             
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

#define MAX_MI              500
#define MAX_MJ              40
#define MAX_G               1
/*
#define MAX_MI              300
#define MAX_MJ              200
#define MAX_G               1
*/


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Types                                                         
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

//typedef     std::valarray< std::valarray< KScalar > > Matrix2;
//typedef     std::valarray< std::valarray< std::valarray< KScalar > > > Matrix3;
//typedef     std::valarray< std::valarray< KScalar > >                  KVector;
//typedef     std::valarray< KScalar >                                   KVector1;
//typedef     std::valarray< std::valarray< std::valarray< KScalar > > > KArray;
//typedef     std::valarray< std::valarray< std::valarray< KScalar > > > KArray2;
//typedef     std::valarray< KScalar >                                   KGenotype1;
//typedef     std::valarray< std::valarray< KScalar > >                  KGenotype2;
//typedef     std::valarray< std::valarray< std::valarray< KScalar > > > KGenotype3;

typedef     int             KInt;
typedef     double          KScalar;
typedef     KScalar         KVector     [MAX_MI+1][MAX_G];
typedef     KScalar         KVector1    [MAX_MI+1];
typedef     KScalar         KArray      [MAX_MI+1][MAX_MJ+1][MAX_G];
typedef     KScalar         KArray2     [MAX_MI+1][MAX_MJ+1];
typedef     KScalar         KGenotype1  [MAX_G];
typedef     KScalar         KGenotype2  [MAX_G][MAX_G];
typedef     KScalar         KGenotype3  [MAX_G][MAX_G][MAX_G];


/////////////////////////////////////////////////////////////////
// Note that the KConfig type is a pointer to struct_KConfig     
/////////////////////////////////////////////////////////////////

typedef     struct struct_KConfig*      KConfig;

struct      struct_KConfig     {
        KInt        is_lethal;    // 1 if mutation class is lethal
        KInt        createlethal; // 1 if lethal and should create progeny with j>0
        // generation
        KInt        generation;
        // type of mating system to be used
        KInt        mating_outcross;    // nonzero if outcrossing
        KInt        mating_self;        // nonzero if selfing
        KInt        mating_apomixis;    // nonzero if apomixis
        // parameters of the model
        KScalar     U;   // genomic mutation rate per generation
        KGenotype1  S;   // selfing rate
        KGenotype1  D_S; // discount from selfing
        KGenotype1  A;   // apomixis rate
        KGenotype1  D_A; // discount from apomixis
        KGenotype1  O;   // outcrossing rate = 1 - (S + A)
        // fitness
        KInt        fitness_function;
        KArray      fitness_precomputed;  // computed ahead of time
        KScalar     fit_s;   // multiplicative selection coefficient
        KScalar     fit_h;   // multiplicative dominance
        KInt        fit_k;   // Kondrashov's truncation point
        KScalar     fit_d;   // Kondrashov's dominance
        KScalar     fit_alpha;  // Kondrashov's alpha for shape of fitness
        // numbers of the classes
        KInt        genotypes;  // number of genotypes minus 1
        KInt        MI;         // max heterozygous mutations; 0 also a class
        KInt        MJ;         // max homozygous mutations; 0 also a class
        // arrays used to determine genotype transformations
        KGenotype2  trfm_S;
        KGenotype2  trfm_A;
        KGenotype3  trfm_O;
        // arrays used during model iterations
        KGenotype1  initial_genotype_freqs;  // freqs at start of model run
        KArray      x;      // adult frequencies
        KArray      xp;     // post-mutation frequencies
        KVector1    fgam;   // should be KVector2 or something...
        KVector1    mgam;   // should be KVector2 or something...
        KArray      xppo;   // outcrossed progeny frequencies
        KArray      xpps;   // selfed progeny frequencies
        KArray      xppa;   // apomictic progeny frequencies
        KArray      xpp;    // summed progeny frequencies
        KArray      X;      // post-selection frequencies
        // values used to determine if equilibrium is reached
        KArray      x_prevgen;  // " " genotype frequencies
        KScalar     epsilon;    // fabs(q[i,j]-q_prevgen[i,j]) must be <=
        // scalars and arrays used by mating-system code 
        KGenotype1  rsrc_SO;    // portion of resources to selfed ovules
        KGenotype1  rsrc_AO;    // portion of resources to apomictic ovules
        KGenotype1  rsrc_OO;    // portion of resources to outcross ovules
        KGenotype1  rsrc_OP;    // portion of resources to outcross pollen
        KScalar     F_female;
        KScalar     F_male;
        KArray2     gamma;
        KArray2     alpha;
        KArray2     beta;
        KArray2     rho;
        KArray      SO;
        KArray      AO;
        KArray      OO;
        KArray      OP;
        // arrays used to speed up model operation
        KVector1    mut_term;
        // arrays and values used to keep track of statistics
        KGenotype1  self_fitness;
        KGenotype1  apomixis_fitness;
        KGenotype1  outcross_fitness;

        // options and extras */

        // if != 0, then set all elements of K->x that are < LOADCLASS_TRUNCATE
        // to 0.0 at the beginning of each generation
        int         option_truncate;   
        // if != 0, then print stats in table format
        int         option_table;
        // if != 0, then do not shortcut computations if lethals involved
        int         option_nolethal;
        // bits are set according to the debug flags in effect
        int         debug_flags;
        // if != 0, create a savefile
        int         save_savefile;
        char*       save_savefile_name;
        // if != 0, load a savefile
        int         load_savefile;
        char*       load_savefile_name;
};

int main_unnested(int argc, char* argv[]);

#include "K_cmdline.h"
#include "K_debug.h"
#include "K_equilibrium.h"
#include "K_genotypes.h"
#include "K_initiate.h"
#include "K_math.h"
#include "K_mutation.h"
#include "K_nextgen.h"
#include "K_reproduction.h"
#include "K_savefile.h"
#include "K_selection.h"
#include "K_stats.h"
#include "K_util.h"
// #include "K_n.h"

