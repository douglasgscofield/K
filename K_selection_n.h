/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Fitness model operations - for nested mutation classes        */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

#define     fitness_function_VALS_n     40

#define     FITNESS_INVALID_n           0
#define     FITNESS_MULTIPLICATIVE_n    1
#define     FITNESS_KONDRASHOV_n        2

typedef     KScalar     (*KPtr_ff_n) (KConfig_n KN, KMutClass m,
                                      KInt i, KInt j);

struct  struct_KFitness_n {
    KInt        thisnum;  /* the registered number 
                          ** of this fitness function */
    KPtr_ff_n   func;     /* pointer to the fitness function */
    char*       name;     /* a name for the function */
    int         mustcompute;  /* is != 0 if it must be computed 
                              ** each time */
};

extern struct struct_KFitness_n fitness_functions_n[fitness_function_VALS_n+1];

void        compute_selection_n (KConfig_n KN);
void        apply_selection_n   (KConfig_n KN, 
                                 KArray_n& to, KArray_n& from);
KScalar     mfitness_computed_n (KConfig_n KN, KMutClass m,
                                 KInt i, KInt j);
KScalar     mfitness_n          (KConfig_n KN, KMutClass m,
                                 KInt i, KInt j);
KScalar     fitness_from_mfitness   (KConfig_n KN, 
                                     KScalar_n mw);
KScalar     fitness_computed_n  (KConfig_n KN, 
                                 KInt i0, KInt j0, 
                                 KInt i1, KInt j1);
KScalar     fitness_n           (KConfig_n KN, 
                                 KInt i0, KInt j0, 
                                 KInt i1, KInt j1);
KScalar     mean_fitness_n      (KConfig_n KN, KArray_n& a);
KScalar     mean_fitness_allprogeny_n   (KConfig_n KN, KArray_n& a);
KScalar     cumulative_fitness_n(KConfig_n KN, KArray_n& a);
KScalar     mean_fitness_outcross_n (KConfig_n KN, KVector_n& out);
KScalar     cumulative_fitness_outcross_n   (KConfig_n KN,
                                             KVector_n& out);
void        initiate_fitness_precomputed_n  (KConfig_n KN);

KInt        register_fitness_function_n (KPtr_ff_n pff_n,
                                         char* ff_name);
void        must_compute_fitness_function_n(KInt ff);
KPtr_ff_n   get_fitness_function_n  (KInt ff);
char*       get_fitness_function_name_n (KInt ff);

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Fitness functions - nested mutation classes                   */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/* The fitness_functions_n[] array is initialized with valid
** entries for the following functions.
*/

KScalar     fitness_kondrashov_n    (KConfig_n KN, KMutClass m,
                                     KInt i, KInt j);
KScalar     fitness_multiplicative_n(KConfig_n KN, KMutClass m,
                                     KInt i, KInt j);

