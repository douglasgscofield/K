/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Fitness model operations                                      */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

#define     fitness_function_VALS      40

#define        FITNESS_INVALID         0
#define        FITNESS_MULTIPLICATIVE  1
#define        FITNESS_KONDRASHOV        2

typedef     KScalar     (*KPtr_ff)  (KConfig K, KInt i, KInt j, 
                                     KInt g);

struct  struct_KFitness {
    KInt        thisnum;  /* the registered number 
                          ** of this fitness function */
    KPtr_ff     func;     /* pointer to the fitness function */
    char*       name;     /* a name for the function */
    int         mustcompute;  /* is != 0 if it must be computed 
                              ** each time */
};

extern struct struct_KFitness fitness_functions[fitness_function_VALS+1];

/*///////////////////////////////////////////////////////////////*/

void        compute_selection       (KConfig K);
void        apply_selection         (KConfig K, KArray& to,
                                     KArray& from);

KScalar     fitness_computed        (KConfig K, KInt i, KInt j, 
                                     KInt g);
KScalar     fitness                 (KConfig K, KInt i, KInt j, 
                                     KInt g);
KScalar     mean_fitness            (KConfig K, KArray& a);
KScalar     mean_fitness_allprogeny (KConfig K, KArray& a);
KScalar     cumulative_fitness      (KConfig K, KArray& a);
void        initiate_fitness_precomputed    (KConfig K);

KInt        register_fitness_function   (KPtr_ff pff,
                                         char* ff_name);
void        must_compute_fitness_function   (KInt ff);
KPtr_ff     get_fitness_function    (KInt fm);
char*       get_fitness_function_name   (KInt ff);

/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Fitness functions                                             */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/* The fitness_functions[] array is initialized with valid
** entries for the following functions.
*/

KScalar     fitness_multiplicative  (KConfig K, KInt i, KInt j, 
                                     KInt g);
KScalar     fitness_kondrashov      (KConfig K, KInt i, KInt j, 
                                     KInt g);
                                     
