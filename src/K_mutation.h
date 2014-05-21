/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/
/* Methods for handling mutation                                 */
/*///////////////////////////////////////////////////////////////*/
/*///////////////////////////////////////////////////////////////*/

/* The maximum number of mutation rates served by mut_term_general */
#define     pois_general_MAX        20

void        compute_mutation        (KConfig K);
void        apply_mutation          (KConfig K, KArray& to, KArray& from);
void        apply_mutation_general  (KConfig K, KScalar U, 
                                     KArray& to, KArray& from);
KScalar     mut_term_general        (KScalar U, KInt x);
KScalar     mut_term                (KConfig K, KInt x);
void        initiate_mut_term       (KConfig K);
KScalar     pois_term               (KScalar U, KInt term);
