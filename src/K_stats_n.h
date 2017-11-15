struct  struct_KStats_n {
    KScalar_n mean_hetloci;
    KScalar_n mean_homloci;
    KScalar_n mean_totmuts;
    KScalar_n var_hetloci;
    KScalar_n var_homloci;
    KScalar_n var_totmuts;
    KScalar_n var_to_mean_totmuts_ratio;
    KScalar   mean_fitness_self_progeny;
    KScalar   mean_fitness_apomixis_progeny;
    KScalar   mean_fitness_outcross_progeny;
    KScalar   population_mean_fitness;
    KScalar   inbreeding_depression;
    KScalar   secondary_selfing_rate;
};


void        stats_all_n             (KConfig_n KN);
void        stats_print_n           (KConfig_n KN);
void        stats_print_table_heading_n (KConfig_n KN);
void        stats_print_table_n     (KConfig_n KN);
void        stats_print_verbose_n   (KConfig_n KN);

void        stats_muts_n    (KConfig_n KN, KMutClass m);

KScalar     stats_mean_fitness_self_progeny_n   (KConfig_n KN);
KScalar     stats_mean_fitness_apomixis_progeny_n   (KConfig_n KN);
KScalar     stats_mean_fitness_outcross_progeny_n   (KConfig_n KN);
KScalar     stats_population_mean_fitness_n (KConfig_n KN);
KScalar     stats_inbreeding_depression_n   (KConfig_n KN);
KScalar     stats_inbreeding_depression_old_n    (KConfig_n KN);
KScalar     stats_secondary_selfing_rate_n  (KConfig_n KN);
