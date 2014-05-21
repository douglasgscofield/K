struct  struct_KStats {
    KScalar   mean_hetloci;
    KScalar   mean_homloci;
    KScalar   mean_totmuts;
    KScalar   var_hetloci;
    KScalar   var_homloci;
    KScalar   var_totmuts;
    KScalar   var_to_mean_totmuts_ratio;
    KScalar   mean_fitness_self_progeny;
    KScalar   mean_fitness_apomixis_progeny;
    KScalar   mean_fitness_outcross_progeny;
    KScalar   population_mean_fitness;
    KScalar   inbreeding_depression;
    KScalar   secondary_selfing_rate;
};


void        stats_all               (KConfig K);
void        stats_print             (KConfig K);
void        stats_print_table_heading   (KConfig K);
void        stats_print_table       (KConfig K);
void        stats_print_verbose     (KConfig K);
void        stats_muts      (KConfig K);
KScalar     stats_mean_fitness_self_progeny (KConfig K);
KScalar     stats_mean_fitness_apomixis_progeny (KConfig K);
KScalar     stats_mean_fitness_outcross_progeny (KConfig K);
KScalar     stats_population_mean_fitness   (KConfig K);
KScalar     stats_inbreeding_depression (KConfig K);
KScalar     stats_inbreeding_depression_old	(KConfig K);
KScalar     stats_secondary_selfing_rate(KConfig K);
