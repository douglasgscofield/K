KConfig     initiate_KConfig        (void);
void        initiate_load_classes   (KConfig K, KInt MI, KInt MJ);
void        initiate_genotypes      (KConfig K, KInt g);
KConfig     initiate_quick          (KInt MI, KInt MJ, KInt g,
                                     KScalar U,
                                     KScalar s, KScalar h,
                                     KScalar S);
void        initiate_model_state    (KConfig K);
void        compute_adults_initial  (KConfig K);
