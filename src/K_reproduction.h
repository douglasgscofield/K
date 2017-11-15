void        compute_self_progeny    (KConfig K);
void        apply_self_progeny  (KConfig K, KArray& to, KArray& from);
void        apply_self_progeny_stats    (KConfig K, 
                                         KArray& to, KScalar fromval,
                                         KInt fi, KInt fj, KInt fg);
void        compute_apomixis_progeny    (KConfig K);
void        apply_apomixis_progeny  (KConfig K, KArray& to,
                                                KArray& from);
void        compute_gametes         (KConfig K);
void        apply_gametes       (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from);
void        apply_gametes_full  (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from);
void        adjust_gametes      (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KArray& from);
void        apply_gametes_stats (KConfig K, 
                                 KVector1& mgam, KVector1& fgam,
                                 KScalar fromval,
                                 KInt fi, KInt fj, KInt fg);
void        compute_zygotes         (KConfig K);
void        apply_zygotes       (KConfig K, KArray& to,
                                 KVector1& mgam, KVector1& fgam);
void        compute_summed_progeny  (KConfig K);
void        apply_summed_progeny    (KConfig K, KArray& to,
                                     KArray& froms,
                                     KArray& froma,
                                     KArray& fromo);
void        compute_auxiliary_values    (KConfig K);
void        set_repro_resources     (KConfig K);
void        set_repro               (KConfig K, 
                                     KScalar s, KScalar ds,
                                     KScalar a, KScalar da);
KScalar     s_self                  (KInt i, KInt j, 
                                     KInt n, KInt v);
KScalar     a_apomixis              (KInt i, KInt j,
                                     KInt n, KInt v);
KScalar     o_outcross              (KInt i, KInt j, 
                                     KInt n, KInt v, 
                                     KInt l, KInt lam);

void        set_repro_resources_kondrashov  (KConfig K);
void        compute_auxiliary_values_kondrashov (KConfig K);
void        compute_self_progeny_kondrashov (KConfig K);
void        compute_apomixis_progeny_kondrashov (KConfig K);
void        compute_outcross_progeny_kondrashov (KConfig K);
