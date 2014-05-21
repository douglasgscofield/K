void        compute_self_progeny_n  (KConfig_n KN);
void        apply_self_progeny_n    (KConfig_n KN,
                                     KArray_n& to, KArray_n& from);
void        compute_apomixis_progeny_n  (KConfig_n KN);
void        apply_apomixis_progeny_n    (KConfig_n KN,
                                         KArray_n& to, KArray_n& from);
void        compute_gametes_n       (KConfig_n KN);
void        apply_gametes_n     (KConfig_n KN, 
                                 KVector_n& mgam, KVector_n& fgam,
                                 KArray_n& from);
void        apply_gametes_full_n    (KConfig_n KN, 
                                     KVector_n& mgam, KVector_n& fgam,
                                     KArray_n& from);
void        adjust_gametes_n    (KConfig_n KN, 
                                 KVector_n& mgam, KVector_n& fgam,
                                 KArray_n& from);
void        compute_zygotes_n   (KConfig_n KN);
void        apply_zygotes_n     (KConfig_n KN, KArray_n& to,
                                 KVector_n& mgam, KVector_n& fgam);
void        compute_summed_progeny_n    (KConfig_n KN);
void        apply_summed_progeny_n  (KConfig_n KN,
                                     KArray_n& to, KArray_n& from);
void        set_repro_n     (KConfig_n KN, KScalar s, KScalar a);
KScalar     s_self_n        (KInt i0, KInt j0, KInt i1, KInt j1, 
                             KInt n0, KInt v0, KInt n1, KInt v1);
KScalar     a_apomixis_n    (KInt i0, KInt j0, KInt i1, KInt j1, 
                             KInt n0, KInt v0, KInt n1, KInt v1);
