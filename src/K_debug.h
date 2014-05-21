extern      int                     debug_flags;

void        not_implemented         (const char* function,
                                     const char* msg);
void        check_normalization     (KConfig K, KArray& a,
                                     const char* caller,
                                     const char* array_name);
void        dump_KArray             (KConfig K, KArray& a,
                                     KInt mi, KInt mj, KInt mg);
void        dump_values_KArray      (KConfig KN,
                                     FILE* fp,
                                     KArray& a,
                                     KInt mi, KInt mj, KInt mg);
void        dump_KArray_full        (KConfig K, KArray& a);
void        dump_KVector1           (KConfig K, KVector1& v,
                                     KInt mi, KInt mg);
void        dump_KVector1_full      (KConfig K, KVector1& v);
void        set_debug               (int lvl);

#if !defined(NO_DEBUG)
#define     IF_DEBUG(_lvl_)         if (debug_flags & (0x1 << (_lvl_ - 1)))
#define     DEBUG(_lvl_)            (debug_flags & (0x1 << (_lvl_ - 1)))
#else
#define     IF_DEBUG(_lvl_)         if (0)
#define     DEBUG(_lvl_)            (0)
#endif

#define     DEBUG_TRACE1                1
#define     DEBUG_TRACE2                2
#define     DEBUG_TRACE3                3
#define     DEBUG_EQUILIBRIUM           4
#define     DEBUG_EQUILIBRIUM_DETAIL    5
#define     DEBUG_TRUNCATE              6
#define     DEBUG_NORMALIZATION         7
#define     DEBUG_OUTCROSS              8
#define     DEBUG_FOLLOW_EQUILIBRIUM    9
#define     DEBUG_GENERATIONS           10
#define     DEBUG_LETHALS               11
#define     DEBUG_TRUNCATE_DETAIL       12

