void check_normalization_n(
	KConfig_n KN,
	KArray_n & a,
	const char *caller,
	char *array_name);
void dump_KArray_n(
	KConfig_n KN,
	KArray_n & a,
	KInt mi0,
	KInt mj0,
	KInt mi1,
	KInt mj1);
void dump_values_KArray_n(
	KConfig_n KN,
	FILE * fp,
	KArray_n & a,
	KInt mi0,
	KInt mj0,
	KInt mi1,
	KInt mj1);
void dump_KArray_full_n(
	KConfig_n KN,
	KArray_n & a);
void dump_KVector_n(
	KConfig_n KN,
	KVector_n & v,
	KInt mi0,
	KInt mi1);
void dump_KVector_full_n(
	KConfig_n KN,
	KVector_n & v);
