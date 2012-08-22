KScalar sum_KArray_n(
	KConfig_n KN,
	KArray_n & a);
KScalar sum_KVector_n(
	KConfig_n KN,
	KVector_n & v);
void normalize_KArray_n(
	KConfig_n KN,
	KArray_n & a);
void truncate_KArray_n(
	KConfig_n KN,
	KArray_n & a,
	KScalar v);
void fill_KArray_n(
	KConfig_n KN,
	KArray_n & a,
	KScalar val);
void fill_KVector_n(
	KConfig_n KN,
	KVector_n & v,
	KScalar val);
void copy_KArray_n(
	KConfig_n KN,
	KArray_n & to,
	KArray_n & from);
int isOK_KArray_n(
	KConfig_n KN,
	KArray_n & a);
int isOK_KVector_n(
	KConfig_n KN,
	KVector_n & v);
void *alloc_KArray_n(
	void);
void *alloc_KVector_n(
	void);
void free_KArray_n(
	void *p);
void free_KVector_n(
	void *p);
