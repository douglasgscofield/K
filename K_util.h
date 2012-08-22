KScalar sum_KArray(
	KConfig K,
	KArray & A);
KScalar sum_KVector1(
	KConfig K,
	KVector1 & v);
KScalar sum_KArray2(
	KConfig K,
	KArray2 & A);
void copy_KArray(
	KConfig K,
	KArray & to,
	KArray & from);
int isOK_KArray(
	KConfig K,
	KArray & a);
void normalize_KArray(
	KConfig K,
	KArray & a);
void truncate_KArray(
	KConfig K,
	KArray & a,
	KScalar v);

KConfig new_KConfig(
	KInt MI,
	KInt MJ);
void free_KConfig(
	KConfig K);
void new_genotypes(
	KConfig K,
	KInt g);
void free_genotypes(
	KConfig K);
void fatal(
	const char *msg);
void warning(
	const char *msg);

void fill_KArray(
	KConfig K,
	KArray & a,
	KScalar val);
