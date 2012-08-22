#include "K.h"

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Utility routines for debugging                                
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// Global variables

int debug_flags = 0x0;

/////////////////////////////////////////////////////////////////
void
not_implemented(
	const char *function,
	const char *msg)
{
	char buf[500];
	sprintf(buf, "%s: %s: not implemented", function, msg);
	fatal(buf);
}

/////////////////////////////////////////////////////////////////
void
check_normalization(
	KConfig K,
	KArray & a,
	const char *caller,
	const char *array_name)
{
	const char *thisfunction = "check_normalization";
	IF_DEBUG(DEBUG_NORMALIZATION) {
		if (!isOK_KArray(K, a)) {
			printf("%s: %s not normalized: %f, tolerance=%lg\n", caller,
				   array_name, sum_KArray(K, a), NORMALIZATION_TOLERANCE);
		}
	}
}

/////////////////////////////////////////////////////////////////
void
dump_KArray(
	KConfig K,
	KArray & a,
	KInt mi,
	KInt mj,
	KInt mg)
{
	const char *thisfunction = "dump_KArray";
	KInt i, j, g;
	if (mi > K->MI || mi < 0)
		mi = K->MI;
	if (mj > K->MJ || mj < 0)
		mj = K->MJ;
	if (mg > K->genotypes || mg < 1)
		mg = K->genotypes;
	printf("%s BEGIN KArray @ 0x%08x; max i=%d max j=%d max g=%d\n",
		   thisfunction, a, mi, mj, mg);
	printf("i\tj\tg\ta[i][j][g]\n");
	for (i = 0; i <= mi; i++) {
		for (j = 0; j <= mj; j++) {
			for (g = 0; g < mg; g++) {
				printf("%d\t%d\t%d\t%lg\n", i, j, g, a[i][j][g]);
			}
		}
	}
	printf("%s END for KArray a @ 0x%8x\n", thisfunction, a);
}

/////////////////////////////////////////////////////////////////
void
dump_values_KArray(
	KConfig K,
	FILE * fp,
	KArray & a,
	KInt mi,
	KInt mj,
	KInt mg)
{
	KInt i, j, g;
	for (i = 0; i <= mi; i++) {
		for (j = 0; j <= mj; j++) {
			for (g = 0; g < mg; g++) {
				fprintf(fp, "%d\t%d\t%d\t%lg\n", i, j, g, a[i][j][g]);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////
void
dump_KArray_full(
	KConfig K,
	KArray & a)
{
	const char *thisfunction = "dump_KArray_full";
	printf("%s BEGIN KArray @ 0x%08x; dump_KArray follows\n", thisfunction, a);
	dump_KArray(K, a, K->MI, K->MJ, K->genotypes);
	printf("%s END for KArray a @ 0x%8x\n", thisfunction, a);
}

/////////////////////////////////////////////////////////////////
void
dump_KVector1(
	KConfig K,
	KVector1 & v,
	KInt mi,
	KInt mg)
{
	const char *thisfunction = "dump_KVector1";
	KInt i, g;
	if (K->genotypes != 1) {
		not_implemented(thisfunction, "genotypes != 1");
	}
	if (mi > K->MI || mi < 0)
		mi = K->MI;
	if (mg > K->genotypes || mg < 1)
		mg = K->genotypes;
	printf("%s BEGIN KVector1 @ 0x%08x; max i=%d max g=%d\n", thisfunction, v,
		   mi, mg);
	printf("i\tg\tv[i] ...[g=0]\n");
	for (i = 0; i <= mi; i++) {
		for (g = 0; g < mg; g++) {
			printf("%d\t%d\t%lg\n", i, g, v[i]);
		}
	}
	printf("%s END for KVector1 v @ 0x%8x\n", thisfunction, v);
}

/////////////////////////////////////////////////////////////////
void
dump_KVector1_full(
	KConfig K,
	KVector1 & v)
{
	const char *thisfunction = "dump_KVector1_full";
	printf("%s BEGIN KArray @ 0x%08x; dump_KVector1 follows\n", thisfunction,
		   v);
	dump_KVector1(K, v, K->MI, K->genotypes);
	printf("%s END for KVector1 v @ 0x%8x\n", thisfunction, v);
}

/////////////////////////////////////////////////////////////////
void
set_debug(
	int lvl)
{
	const char *thisfunction = "set_debug";
	if (lvl == 0) {
		debug_flags = 0;
	} else {
		debug_flags |= (0x1 << (lvl - 1));
	}
}
