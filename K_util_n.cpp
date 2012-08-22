#include "K.h"

/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Utility routines for debugging - for nested mutation classes  */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////// */
KScalar
sum_KArray_n(
	KConfig_n KN,
	KArray_n & a)
/*
** Computes sum of all array elements in a according to KN
*/
{
	const char *thisfunction = "sum_KArray_n";
	KInt i0, j0, i1, j1;
	KScalar ans = 0.0;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					ans += a[i0][j0][i1][j1];
				}
			}
		}
	}
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
KScalar
sum_KVector_n(
	KConfig_n KN,
	KVector_n & v)
/*
** Computes sum of all vector elements in v according to KN
*/
{
	const char *thisfunction = "sum_KVector_n";
	KInt i0, i1;
	KScalar ans = 0.0;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (i1 = 0; i1 <= KN->MI1; i1++) {
			ans += v[i0][i1];
		}
	}
	return ans;
}

/* ///////////////////////////////////////////////////////////// */
void
normalize_KArray_n(
	KConfig_n KN,
	KArray_n & a)
/*
** Normalizes a according to KN so that it sums to 1
*/
{
	const char *thisfunction = "normalize_KArray_n";
	KInt i0, j0, i1, j1;
	KScalar sum = sum_KArray_n(KN, a);
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	IF_DEBUG(DEBUG_NORMALIZATION) {
		printf("%s: prenormalization, sum_KArray_n = %lg\n", thisfunction,
			   sum);
	}
	if (sum != 1.0) {
		for (i0 = 0; i0 <= KN->MI0; i0++) {
			for (j0 = 0; j0 <= KN->MJ0; j0++) {
				for (i1 = 0; i1 <= KN->MI1; i1++) {
					for (j1 = 0; j1 <= KN->MJ1; j1++) {
						a[i0][j0][i1][j1] /= sum;
					}
				}
			}
		}
	}
	IF_DEBUG(DEBUG_NORMALIZATION) {
		printf("%s: postnormalization, sum_KArray_n = %lg\n", thisfunction,
			   sum_KArray_n(KN, a));
	}
}

/* ///////////////////////////////////////////////////////////// */
void
truncate_KArray_n(
	KConfig_n KN,
	KArray_n & a,
	KScalar v)
/*
** Truncate a so that any load class-genotype that 
** is < v becomes 0.0
*/
{
	const char *thisfunction = "truncate_KArray_n";
	KInt i0, j0, i1, j1;
	KInt num = 0;
	KScalar sum = 0.0;
	IF_DEBUG(DEBUG_TRACE1) printf("%s\n", thisfunction);
	IF_DEBUG(DEBUG_TRUNCATE) {
		printf("%s: truncating values lower than v=%lg \n", thisfunction, v);
	}
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					if (a[i0][j0][i1][j1] > 0.0 && a[i0][j0][i1][j1] < v) {
						sum += v - a[i0][j0][i1][j1];
						num++;
						a[i0][j0][i1][j1] = 0.0;
						IF_DEBUG(DEBUG_TRUNCATE) {
							printf("%s: truncated a[%d][%d][%d][%d]\n",
								   thisfunction, i0, j0, i1, j1);
						}
					}
				}
			}
		}
	}
	IF_DEBUG(DEBUG_TRUNCATE) {
		printf("%s: truncated %d elements, whose sum is %lg\n", thisfunction,
			   num, sum);
	}
}

/* ///////////////////////////////////////////////////////////// */
void
fill_KArray_n(
	KConfig_n KN,
	KArray_n & a,
	KScalar val)
{
	const char *thisfunction = "fill_KArray_n";
	KInt i0, j0, i1, j1;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					a[i0][j0][i1][j1] = val;
				}
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////// */
void
fill_KVector_n(
	KConfig_n KN,
	KVector_n & v,
	KScalar val)
{
	const char *thisfunction = "fill_KVector_n";
	KInt i0, i1;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (i1 = 0; i1 <= KN->MI1; i1++) {
			v[i0][i1] = val;
		}
	}
}

/* ///////////////////////////////////////////////////////////// */
void
copy_KArray_n(
	KConfig_n KN,
	KArray_n & to,
	KArray_n & from)
{
	const char *thisfunction = "copy_KArray_n";
	KInt i0, j0, i1, j1;
	for (i0 = 0; i0 <= KN->MI0; i0++) {
		for (j0 = 0; j0 <= KN->MJ0; j0++) {
			for (i1 = 0; i1 <= KN->MI1; i1++) {
				for (j1 = 0; j1 <= KN->MJ1; j1++) {
					to[i0][j0][i1][j1] = from[i0][j0][i1][j1];
				}
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////// */
int
isOK_KArray_n(
	KConfig_n KN,
	KArray_n & a)
/*
** Checks to see if the classes in a sum to 1
*/
{
	const char *thisfunction = "isOK_KArray_n";
	KScalar t1;
	t1 = fabs(sum_KArray_n(KN, a) - 1.0);
	return (t1 < NORMALIZATION_TOLERANCE) ? 1 : 0;
}

/* ///////////////////////////////////////////////////////////// */
int
isOK_KVector_n(
	KConfig_n KN,
	KVector_n & v)
/*
** checks to see if the classes in v sum to 1
*/
{
	const char *thisfunction = "isOK_KVector_n";
	KScalar t1;
	t1 = fabs(sum_KVector_n(KN, v) - 1.0);
	return (t1 < NORMALIZATION_TOLERANCE) ? 1 : 0;
}

/* ///////////////////////////////////////////////////////////// */
void *
alloc_KArray_n(
	void)
{
	return (void *) calloc((size_t) 1, sizeof(KArray_n));
}

/* ///////////////////////////////////////////////////////////// */
void *
alloc_KVector_n(
	void)
{
	return (void *) calloc((size_t) 1, sizeof(KVector_n));
}

/* ///////////////////////////////////////////////////////////// */
void
free_KArray_n(
	void *p)
{
	free(p);
}

/* ///////////////////////////////////////////////////////////// */
void
free_KVector_n(
	void *p)
{
	free(p);
}
