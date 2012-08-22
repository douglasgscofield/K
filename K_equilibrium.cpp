#include "K.h"


/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */
/* Functions for determining equilibrium                         */
/* ///////////////////////////////////////////////////////////// */
/* ///////////////////////////////////////////////////////////// */


/* ///////////////////////////////////////////////////////////// */
int
is_equilibrium(
	KConfig K)
/* 
** return 1 if K satisfies equilibrium constraints as defined by use
*/
{
	const char *thisfunction = "is_equilibrium";
	if (diff_epsilon_KArray(K, K->x, K->x_prevgen, K->epsilon)) {
		IF_DEBUG(DEBUG_TRACE1) printf("is_equilibrium: no equilibrium\n");
		return 0;
	} else {
		IF_DEBUG(DEBUG_TRACE1) printf("is_equilibrium: at equilibrium\n");
		return 1;
	}
}

/* ///////////////////////////////////////////////////////////// */
int
diff_epsilon_KArray(
	KConfig K,
	KArray & a1,
	KArray & a2,
	KScalar epsilon)
/*
** Checks to see if any of the proportions in a1 and a2 differ by
** more than epsilon.  If true, returns 1, else 0.  
*/
{
	const char *thisfunction = "diff_epsilon_KArray";
	KInt i, j, g;
	IF_DEBUG(DEBUG_FOLLOW_EQUILIBRIUM) {
		KInt n, v, xi, l, lam, zeta;
		KScalar diff, maxdiff = 0.0;
		for (n = 0; n < K->MI; n++) {
			for (v = 0; v < K->MJ; v++) {
				for (xi = 0; xi < K->genotypes; xi++) {
					diff = fabs(a1[n][v][xi] - a2[n][v][xi]);
					if (diff > maxdiff) {
						maxdiff = diff;
						l = n;
						lam = v;
						zeta = xi;
					}
				}
			}
		}
		printf
			("diff_epsilon_KArray: epsilon = %lg, max diff @ [%d][%d][%d] = %lg\n",
			 epsilon, l, lam, zeta, maxdiff);
	}
	IF_DEBUG(DEBUG_EQUILIBRIUM) {
		KInt n, v, xi;
		KInt l = 0, lam = 0, zeta = 0;
		KScalar diff, maxdiff = 0.0;
		IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
			printf
				("diff_epsilon_KArray: BEGIN where (a1[][][]-a2[][][] != 0.0)\n");
			printf("i\tj\tg\tdiff\n");
		}
		for (n = 0; n < K->MI; n++) {
			for (v = 0; v < K->MJ; v++) {
				for (xi = 0; xi < K->genotypes; xi++) {
					diff = fabs(a1[n][v][xi] - a2[n][v][xi]);
					if (diff > maxdiff) {
						maxdiff = diff;
						l = n;
						lam = v;
						zeta = xi;
					}
					IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
						if (n < 10 && v < 10) {
							diff = a1[n][v][xi] - a2[n][v][xi];
							if (diff != 0.0) {
								printf("%d\t%d\t%d\t%lg\n", n, v, xi, diff);
							}
						}
					}
				}
			}
		}
		IF_DEBUG(DEBUG_EQUILIBRIUM_DETAIL) {
			printf("diff_epsilon_KArray: END (a1[][][]-a2[][][] != 0.0)\n");
		}
		printf
			("diff_epsilon_KArray: epsilon = %lg, max diff @ [%d][%d][%d] = %lg\n",
			 epsilon, l, lam, zeta, maxdiff);
	}
	/* TODO: sort out whether these should be "<=" or "<" */
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			for (g = 0; g < K->genotypes; g++) {
				if (fabs(a1[i][j][g] - a2[i][j][g]) > epsilon) {
					/* only need one element that exceeds epsilon */
					return 1;
				}
			}
		}
	}
	return 0;
}
