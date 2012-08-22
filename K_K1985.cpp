#include "K.h"


// Trying to reconstruct the distant past... I believe these are routines
// that mimic Kondrashov's original 1985 model.

/* ///////////////////////////////////////////////////////////////////// */
KConfig
initiate_K1985(
	KInt MI,
	KInt MJ,
	KScalar U,
	KScalar d,
	KInt k,
	KScalar alpha)
{
	KConfig K;
	K = initiate_KConfig(MI, MJ);
	K->type = TYPE_KONDRASHOV;
	/* K->mating = MATING_OUTCROSS_SELF; */
	K->U = U;
	K->d = d;
	K->k = k;
	K->alpha = alpha;
	initiate_setup(K);
	return K;
}

void
set_fitness_K1985(
	KConfig K,
	KScalar U,
	KScalar d,
	KInt k,
	KScalar alpha)
{
}

void
set_repro_K1985(
	KConfig K,
	KInt genot,
	KScalar y[],
	KScalar z[])
{
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_q_initial_K1985(
	KConfig K)
/*
** Initiate adult frequencies
*/
{
	KInt i, j, g;
	KScalar sum;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			if (i == 0 && j == 0) {
				K->q[i][j] = 1;
			} else {
				K->q[i][j] = 0;
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_q_K1985_from_x(
	KConfig K,
	KArray2 qk,
	KArray x)
/*
** Compute KArray2 members ("q" arrays) as described in Kondrashov 1985 
** by summing KArray "x" elements, as used here, across g for each i,j.
*/
{
	KInt i, j, g;
	KScalar sum;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			sum = 0.0;
			for (g = 0; g < K->genotypes; g++) {
				sum += x[g][i][j];
			}
			qk[i][j] = sum;
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_x_K1985_from_x(
	KConfig K,
	KArray xk,
	KArray x)
/*
** Compute KArray members ("x" arrays) as described in Kondrashov 1985 
** by normalizing KArray "x" elements, as used here, across g for each i,j.
*/
{
	KInt i, j, g;
	KScalar sum;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			sum = 0.0;
			for (g = 0; g < K->genotypes; g++) {
				sum += x[g][i][j];
			}
			for (g = 0; g < K->genotypes; g++) {
				xk[g][i][j] = x[g][i][j] / sum;
			}

		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_q_K1985(
	KConfig K)
/*
** Compute K->q array equivalent to Kondrashov's q.
** K->q is correct upon exit.
*/
{
	derive_q_K1985_from_x(K, K->q, K->x);
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_x_K1985(
	KConfig K)
/*
** Compute K->x array equivalent to Kondrashov's q.
** K->q is correct upon exit.
*/
{
	derive_q_K1985_from_x(K, K->q, K->x);
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_qp_K1985(
	KConfig K)
/*
** Compute K->qp array equivalent to Kondrashov's q'.
** K->qp is correct upon exit.
*/
{
	derive_q_K1985_from_x(K, K->qp, K->xp);
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_qpp_K1985(
	KConfig K)
/*
** Compute K->qpp array equivalent to the summed result of
** Kondrashov's equation for q".
** K->qpp is correct upon exit.
*/
{
	derive_q_K1985_from_x(K, K->qpp, K->xpp);
}

/* ///////////////////////////////////////////////////////////////////// */
void
derive_Q_K1985(
	KConfig K)
/*
** Compute K->Q array equivalent to Kondrashov's Q.
** K->Q is correct upon exit.
*/
{
	derive_q_K1985_from_x(K, K->Q, K->X);
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_F_female_K1985(
	KConfig K)
/*
** Proportion of population resources allocated to producing all ovules
** Must be called after compute_xp.
*/
{
	KInt k, i, j;
	KScalar sum, grandsum;
	grandsum = 0.0;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			sum = 0.0;
			for (k = 0; k < K->genotypes; k++) {
				sum += K->xp[k][i][j] * (1 + K->y[k] * K->z[k]);
			}
			grandsum += K->qp[i][j] * sum;
		}
	}
	K->F_female = 0.5 * grandsum;
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_F_male_K1985(
	KConfig K)
/*
** Proportion of population resources allocated to producing all pollen
** Must be called after compute_F_female.
*/
{
	K->F_male = 1 - K->F_female;
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_beta_gamma_rho_gk_ak_hk_K1985(
	KConfig K)
/*
** Arrays required by Kondrashov's original formulation of the model.
** Proportion of F_female that is allocated to outcrossed seeds (beta),
**     to apomictic/selfed seeds (gamma) and to pollen for outcrossing
**     (rho), by genotype
** Proportion of resources of kappa-th genotype in G(i,j) allocated
**     to produce outcrossed seeds (gk), apomictic/selfed seeds (ak),
**     and to pollen (hk) from all resources allocated to same function
**     for all individuals in G(i,j)
*/
{
	KInt g, i, j;
	KScalar betasum, gammasum, rhosum;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			betasum = gammasum = rhosum = 0.0;
			for (g = 0; g < K->genotypes; g++) {
				betasum += K->xp[g][i][j] * (1 - K->z[g]);
				gammasum += K->xp[g][i][j] * (K->z[g] + K->y[g] * K->z[g]);
				rhosum += K->xp[g][i][j] * (1 - K->y[g] * K->z[g]);
			}
			K->beta[i][j] = ((0.5 * K->qp[i][j]) / K->F_female) * betasum;
			K->gamma[i][j] = ((0.5 * K->qp[i][j]) / K->F_female) * gammasum;
			K->rho[i][j] = ((0.5 * K->qp[i][j]) / K->F_male) * rhosum;
			for (g = 0; g < K->genotypes; g++) {
				K->gk[g][i][j] = K->xp[g][i][j] * (1 - K->z[g]) / betasum;
				K->ak[g][i][j] =
					K->xp[g][i][j] * (K->z[g] + K->y[g] * K->z[g]) / gammasum;
				K->hk[g][i][j] =
					K->xp[g][i][j] * (1 - K->y[g] * K->z[g]) / rhosum;
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_progeny_K1985(
	KConfig K)
/*
** Compute proportions of genotypes produced by population.
** Must be run after compute_qpp
** This implements Kondrashov 1985, p.640, eq. (2), x"i(k) part
** K->xpp is correct on exit
*/
{
	KInt g, i, j, n, v, l, lam, xi, zeta;
	KScalar ksum, leftsum, rightsum, t1, t2;
	for (g = 0; g < K->genotypes; g++) {
		for (i = 0; i <= K->MI; i++) {
			for (j = 0; j <= K->MJ; j++) {
				leftsum = 0.0;
				rightsum = 0.0;
				for (n = 0; n <= K->MI; n++) {
					for (v = 0; n <= K->MJ; v++) {
						/* implements right side of "=", term left of "+" */
						for (l = 0; l <= K->MI; l++) {
							for (lam = 0; lam <= K->MJ; lam++) {
								t1 = K->beta[n][v] * K->rho[l][lam] * b(i, j,
																		n, v,
																		l,
																		lam);
								ksum = 0.0;
								for (xi = 0; xi < K->genotypes; xi++) {
									for (zeta = 0; zeta < K->genotypes; zeta++) {
										t2 = K->gk[xi][n][v] *
											K->hk[zeta][l][lam];
										ksum +=
											t2 * O_transform(K, g, xi, zeta);
									}
								}
								leftsum += t1 * ksum;
							}
						}
						/* implements right side of "=", term right of "+" */
						if (!K->mating_self && !K->mating_apomixis) {
							rightsum = 0.0;
						} else {
							ksum = 0.0;
							for (xi = 0; xi < K->genotypes; xi++) {
								if (K->mating_apomixis) {
									t1 = A_transform(K, g, xi);
								} else if (K->mating_self) {
									t1 = S_transform(K, g, xi);
								}
								ksum += K->ak[xi][n][v] * t1;
							}
							if (K->mating_apomixis) {
								t1 = a_apomixis(i, j, n, v);
							} else if (K->mating_self) {
								t1 = s_self(i, j, n, v);
							} else {
								nrerror
									("compute_progeny_K1985: unrecognized mating");
							}
							rightsum += (K->gamma[n][v] * t1 * ksum);
						}
					}
				}
				/* final assignment to x"i(g) */
				K->xpp[g][i][j] = (leftsum + rightsum) / K->qpp[i][j];
				/* it might be true that rightsum is not divided; paper is unclear */
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_qpp_K1985(
	KConfig K)
/*
** Add population proportions of selfed and outcrossed zygotes.
** Proportions in K->qppx and K->qpps and K->qppa are summed.
** The algorithm is O(K->MI * K->MJ).
** Must be called after compute_qppx(), compute_qpps(), compute_qppa()
** K->qpp is correct upon exit.
*/
{
	KInt i, j;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			K->qpp[i][j] = K->qpps[i][j] + K->qppa[i][j] + K->qppx[i][j];
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_Q_K1985(
	KConfig K)
/*
** Apply selection using the method appropriate to the model type.
** Kondrashov and Charlesworth et al. apparently do it differently.
** K->Q is correct upon exit.
*/
{
	KInt i, j;
	KScalar k_wm = w_mean(K, K->qpp);
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			switch (K->type) {
			case TYPE_KONDRASHOV:
				K->Q[i][j] = (K->qpp[i][j] * w(K, i, j)) / k_wm;
				break;
			case TYPE_CHARLESWORTH_ET_AL:
				K->Q[i][j] = K->qpp[i][j] * w(K, i, j);
				break;
			default:
				nrerror("compute_Q: unimplemented model type");
				break;
			}
		}
	}
}

/* ///////////////////////////////////////////////////////////////////// */
void
compute_q_nextgen_K1985(
	KConfig K)
/*
** Create reproductive adults that will begin the next generation
** K->q is correct and ready for the next iteration upon exit.
*/
{
	KInt i, j;
	for (i = 0; i <= K->MI; i++) {
		for (j = 0; j <= K->MJ; j++) {
			K->q_prevgen[i][j] = K->q[i][j];
			K->q[i][j] = K->Q[i][j];
		}
	}
}
