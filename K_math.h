#define     MIN(_a_,_b_)            ((_a_)<(_b_) ? _a_ : _b_)
#define     MAX(_a_,_b_)            ((_a_)>(_b_) ? _a_ : _b_)

#define     pow_half_VALS           (MAX_MI * 4 + 4)
#define     factorial_VALS          (MAX_MI * 2)
#define     binomial_K_VALS         (MAX_MI)
#define     binomial_N_VALS         (binomial_K_VALS * 2)

void init_math(
	void);
KScalar pow_half(
	KInt n);
KScalar lnpow_half(
	KInt n);
KScalar factorial(
	KInt n);
KScalar lnfactorial(
	KInt n);
KScalar binomial(
	KInt n,
	KInt k);
KScalar lnbinomial(
	KInt n,
	KInt k);
