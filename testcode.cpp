#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "K.h"

KScalar f(
	KConfig K,
	KInt i,
	KInt j,
	KInt g);

int
main(
	int argc,
	char *argv[])
{
	KInt i, j;
	KScalar t1;
	KArray a;
	KConfig K = (KConfig) calloc(sizeof(struct struct_KConfig), 1);
	void *aPtr = (void *) calloc(sizeof(KArray), (size_t) 1);
	K->MI = 20;
	K->MJ = 10;
	K->genotypes = 1;
	blah_KArray(K, a, -1.0);
	blah_KArray(K, aPtr, -2.0);
	dump_KArray typedef KScalar(
	*KSelModelPtr) (
	KConfig K,
	KInt i,
	KInt j,
	KInt g);
	KSelModelPtr kp;
	printf("sizeof(KSelModelPtr)=%d\n", sizeof(KSelModelPtr));
	printf("sizeof(kp)=%d\n", sizeof(kp));
	printf("kp=0x08%x\n", kp);
	printf("f=0x08%x\n", f);
	printf("&f=0x08%x\n", &f);
	kp = f;
	printf("kp=0x08%x\n", kp);


	/*
	   KConfig K = new_KConfig(200,40);
	   new_genotypes(K, 1);
	   K->U = 100.0;
	   /* set_debug(K, DEBUG_NORMALIZATION); */
	/* initiate_mut_term(K); */
	/*
	   for (i=0; i <= K->MI; i++) {
	   t1 = exp(-K->U) * pow(K->U, i) / factorial(i);
	   if (K->mut_term[i] != t1) {
	   printf("mut_term[%d] %g != direct %g\n", i, K->mut_term[i], t1);
	   }
	   }
	 */
}

KScalar
f(
	KConfig K,
	KInt i,
	KInt j,
	KInt g)
{
}
