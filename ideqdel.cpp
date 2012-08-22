
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#define CRIT 0.00000001

FILE *fopen(), *outfile, *outfile2;

double calcd,calch,calcr,s, U, tmp,tmp1,tmp2,tmp3,n, q, H, D, R, mu, h, wbarout,
     wbarself, delta, rstar, r, Hstar, Dstar, Rstar, Hstar2, Dstar2, Rstar2;
long gens;
int done;

main(argc,argv)
 int argc;
 char *argv[];

{
int i;
void calcdel();
void nextgen();
void nrerror();

    //outfile = fopen("outgens", "w");
    //outfile2 = fopen("finalideq", "w");

    if (argc != 5) nrerror("Usage: del r U s h");

    sscanf(argv[1], "%lf", &r);
    sscanf(argv[2], "%lf", &U);
    sscanf(argv[3], "%lf", &s);
    sscanf(argv[4], "%lf", &h);

    if (r > 1 || r < 0) nrerror("Selfing rate out of bounds");
    if (s > 1 || s < 0) nrerror("Selection coefficient out of bounds");
    if (h > 0.5 || h < 0) nrerror("Dominance out of bounds");

	D = 1;
	H = 0;
	R = 0;
	n = 100000000;
	mu = U/(2*n);
	
	/*
	fprintf(outfile, "Gens   D     2H      delta     rstar     meanload   q  \n");

        */

        done = 0;
        gens = 0;

	while (done != 1) {
          gens += 1;
          tmp1 = delta;
          tmp2 = 2*n*H;
          tmp3 = 2*n*R;
	  calcdel();
	  nextgen();
	  /*
	  fprintf(outfile, "%ld %lf %lf %lf %lf %lf %lf\n", gens, D, 2*H,  delta, rstar, 2*n*H, H);
	  */
         if (gens > 20) if (fabs(delta/tmp1-1) < CRIT) if (fabs(tmp2/(2*n*H)-1) < CRIT) if (fabs(tmp3/(2*n*R) -1) < CRIT) done = 1;
        if (gens > 100000) done = 1; 
	}
// fprintf(outfile2, "Gens   D       2H     meanout meanself  delta    rstar   meanload     q\n");
//	  fprintf(outfile2, "%ld %lf %g %lf %lf %lf %lf %lf %lf\n", gens, D, 2*H, wbarout,wbarself,delta, rstar, 2*n*H, H);
        printf("r = %lf, ", r);
        printf("mu = %g, ", mu);
        printf("s = %lf, ", s);
        printf("h = %lf, ", h);
        printf("No. loci = %0.0lf\n", n);
	printf("Gens		%ld\nD		%lf\n2H		%g\nR		%g\nmeanout		%lf\nmeanself	%lf\nmeanfit		%lf\ndelta		%lf  %lf\nrstar		%lf\nmeanhet		%lf  %lf\nmeanhom		%lf  %lf\nq		%g\nCharl. U	%lf    %lf\n", gens, D, 2*H, R, wbarout,wbarself, r*wbarself+(1-r)*wbarout,delta, -2*log(1-delta),rstar, 2*n*H, 2*n*Hstar2, n*R, n*Rstar2, H+R, -2*log(1-delta)/(1-2*h), 4*h*log(1-delta)/(2*h-1));

	calcd = (1-2*mu)/(1+3*mu);
	calch = 2*mu*(1-2*mu)/(1+2*mu-3*mu*mu);
	calcr = mu/(1-mu);
	printf("Calcd = %g, 2CalcH = %g, calcr = %g\n",calcd,2*calch,calcr);
}

void calcdel()

{
  q = H + R + mu*(1-H-R);
  wbarout = exp(n*log((-2*q*s+2*q*q*s)*h-q*q*s+1));
  tmp = (-h*s+h*s*mu*mu+3*s/2-s*mu-0.5*s*mu*mu)*H-s-2*D*mu*h*s+2*D*mu*mu*h*s+1+s*D-s*D*mu*mu;
  wbarself = exp(n*log(tmp));
  delta = 1 - wbarself/wbarout;

  rstar = r*wbarself/(r*wbarself+(1-r)*wbarout);

}


void nextgen()

{
double sum;

  Dstar = rstar*(D+H/2)+(1-rstar)*(1-H-R)*(1-H-R);
  Hstar = rstar*H/2 + (1-rstar)*(H+R)*(1-H-R);
  Rstar = 1-Dstar-2*Hstar;
  
  Dstar2 = (1-mu)*(1-mu)*Dstar;
  Hstar2 = (1-mu)*Hstar + mu*(1-mu)*Dstar;
  Rstar2 = 1-Dstar2-2*Hstar2;
  
  sum = Dstar2+(1-s*h)*2*Hstar2+(1-s)*Rstar2;

  D = Dstar2/sum;
  H = (1-s*h)*Hstar2/sum;
  R = (1-s)*Rstar2/sum;

}

void nrerror(error_text)
char error_text[];

{
  void exit();
  fprintf(stderr,"%s\n",error_text);
  exit(1);
}
