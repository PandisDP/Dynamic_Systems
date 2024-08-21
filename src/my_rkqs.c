#include <math.h>
#include "nrutil.h"
#include "datatype.h"
#include "my_nr.h"
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

int rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
	double yscal[], double *hdid, double *hnext,	void (*derivs)(double, double [], double [], double[], double[]), double C[], double AUX[])

{
 	int i;
	double errmax,h,htemp,xnew,*yerr,*ytemp;
	yerr=dvector(1,n);
	ytemp=dvector(1,n);
	h=htry;
	for (;;) {
		rkck(y,dydx,n,*x,h,ytemp,yerr,derivs, C, AUX);
		errmax=0.0;
		for (i=1;i<=n;i++) errmax=DMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? DMAX(htemp,0.1*h) : DMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) 
		  { free_dvector(ytemp,1,n);
		    free_dvector(yerr,1,n);
			 nrerror("stepsize underflow in rkqs"); return 1; }
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=1;i<=n;i++) y[i]=ytemp[i];
	free_dvector(ytemp,1,n);
	free_dvector(yerr,1,n);
	return 0;
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
/* --------------- */ 
