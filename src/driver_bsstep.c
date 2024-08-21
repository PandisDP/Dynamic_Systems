/* Driver for routine bsstep */

#include <stdio.h>
#define NRANSI
#include "my_nr.h"
#include "nrutil.h"

#define N 4

double dxsav,*xp,**py;   /* defining declarations */
int kmax,kount;

int nrhs;   /* counts function evaluations */

void derivs(double x,double y[],double dydx[], double C[], double AUX[])
{
	nrhs++;
	dydx[1] = -y[2];
	dydx[2]=y[1]-(1.0/x)*y[2];
	dydx[3]=y[2]-(2.0/x)*y[3];
	dydx[4]=y[3]-(3.0/x)*y[4];
}

int main(void)
{
	int i,nbad,nok;
	double eps=1.0e-4,h1=0.1,hmin=0.0,x1=1.0,x2=10.0,*y;
   double *C, *AUX;
	y=dvector(1,N);
	xp=dvector(1,200);
	py=dmatrix(1,10,1,200);
   C = dvector(1,20);
   AUX = dvector(1,20);
	y[1]=(double)bessj0(x1);
	y[2]=(double)bessj1(x1);
	y[3]=(double)bessj(2,x1);
	y[4]=(double)bessj(3,x1);
	nrhs=0;
	kmax=100;
	dxsav=(x2-x1)/20.0;
	odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep_rz, C, AUX);
	printf("\n%s %13s %3d\n","successful steps:"," ",nok);
	printf("%s %20s %3d\n","bad steps:"," ",nbad);
	printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
	printf("\n%s %3d\n","stored intermediate values:    ",kount);
	printf("\n%8s %18s %15s\n","x","integral","bessj(3,x)");
	for (i=1;i<=kount;i++)
		printf("%10.4f %16.6f %14.6f\n",
			xp[i],py[4][i],bessj(3,xp[i]));
	free_dmatrix(py,1,10,1,200);
	free_dvector(xp,1,200);
	free_dvector(y,1,N);
	return 0;
}
#undef NRANSI
/* gcc driver_bsstep.c my_bsstep.c my_mmid.c nrutil.c bessj.c bessj1.c bessj0.c my_odeint.c my_pzextr.c -o ch2 -lm */
