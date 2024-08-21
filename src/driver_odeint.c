/* Driver for routine odeint */

#include <stdio.h>
#include <math.h>
#define NRANSI
#include "my_nr.h"
#include "nrutil.h"

#define N 4
//typedef struct {char name[50];  float max, min; } eqfig;

double dxsav,*xp,**py;  /* defining declarations */

int kmax,kount;

int nrhs;   /* counts function evaluations */

void load_data_figures(float **FIGS, eqfig figdata[], int nfigs, int j, double y[], double C[], double AUX[], double x)

{

}

void derivs(double x, double y[], double dydx[], double C[], double AUX[])
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
	int NPOINTS = 5000;
	double eps=1.0e-4,h1=0.1,hmin=0.0,x1=1.0,x2=1000.0,*ystart, *C, *AUX;
	float **FIGS;   
	eqfig figdata[20]; int nfigs;
  	ystart=dvector(1,N);
	xp=dvector(1,NPOINTS);
	C =dvector(1,20);
	AUX=dvector(1,20);
	py=dmatrix(1,10,1,NPOINTS);

	ystart[1]=bessj0(x1);
	ystart[2]=bessj1(x1);
	ystart[3]=bessj(2,x1);
	ystart[4]=bessj(3,x1);

	nrhs=0;
	kmax=2000;
	dxsav=(x2-x1)/(float)NPOINTS;
	int res;
   res=odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs, rkqs, C, AUX, FIGS, figdata, nfigs);
	/*	printf("\n%s %13s %3d\n","successful steps:"," ",nok);
	printf("%s %20s %3d\n","bad steps:"," ",nbad);
	printf("%s %9s %3d\n","function evaluations:"," ",nrhs);
	printf("\n%s %3d\n","stored intermediate values:    ",kount);
	printf("\n%8s %18s %15s\n","x","integral","bessj(3,x)"); */
	for (i=1;i<=kount;i++)
		printf("%10.4f %16.6f %14.6f\n",
			xp[i],py[4][i],bessj(3,xp[i]));
	if (res > 0) printf("error in odeint ");
	free_dmatrix(py,1,10,1,NPOINTS);
	free_dvector(xp,1,NPOINTS);
	free_dvector(ystart,1,N);


	return 0;
}
/*  gcc driver_odeint.c my_rkck.c my_rkqs.c my_odeint.c nrutil.c bessj.c bessj1.c bessj0.c -o borra1 -lm */
#undef NRANSI
