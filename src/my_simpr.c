#include "nrutil.h"
#include "my_nr.h"

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
			  void (*derivs)(double, double [], double [],double C[], double AUX[]),double C[], double AUX[] )  
{

	int i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;

	indx=ivector(1,n);
	a=dmatrix(1,n,1,n);
	del=dvector(1,n);
	ytemp=dvector(1,n);
	h=htot/nstep;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
		++a[i][i];
	}
	ludcmp(a,n,indx,&d);
	for (i=1;i<=n;i++)
		yout[i]=h*(dydx[i]+h*dfdx[i]);
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	(*derivs)(x,ytemp,yout, C, AUX);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
			yout[i]=h*yout[i]-del[i];
		lubksb(a,n,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		(*derivs)(x,ytemp,yout, C, AUX);
	}
	for (i=1;i<=n;i++)
		yout[i]=h*yout[i]-del[i];
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_dvector(ytemp,1,n);
	free_dvector(del,1,n);
	free_dmatrix(a,1,n,1,n);
	free_ivector(indx,1,n);
}

