#include <math.h>
#include "nrutil.h"
#include "my_nr.h"
#include "datatype.h"
#define MAXSTP 1000000 /* it was 10000 */
#define TINY 1.0e-30

extern int kmax, kount;
extern double *xp,**py, dxsav;

int odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], double C[], double AUX[]),
	int (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [], double [], double[]), 
        double [], double []), double C[], double AUX[], double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, stream  myfile)
{
  
  int nstp,i,j, resoper; int resfig;
	double xsav,x,hnext,hdid,h;
	double *yscal,*y,*dydx;
	yscal=dvector(1,nvar);
	y=dvector(1,nvar);
	dydx=dvector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) = kount = 0; j=0;
	for (i=1;i<=nvar;i++) y[i]=ystart[i];
	if (kmax > 0) xsav=x-dxsav*2.0;
	for (nstp=1;nstp<=MAXSTP;nstp++) 
	  {
		 (*derivs)(x,y,dydx, C, AUX);
		 for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		 if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) 
			{
			  xp[++kount]=x;
			  for (i=1;i<=nvar;i++) py[i][kount]=y[i];
			  xsav=x;		
			  resfig=load_data_figures(FIGS, figdata, nfigs, kount, y, C, AUX, x, myfile);
			  if (resfig > 0) /* error */
				 {
					free_dvector(dydx,1,nvar);
					free_dvector(y,1,nvar);
					free_dvector(yscal,1,nvar);
					nrerror( "Overflow error!");
					return 4;
				 }
			}
		 if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		resoper=(*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs, C, AUX);
		if (resoper > 0)
		  {
			 free_dvector(dydx,1,nvar);
			 free_dvector(y,1,nvar);
			 free_dvector(yscal,1,nvar);
			 nrerror("Error in internal operation of evaluator -intode");
			 return 3; 
		  }
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
 			load_data_figures(FIGS, figdata, nfigs, kount, y, C, AUX, x, myfile);
				if (kmax) {
				xp[++kount]=x;
				for (i=1;i<=nvar;i++) py[i][kount]=y[i];
				resfig = load_data_figures(FIGS, figdata, nfigs, kount, y, C, AUX, x, myfile);
				if (resfig > 0) /* error */
			  {
				 free_dvector(dydx,1,nvar);
				 free_dvector(y,1,nvar);
				 free_dvector(yscal,1,nvar);
				 nrerror( "Overflow error!");
				 return 4;
			  }	
				}				
			free_dvector(dydx,1,nvar);
			free_dvector(y,1,nvar);
			free_dvector(yscal,1,nvar);
			return 0;
		}
		if (fabs(hnext) <= hmin) {  
		  free_dvector(dydx,1,nvar);
		  free_dvector(y,1,nvar);
		  free_dvector(yscal,1,nvar);
		  nrerror("Step size too small in odeint");
		  return 1;
		  }
		h=hnext;
	}
	free_dvector(dydx,1,nvar);
	free_dvector(y,1,nvar);
	free_dvector(yscal,1,nvar);
	nrerror("--Too many steps in routine odeint");
	return 2;
}
#undef MAXSTP
#undef TINY
#undef NRANSI
