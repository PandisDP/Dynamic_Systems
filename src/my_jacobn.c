#include <math.h>
#include "my_nr.h"
#include "nrutil.h"

#define EPSJ 1e-6
	

/* void fdjac(int n, float x[], float fvec[], float **df,void (*vecfunc)(int, float [], float [])) */
void jacobn(double n1, double x[], double fvec[], double **df, int n,  double C[], double AUX[])
{
   
   int i,j;
   double h,temp,*f;
   f=dvector(1,n);
   for (j=1;j<=n;j++)
     {
	
        temp=x[j];
        h=EPSJ*fabs(temp);
        if (h == 0.0) h=EPSJ;
        x[j]=temp+h;
        h=x[j]-temp;
        derivs(n,x,f,C,AUX);
        x[j]=temp;
        for (i=1;i<=n;i++) df[i][j]=(f[i]-fvec[i])/h;
     }
   
   free_dvector(f,1,n);
}


