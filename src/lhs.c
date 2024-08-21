
#include "f0.h"
#include "ran2.h"
#include <string.h>
/* The next three subroutines take the upper and lower limits of the          */ 
/* distributions, divide them into equiprobable regions, sample the midpoints */ 
/* of the regions,  randomize them and put them in the LHC matrix             */ 
 
void FillHypercube(int element, double LHSMatrix[][MAX_LHS_MATRIX], st_struct pArray[], int NumSim ) 
{ 
     int      i,loc; 
     double   temp; 
     if (NumSim >= MAX_LHS_MATRIX) error(" NumSim variable too big, or greather than 20");
     int my_element = pArray[element].lhsrun;
     if (pArray[element].type_range == "U")      
	 FillUniformHypercube( pArray[element].begin_range,pArray[element].end_range, my_element, LHSMatrix, NumSim); 
     else if (pArray[element].type_range == "C")  
       FillIntervalHypercube(pArray[element].begin_range,pArray[element].end_range, my_element, LHSMatrix, NumSim);
     else if (pArray[element].type_range == "T")  
       FillTriangularHypercube(pArray[element].begin_range,pArray[element].mid_range,pArray[element].end_range,my_element, LHSMatrix, NumSim);  
     else if (pArray[element].type_range == "L")  
       FillLogHypercube(pArray[element].begin_range,pArray[element].end_range,my_element, LHSMatrix, NumSim); 
     else  error("This error should not happen ever !!");   
	  
      	      
} 
 
void FillTriangularHypercube(double min, double peak, double max,int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim) 
{ 
  double height, volumeSlice, ximin, m1, m2, m, ximax, b1, b2, b, xintoInverse, discrim, LeftVol; 
  int i, needtoswitch; 
 
  height = 2.0/(max-min); 
  volumeSlice = 1.0/NumSim; 
  ximin = min; 
  needtoswitch = 1;  /* initially assume that we do need to switch from upline to downline */ 
   
  if (peak > min) {  /* so there is an upline */ 
    m1 = height/(peak - min); 
    b1 = -height*min/(peak-min); 
    m = m1; 
    b = b1; 
 
    if (peak == max) /* if peak is at max, no downline, so no need to switch */ 
      needtoswitch = 0; 
  } 
 
  if (peak < max) {  /* there is a downline */ 
    m2 = -height/(max - peak); 
    b2 = height*max/(max-peak); 
 
    if (peak == min) {    /* if peak is at min, no upline */ 
      needtoswitch = 0;   /* so no need to switch */ 
      m = m2;             /* and set m, b for downline */ 
      b = b2; 
    } 
  } 
 
  for(i = 0; i< NumSim; i++) { 
    xintoInverse = m*ximin*ximin/2 + b*ximin + volumeSlice; 
    discrim = b*b + 2*m*xintoInverse; 
    if (discrim < 0) 
      discrim = 0; 
 
    ximax = (-b + sqrt(discrim))/m; 
     
    if (needtoswitch)     /* if there is a need to switch from upline to downline */ 
      if(ximax > peak) {  /* check whether we've passed the peak */ 
	 
	LeftVol = m*peak*peak/2 + b*peak - m*ximin*ximin/2 - b*ximin;  /* vol from ximin to peak */ 
 
	m = m2;           /* switch m, b for downline */ 
	b = b2; 
 
	xintoInverse = m*peak*peak/2 + b*peak + (volumeSlice - LeftVol); 
	ximax = (-b + sqrt(b*b + 2*m*xintoInverse))/m;    /* calculate ximax s.t. volume from peak to ximax = volumeSlice - LeftVol */ 
     
	peak = 2*max; 
      } 
    LHSMatrix[element][i] = (ximin + ximax)/2; 
    ximin = ximax; 
  } 
} 
 
 
 /* This subroutine fills in one row of the LHSMatrix with sampled values       */ 
 /* for the current parameter, which has a uniform distribution. Itakes the     */ 
 /* lower and upper limits of the distribution and determines what the          */ 
 /* midpoint of each of region would be if the distribution were divided into   */ 
 /* NumSim equiprobable regions. These midpoints are the sampled values that    */ 
 /* are initially put into the LHS matrix in order from smallest to largest     */ 
 
void FillUniformHypercube( double distmin, double distmax, int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim)
{
      double incr,min; 
      int    k; 
 
      incr = (distmax-distmin)/NumSim; 
      min=distmin; 
      min += incr/2; 
      for(k = 0; k< NumSim; k++) 
      { 
	 LHSMatrix[element][k] = min; 
	 min += incr; 
      } 
 } 
 
 
 /* This subroutine sampes evenly across the entire interval instead of just */ 
 /* at the midpoints of equiprobable regions                                 */ 
 
 void FillIntervalHypercube( double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim ) 
 { 
      double incr,min; 
      int    k; 
 
      incr = (distmax-distmin)/(NumSim-1);   /* width of each subinterval is max-min / number of simulations */       
      min=distmin; 
      for(k = 0; k< NumSim; k++)   
      { 
	 LHSMatrix[element][k] = min;   /* store the last midpoint value in the LHSMatrix */      
	 min += incr;			     /* next midpoint = current midpoint + size of the interval */	 
      } 
 } 
 
 
 void FillLogHypercube(double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim ) 
 { 
     double increment,interval,min; 
     int k; 
     double logdistmin,logdistmax; 
 
     logdistmax=log10(distmax); 
     logdistmin=log10(distmin); 
     increment = (logdistmax-logdistmin)/(NumSim-1); 
     min=logdistmin; 
     for (k=0;k<NumSim;k++) 
     { 
	 LHSMatrix[element][k] = exp(min*log(10.0));    
	 min += increment;   
     } 
 } 
 
int rand_int(int lo,  int hi, long *midum)
{ 
  /*  j = (int)n1+(int) (n2*ran2(midum)/(RNMX_INT+1.0));  incorrect */
  return lo + (int)(ran2(midum)*hi) % (hi-lo+1);
}

