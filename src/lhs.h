#include "st.h"
#include "f0.h"

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void FillHypercube(int element, double LHSMatrix[][MAX_LHS_MATRIX], struct_st pArray[], int NumSim);
void FillTriangularHypercube(double min, double peak, double max,int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim) ;
void FillUniformHypercube( double distmin, double distmax, int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim);
void FillIntervalHypercube( double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim ) ;
void FillLogHypercube(double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim ) ;
int rand_int(int n1,  int n2, long *idum);
float ran2(long *idum);
