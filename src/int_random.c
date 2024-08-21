#include <limits.h>
#include <math.h>
#define RNMX_INT 0.99999999999999

float ran2(long *idum);

int rand_int(int lo,  int hi, long *midum)
{ 
  /*  j = (int)n1+(int) (n2*ran2(midum)/(RNMX_INT+1.0));  incorrect */
  return lo + (int)(ran2(midum)*hi) % (hi-lo+1);
}
#undef RNMX_INT
