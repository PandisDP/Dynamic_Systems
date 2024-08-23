/*  Program generated by teval1, a translator and LHS/PRCC evaluator written by Jose L. Segovia,  2002-03, for the Kirschner Lab.  */ 
#include <stdio.h> 
#include <stdlib.h> 
#include <math.h>  
#include <time.h>  
#include <string.h> 
#include "datatype.h"
#define NRANSI 
#include "my_nr.h"
#include "nrutil.h"
double dxsav,*xp,**py;  /* defining declarations */ 
int kmax,kount; 
int nrhs;   /* counts function evaluations */ 
#define NINC 11
#define N 6
#define NFigs 6
	 int myreturnsyscall = 0;

void derivs(double TimeAct, double y[], double dydx[], double C[], double AUX[]) 
{
	nrhs++; 
	AUX[1]=  C[9];  		 /* trucks_in=w1; */
	AUX[2]=  C[9];  		 /* queue_line_a=w1; */
	AUX[3]=  C[9];  		 /* queue_line_b=w1; */
	AUX[4]=  C[9];  		 /* flow_trucks=w1; */
}
void load_constants(double C[]) {  
	 /* Filling initial conditions */ 
	 C[1] =  10 ; 		 /* trucks_init=10 5 10 15 $U ; */
	 C[2] =  50 ; 		 /* percentage_dis=50 30 50 90 $U ; */
	 C[3] =  2 ; 		 /* distance_speed=2 2 6 10 $U ; */
	 C[4] =  1 ; 		 /* flow_left=1 1 10 50 $U ; */
	 C[5] =  1 ; 		 /* flow_right=1 1 10 50 $U ; */
	 C[6] =  5 ; 		 /* flow_up=5 1 10 50 $U ; */
	 C[7] =  5 ; 		 /* flow_down=5 1 10 50 $U ; */
	 C[8] =  0.1 ; 		 /* percentage_minimum_dist=0.1 0.1 0.3 0.5 $U ; */
	 C[9] =  0 ; 		 /* w1=0 ; */
	 } 
/* End of loading constants */ 
void load_auxlevel(double AUX[], double y[], double C[]) {
	 /* Filling leval and aux variables with initial conditions */ 
	 AUX[1] =  0.0; 	 /* trucks_in = */ 
	 AUX[2] =  0.0; 	 /* queue_line_a = */ 
	 AUX[3] =  0.0; 	 /* queue_line_b = */ 
	 AUX[4] =  0.0; 	 /* flow_trucks = */ 
} 


void update_variables_lhs(int c, double C[]) { 
	 /* Updates the following variables at point given by c */ 
	 
	 if (c == 0)   { 
 		 C[1] = 10.75; 
 		 C[2] = 72.9; 
 		 C[3] = 5.24; 
 		 C[4] = 38.485; 
 		 C[5] = 33.095; 
 		 C[6] = 21.825; 
 		 C[7] = 9.085; 
 		 C[8] = 0.41; 
	 } 
	 if (c == 1)   { 
 		 C[1] = 12.05; 
 		 C[2] = 45.3; 
 		 C[3] = 3.16; 
 		 C[4] = 17.415; 
 		 C[5] = 14.475; 
 		 C[6] = 28.195; 
 		 C[7] = 38.485; 
 		 C[8] = 0.322; 
	 } 
	 if (c == 2)   { 
 		 C[1] = 7.35; 
 		 C[2] = 51.3; 
 		 C[3] = 8.04; 
 		 C[4] = 13.005; 
 		 C[5] = 34.565; 
 		 C[6] = 19.375; 
 		 C[7] = 24.275; 
 		 C[8] = 0.418; 
	 } 
	 if (c == 3)   { 
 		 C[1] = 6.65; 
 		 C[2] = 50.7; 
 		 C[3] = 5.88; 
 		 C[4] = 14.965; 
 		 C[5] = 25.255; 
 		 C[6] = 5.165; 
 		 C[7] = 45.345; 
 		 C[8] = 0.358; 
	 } 
	 if (c == 4)   { 
 		 C[1] = 8.15; 
 		 C[2] = 86.7; 
 		 C[3] = 3.8; 
 		 C[4] = 9.575; 
 		 C[5] = 22.315; 
 		 C[6] = 40.445; 
 		 C[7] = 20.355; 
 		 C[8] = 0.21; 
	 } 
	 if (c == 5)   { 
 		 C[1] = 5.35; 
 		 C[2] = 39.9; 
 		 C[3] = 8.2; 
 		 C[4] = 21.825; 
 		 C[5] = 5.165; 
 		 C[6] = 2.225; 
 		 C[7] = 4.675; 
 		 C[8] = 0.254; 
	 } 
	 if (c == 6)   { 
 		 C[1] = 7.55; 
 		 C[2] = 51.9; 
 		 C[3] = 7.88; 
 		 C[4] = 49.265; 
 		 C[5] = 38.975; 
 		 C[6] = 48.285; 
 		 C[7] = 10.065; 
 		 C[8] = 0.434; 
	 } 
	 if (c == 7)   { 
 		 C[1] = 14.45; 
 		 C[2] = 44.1; 
 		 C[3] = 9.24; 
 		 C[4] = 34.565; 
 		 C[5] = 6.145; 
 		 C[6] = 37.995; 
 		 C[7] = 33.095; 
 		 C[8] = 0.334; 
	 } 
	 if (c == 8)   { 
 		 C[1] = 5.85; 
 		 C[2] = 75.9; 
 		 C[3] = 4.04; 
 		 C[4] = 39.955; 
 		 C[5] = 48.285; 
 		 C[6] = 18.395; 
 		 C[7] = 26.725; 
 		 C[8] = 0.194; 
	 } 
	 if (c == 9)   { 
 		 C[1] = 9.15; 
 		 C[2] = 63.9; 
 		 C[3] = 9.16; 
 		 C[4] = 20.845; 
 		 C[5] = 42.895; 
 		 C[6] = 33.585; 
 		 C[7] = 32.605; 
 		 C[8] = 0.166; 
	 } 
	 if (c == 10)   { 
 		 C[1] = 9.95; 
 		 C[2] = 89.7; 
 		 C[3] = 8.28; 
 		 C[4] = 12.025; 
 		 C[5] = 3.695; 
 		 C[6] = 46.325; 
 		 C[7] = 36.035; 
 		 C[8] = 0.154; 
	 } 
	 if (c == 11)   { 
 		 C[1] = 6.55; 
 		 C[2] = 85.5; 
 		 C[3] = 6.44; 
 		 C[4] = 22.805; 
 		 C[5] = 9.085; 
 		 C[6] = 36.525; 
 		 C[7] = 16.925; 
 		 C[8] = 0.162; 
	 } 
	 if (c == 12)   { 
 		 C[1] = 13.65; 
 		 C[2] = 52.5; 
 		 C[3] = 2.28; 
 		 C[4] = 43.875; 
 		 C[5] = 27.705; 
 		 C[6] = 17.415; 
 		 C[7] = 32.115; 
 		 C[8] = 0.282; 
	 } 
	 if (c == 13)   { 
 		 C[1] = 7.25; 
 		 C[2] = 72.3; 
 		 C[3] = 3.72; 
 		 C[4] = 21.335; 
 		 C[5] = 16.925; 
 		 C[6] = 31.625; 
 		 C[7] = 14.475; 
 		 C[8] = 0.466; 
	 } 
	 if (c == 14)   { 
 		 C[1] = 8.55; 
 		 C[2] = 38.1; 
 		 C[3] = 2.36; 
 		 C[4] = 43.385; 
 		 C[5] = 41.425; 
 		 C[6] = 24.275; 
 		 C[7] = 34.075; 
 		 C[8] = 0.45; 
	 } 
	 if (c == 15)   { 
 		 C[1] = 8.75; 
 		 C[2] = 61.5; 
 		 C[3] = 3.64; 
 		 C[4] = 1.735; 
 		 C[5] = 16.435; 
 		 C[6] = 32.115; 
 		 C[7] = 39.465; 
 		 C[8] = 0.474; 
	 } 
	 if (c == 16)   { 
 		 C[1] = 6.95; 
 		 C[2] = 47.1; 
 		 C[3] = 6.36; 
 		 C[4] = 42.405; 
 		 C[5] = 10.555; 
 		 C[6] = 28.685; 
 		 C[7] = 27.705; 
 		 C[8] = 0.394; 
	 } 
	 if (c == 17)   { 
 		 C[1] = 14.75; 
 		 C[2] = 68.1; 
 		 C[3] = 6.84; 
 		 C[4] = 5.655; 
 		 C[5] = 4.185; 
 		 C[6] = 32.605; 
 		 C[7] = 15.455; 
 		 C[8] = 0.498; 
	 } 
	 if (c == 18)   { 
 		 C[1] = 8.05; 
 		 C[2] = 65.1; 
 		 C[3] = 7; 
 		 C[4] = 18.885; 
 		 C[5] = 39.955; 
 		 C[6] = 34.075; 
 		 C[7] = 49.755; 
 		 C[8] = 0.426; 
	 } 
	 if (c == 19)   { 
 		 C[1] = 7.45; 
 		 C[2] = 43.5; 
 		 C[3] = 4.84; 
 		 C[4] = 33.585; 
 		 C[5] = 43.385; 
 		 C[6] = 35.545; 
 		 C[7] = 49.265; 
 		 C[8] = 0.134; 
	 } 
	 if (c == 20)   { 
 		 C[1] = 12.65; 
 		 C[2] = 81.9; 
 		 C[3] = 2.12; 
 		 C[4] = 41.425; 
 		 C[5] = 19.865; 
 		 C[6] = 3.205; 
 		 C[7] = 46.815; 
 		 C[8] = 0.146; 
	 } 
	 if (c == 21)   { 
 		 C[1] = 11.95; 
 		 C[2] = 32.1; 
 		 C[3] = 2.68; 
 		 C[4] = 30.645; 
 		 C[5] = 7.615; 
 		 C[6] = 44.855; 
 		 C[7] = 35.055; 
 		 C[8] = 0.486; 
	 } 
	 if (c == 22)   { 
 		 C[1] = 7.15; 
 		 C[2] = 32.7; 
 		 C[3] = 2.6; 
 		 C[4] = 16.925; 
 		 C[5] = 37.505; 
 		 C[6] = 29.175; 
 		 C[7] = 47.305; 
 		 C[8] = 0.242; 
	 } 
	 if (c == 23)   { 
 		 C[1] = 14.95; 
 		 C[2] = 49.5; 
 		 C[3] = 3.96; 
 		 C[4] = 25.255; 
 		 C[5] = 45.835; 
 		 C[6] = 20.355; 
 		 C[7] = 25.745; 
 		 C[8] = 0.446; 
	 } 
	 if (c == 24)   { 
 		 C[1] = 14.05; 
 		 C[2] = 54.9; 
 		 C[3] = 2.52; 
 		 C[4] = 2.225; 
 		 C[5] = 26.725; 
 		 C[6] = 13.985; 
 		 C[7] = 9.575; 
 		 C[8] = 0.458; 
	 } 
	 if (c == 25)   { 
 		 C[1] = 7.85; 
 		 C[2] = 57.9; 
 		 C[3] = 3.56; 
 		 C[4] = 49.755; 
 		 C[5] = 29.175; 
 		 C[6] = 44.365; 
 		 C[7] = 44.855; 
 		 C[8] = 0.102; 
	 } 
	 if (c == 26)   { 
 		 C[1] = 9.45; 
 		 C[2] = 48.3; 
 		 C[3] = 7.72; 
 		 C[4] = 18.395; 
 		 C[5] = 24.275; 
 		 C[6] = 12.025; 
 		 C[7] = 4.185; 
 		 C[8] = 0.49; 
	 } 
	 if (c == 27)   { 
 		 C[1] = 12.45; 
 		 C[2] = 46.5; 
 		 C[3] = 9.96; 
 		 C[4] = 26.235; 
 		 C[5] = 7.125; 
 		 C[6] = 47.305; 
 		 C[7] = 41.915; 
 		 C[8] = 0.13; 
	 } 
	 if (c == 28)   { 
 		 C[1] = 9.35; 
 		 C[2] = 59.1; 
 		 C[3] = 8.92; 
 		 C[4] = 23.785; 
 		 C[5] = 11.535; 
 		 C[6] = 24.765; 
 		 C[7] = 48.775; 
 		 C[8] = 0.23; 
	 } 
	 if (c == 29)   { 
 		 C[1] = 10.85; 
 		 C[2] = 45.9; 
 		 C[3] = 7.4; 
 		 C[4] = 17.905; 
 		 C[5] = 28.195; 
 		 C[6] = 22.315; 
 		 C[7] = 23.785; 
 		 C[8] = 0.342; 
	 } 
	 if (c == 30)   { 
 		 C[1] = 9.75; 
 		 C[2] = 33.3; 
 		 C[3] = 3.32; 
 		 C[4] = 14.475; 
 		 C[5] = 28.685; 
 		 C[6] = 23.785; 
 		 C[7] = 26.235; 
 		 C[8] = 0.198; 
	 } 
	 if (c == 31)   { 
 		 C[1] = 13.95; 
 		 C[2] = 35.7; 
 		 C[3] = 5.48; 
 		 C[4] = 27.705; 
 		 C[5] = 18.395; 
 		 C[6] = 27.705; 
 		 C[7] = 7.615; 
 		 C[8] = 0.258; 
	 } 
	 if (c == 32)   { 
 		 C[1] = 5.65; 
 		 C[2] = 48.9; 
 		 C[3] = 6.28; 
 		 C[4] = 31.625; 
 		 C[5] = 46.815; 
 		 C[6] = 9.085; 
 		 C[7] = 34.565; 
 		 C[8] = 0.354; 
	 } 
	 if (c == 33)   { 
 		 C[1] = 9.85; 
 		 C[2] = 69.9; 
 		 C[3] = 3.4; 
 		 C[4] = 16.435; 
 		 C[5] = 31.625; 
 		 C[6] = 4.675; 
 		 C[7] = 21.825; 
 		 C[8] = 0.31; 
	 } 
	 if (c == 34)   { 
 		 C[1] = 5.45; 
 		 C[2] = 62.7; 
 		 C[3] = 6.92; 
 		 C[4] = 2.715; 
 		 C[5] = 44.365; 
 		 C[6] = 37.015; 
 		 C[7] = 3.205; 
 		 C[8] = 0.438; 
	 } 
	 if (c == 35)   { 
 		 C[1] = 5.25; 
 		 C[2] = 67.5; 
 		 C[3] = 8.6; 
 		 C[4] = 6.635; 
 		 C[5] = 3.205; 
 		 C[6] = 8.105; 
 		 C[7] = 47.795; 
 		 C[8] = 0.19; 
	 } 
	 if (c == 36)   { 
 		 C[1] = 11.85; 
 		 C[2] = 63.3; 
 		 C[3] = 3.48; 
 		 C[4] = 10.555; 
 		 C[5] = 40.445; 
 		 C[6] = 13.005; 
 		 C[7] = 42.405; 
 		 C[8] = 0.246; 
	 } 
	 if (c == 37)   { 
 		 C[1] = 10.25; 
 		 C[2] = 34.5; 
 		 C[3] = 6.12; 
 		 C[4] = 29.665; 
 		 C[5] = 2.225; 
 		 C[6] = 20.845; 
 		 C[7] = 17.415; 
 		 C[8] = 0.262; 
	 } 
	 if (c == 38)   { 
 		 C[1] = 7.05; 
 		 C[2] = 75.3; 
 		 C[3] = 9; 
 		 C[4] = 37.995; 
 		 C[5] = 21.335; 
 		 C[6] = 7.125; 
 		 C[7] = 15.945; 
 		 C[8] = 0.39; 
	 } 
	 if (c == 39)   { 
 		 C[1] = 9.65; 
 		 C[2] = 41.1; 
 		 C[3] = 9.32; 
 		 C[4] = 15.455; 
 		 C[5] = 13.985; 
 		 C[6] = 41.915; 
 		 C[7] = 40.935; 
 		 C[8] = 0.298; 
	 } 
	 if (c == 40)   { 
 		 C[1] = 11.35; 
 		 C[2] = 59.7; 
 		 C[3] = 8.12; 
 		 C[4] = 25.745; 
 		 C[5] = 23.295; 
 		 C[6] = 47.795; 
 		 C[7] = 20.845; 
 		 C[8] = 0.326; 
	 } 
	 if (c == 41)   { 
 		 C[1] = 8.35; 
 		 C[2] = 60.3; 
 		 C[3] = 5.56; 
 		 C[4] = 32.605; 
 		 C[5] = 13.495; 
 		 C[6] = 43.385; 
 		 C[7] = 13.005; 
 		 C[8] = 0.182; 
	 } 
	 if (c == 42)   { 
 		 C[1] = 7.95; 
 		 C[2] = 78.9; 
 		 C[3] = 5.08; 
 		 C[4] = 45.345; 
 		 C[5] = 21.825; 
 		 C[6] = 30.645; 
 		 C[7] = 14.965; 
 		 C[8] = 0.398; 
	 } 
	 if (c == 43)   { 
 		 C[1] = 10.15; 
 		 C[2] = 56.1; 
 		 C[3] = 9.08; 
 		 C[4] = 40.445; 
 		 C[5] = 15.945; 
 		 C[6] = 8.595; 
 		 C[7] = 19.865; 
 		 C[8] = 0.106; 
	 } 
	 if (c == 44)   { 
 		 C[1] = 13.45; 
 		 C[2] = 83.7; 
 		 C[3] = 6.52; 
 		 C[4] = 6.145; 
 		 C[5] = 8.595; 
 		 C[6] = 41.425; 
 		 C[7] = 6.635; 
 		 C[8] = 0.346; 
	 } 
	 if (c == 45)   { 
 		 C[1] = 10.35; 
 		 C[2] = 42.9; 
 		 C[3] = 6.6; 
 		 C[4] = 5.165; 
 		 C[5] = 17.905; 
 		 C[6] = 30.155; 
 		 C[7] = 10.555; 
 		 C[8] = 0.362; 
	 } 
	 if (c == 46)   { 
 		 C[1] = 13.35; 
 		 C[2] = 35.1; 
 		 C[3] = 2.76; 
 		 C[4] = 1.245; 
 		 C[5] = 32.115; 
 		 C[6] = 17.905; 
 		 C[7] = 18.395; 
 		 C[8] = 0.294; 
	 } 
	 if (c == 47)   { 
 		 C[1] = 7.65; 
 		 C[2] = 57.3; 
 		 C[3] = 7.96; 
 		 C[4] = 13.495; 
 		 C[5] = 30.155; 
 		 C[6] = 34.565; 
 		 C[7] = 23.295; 
 		 C[8] = 0.15; 
	 } 
	 if (c == 48)   { 
 		 C[1] = 14.25; 
 		 C[2] = 53.7; 
 		 C[3] = 5.16; 
 		 C[4] = 3.205; 
 		 C[5] = 15.455; 
 		 C[6] = 1.245; 
 		 C[7] = 35.545; 
 		 C[8] = 0.11; 
	 } 
	 if (c == 49)   { 
 		 C[1] = 10.45; 
 		 C[2] = 76.5; 
 		 C[3] = 5.4; 
 		 C[4] = 36.035; 
 		 C[5] = 37.015; 
 		 C[6] = 42.895; 
 		 C[7] = 8.105; 
 		 C[8] = 0.314; 
	 } 
	 if (c == 50)   { 
 		 C[1] = 9.55; 
 		 C[2] = 89.1; 
 		 C[3] = 6.04; 
 		 C[4] = 48.775; 
 		 C[5] = 42.405; 
 		 C[6] = 2.715; 
 		 C[7] = 36.525; 
 		 C[8] = 0.442; 
	 } 
	 if (c == 51)   { 
 		 C[1] = 12.95; 
 		 C[2] = 79.5; 
 		 C[3] = 8.84; 
 		 C[4] = 26.725; 
 		 C[5] = 20.845; 
 		 C[6] = 16.925; 
 		 C[7] = 5.655; 
 		 C[8] = 0.482; 
	 } 
	 if (c == 52)   { 
 		 C[1] = 8.95; 
 		 C[2] = 53.1; 
 		 C[3] = 2.04; 
 		 C[4] = 29.175; 
 		 C[5] = 8.105; 
 		 C[6] = 40.935; 
 		 C[7] = 46.325; 
 		 C[8] = 0.43; 
	 } 
	 if (c == 53)   { 
 		 C[1] = 12.55; 
 		 C[2] = 58.5; 
 		 C[3] = 7.32; 
 		 C[4] = 7.125; 
 		 C[5] = 23.785; 
 		 C[6] = 39.465; 
 		 C[7] = 13.495; 
 		 C[8] = 0.186; 
	 } 
	 if (c == 54)   { 
 		 C[1] = 5.55; 
 		 C[2] = 80.7; 
 		 C[3] = 4.92; 
 		 C[4] = 27.215; 
 		 C[5] = 35.055; 
 		 C[6] = 6.145; 
 		 C[7] = 48.285; 
 		 C[8] = 0.406; 
	 } 
	 if (c == 55)   { 
 		 C[1] = 11.65; 
 		 C[2] = 66.3; 
 		 C[3] = 5; 
 		 C[4] = 22.315; 
 		 C[5] = 45.345; 
 		 C[6] = 22.805; 
 		 C[7] = 38.975; 
 		 C[8] = 0.234; 
	 } 
	 if (c == 56)   { 
 		 C[1] = 13.15; 
 		 C[2] = 80.1; 
 		 C[3] = 7.48; 
 		 C[4] = 30.155; 
 		 C[5] = 4.675; 
 		 C[6] = 46.815; 
 		 C[7] = 21.335; 
 		 C[8] = 0.454; 
	 } 
	 if (c == 57)   { 
 		 C[1] = 6.85; 
 		 C[2] = 36.3; 
 		 C[3] = 5.8; 
 		 C[4] = 37.505; 
 		 C[5] = 19.375; 
 		 C[6] = 31.135; 
 		 C[7] = 28.195; 
 		 C[8] = 0.226; 
	 } 
	 if (c == 58)   { 
 		 C[1] = 10.65; 
 		 C[2] = 44.7; 
 		 C[3] = 9.72; 
 		 C[4] = 11.535; 
 		 C[5] = 1.245; 
 		 C[6] = 4.185; 
 		 C[7] = 2.715; 
 		 C[8] = 0.35; 
	 } 
	 if (c == 59)   { 
 		 C[1] = 13.85; 
 		 C[2] = 83.1; 
 		 C[3] = 5.72; 
 		 C[4] = 28.685; 
 		 C[5] = 41.915; 
 		 C[6] = 26.725; 
 		 C[7] = 22.315; 
 		 C[8] = 0.462; 
	 } 
	 if (c == 60)   { 
 		 C[1] = 11.45; 
 		 C[2] = 60.9; 
 		 C[3] = 7.56; 
 		 C[4] = 10.065; 
 		 C[5] = 33.585; 
 		 C[6] = 6.635; 
 		 C[7] = 11.045; 
 		 C[8] = 0.386; 
	 } 
	 if (c == 61)   { 
 		 C[1] = 10.95; 
 		 C[2] = 71.7; 
 		 C[3] = 4.76; 
 		 C[4] = 23.295; 
 		 C[5] = 10.065; 
 		 C[6] = 19.865; 
 		 C[7] = 42.895; 
 		 C[8] = 0.158; 
	 } 
	 if (c == 62)   { 
 		 C[1] = 5.15; 
 		 C[2] = 66.9; 
 		 C[3] = 5.96; 
 		 C[4] = 47.305; 
 		 C[5] = 17.415; 
 		 C[6] = 16.435; 
 		 C[7] = 41.425; 
 		 C[8] = 0.366; 
	 } 
	 if (c == 63)   { 
 		 C[1] = 6.75; 
 		 C[2] = 77.1; 
 		 C[3] = 9.4; 
 		 C[4] = 9.085; 
 		 C[5] = 46.325; 
 		 C[6] = 10.555; 
 		 C[7] = 6.145; 
 		 C[8] = 0.29; 
	 } 
	 if (c == 64)   { 
 		 C[1] = 11.05; 
 		 C[2] = 42.3; 
 		 C[3] = 6.76; 
 		 C[4] = 48.285; 
 		 C[5] = 1.735; 
 		 C[6] = 45.835; 
 		 C[7] = 12.515; 
 		 C[8] = 0.126; 
	 } 
	 if (c == 65)   { 
 		 C[1] = 12.35; 
 		 C[2] = 74.1; 
 		 C[3] = 3.08; 
 		 C[4] = 45.835; 
 		 C[5] = 26.235; 
 		 C[6] = 14.965; 
 		 C[7] = 1.735; 
 		 C[8] = 0.138; 
	 } 
	 if (c == 66)   { 
 		 C[1] = 8.85; 
 		 C[2] = 77.7; 
 		 C[3] = 3; 
 		 C[4] = 37.015; 
 		 C[5] = 38.485; 
 		 C[6] = 25.745; 
 		 C[7] = 11.535; 
 		 C[8] = 0.222; 
	 } 
	 if (c == 67)   { 
 		 C[1] = 13.75; 
 		 C[2] = 70.5; 
 		 C[3] = 7.24; 
 		 C[4] = 8.105; 
 		 C[5] = 2.715; 
 		 C[6] = 38.485; 
 		 C[7] = 16.435; 
 		 C[8] = 0.218; 
	 } 
	 if (c == 68)   { 
 		 C[1] = 14.85; 
 		 C[2] = 78.3; 
 		 C[3] = 4.6; 
 		 C[4] = 24.765; 
 		 C[5] = 49.755; 
 		 C[6] = 23.295; 
 		 C[7] = 8.595; 
 		 C[8] = 0.494; 
	 } 
	 if (c == 69)   { 
 		 C[1] = 11.75; 
 		 C[2] = 31.5; 
 		 C[3] = 4.36; 
 		 C[4] = 38.975; 
 		 C[5] = 31.135; 
 		 C[6] = 11.045; 
 		 C[7] = 33.585; 
 		 C[8] = 0.278; 
	 } 
	 if (c == 70)   { 
 		 C[1] = 8.25; 
 		 C[2] = 86.1; 
 		 C[3] = 8.44; 
 		 C[4] = 11.045; 
 		 C[5] = 20.355; 
 		 C[6] = 10.065; 
 		 C[7] = 43.875; 
 		 C[8] = 0.414; 
	 } 
	 if (c == 71)   { 
 		 C[1] = 12.25; 
 		 C[2] = 33.9; 
 		 C[3] = 4.28; 
 		 C[4] = 34.075; 
 		 C[5] = 48.775; 
 		 C[6] = 13.495; 
 		 C[7] = 22.805; 
 		 C[8] = 0.274; 
	 } 
	 if (c == 72)   { 
 		 C[1] = 8.65; 
 		 C[2] = 47.7; 
 		 C[3] = 6.2; 
 		 C[4] = 40.935; 
 		 C[5] = 37.995; 
 		 C[6] = 27.215; 
 		 C[7] = 37.505; 
 		 C[8] = 0.214; 
	 } 
	 if (c == 73)   { 
 		 C[1] = 14.55; 
 		 C[2] = 82.5; 
 		 C[3] = 2.2; 
 		 C[4] = 3.695; 
 		 C[5] = 25.745; 
 		 C[6] = 38.975; 
 		 C[7] = 31.625; 
 		 C[8] = 0.27; 
	 } 
	 if (c == 74)   { 
 		 C[1] = 5.75; 
 		 C[2] = 54.3; 
 		 C[3] = 6.68; 
 		 C[4] = 12.515; 
 		 C[5] = 22.805; 
 		 C[6] = 29.665; 
 		 C[7] = 13.985; 
 		 C[8] = 0.286; 
	 } 
	 if (c == 75)   { 
 		 C[1] = 10.05; 
 		 C[2] = 40.5; 
 		 C[3] = 4.2; 
 		 C[4] = 44.365; 
 		 C[5] = 36.525; 
 		 C[6] = 33.095; 
 		 C[7] = 27.215; 
 		 C[8] = 0.382; 
	 } 
	 if (c == 76)   { 
 		 C[1] = 12.85; 
 		 C[2] = 74.7; 
 		 C[3] = 3.88; 
 		 C[4] = 36.525; 
 		 C[5] = 39.465; 
 		 C[6] = 15.455; 
 		 C[7] = 17.905; 
 		 C[8] = 0.206; 
	 } 
	 if (c == 77)   { 
 		 C[1] = 13.05; 
 		 C[2] = 38.7; 
 		 C[3] = 8.76; 
 		 C[4] = 4.185; 
 		 C[5] = 43.875; 
 		 C[6] = 14.475; 
 		 C[7] = 24.765; 
 		 C[8] = 0.25; 
	 } 
	 if (c == 78)   { 
 		 C[1] = 8.45; 
 		 C[2] = 64.5; 
 		 C[3] = 9.56; 
 		 C[4] = 13.985; 
 		 C[5] = 34.075; 
 		 C[6] = 25.255; 
 		 C[7] = 7.125; 
 		 C[8] = 0.422; 
	 } 
	 if (c == 79)   { 
 		 C[1] = 6.15; 
 		 C[2] = 56.7; 
 		 C[3] = 7.08; 
 		 C[4] = 32.115; 
 		 C[5] = 35.545; 
 		 C[6] = 35.055; 
 		 C[7] = 18.885; 
 		 C[8] = 0.47; 
	 } 
	 if (c == 80)   { 
 		 C[1] = 11.15; 
 		 C[2] = 36.9; 
 		 C[3] = 7.16; 
 		 C[4] = 19.865; 
 		 C[5] = 49.265; 
 		 C[6] = 9.575; 
 		 C[7] = 29.175; 
 		 C[8] = 0.306; 
	 } 
	 if (c == 81)   { 
 		 C[1] = 6.35; 
 		 C[2] = 50.1; 
 		 C[3] = 4.68; 
 		 C[4] = 41.915; 
 		 C[5] = 36.035; 
 		 C[6] = 12.515; 
 		 C[7] = 44.365; 
 		 C[8] = 0.114; 
	 } 
	 if (c == 82)   { 
 		 C[1] = 10.55; 
 		 C[2] = 30.9; 
 		 C[3] = 2.44; 
 		 C[4] = 47.795; 
 		 C[5] = 9.575; 
 		 C[6] = 3.695; 
 		 C[7] = 2.225; 
 		 C[8] = 0.302; 
	 } 
	 if (c == 83)   { 
 		 C[1] = 11.25; 
 		 C[2] = 87.9; 
 		 C[3] = 2.92; 
 		 C[4] = 19.375; 
 		 C[5] = 40.935; 
 		 C[6] = 21.335; 
 		 C[7] = 3.695; 
 		 C[8] = 0.378; 
	 } 
	 if (c == 84)   { 
 		 C[1] = 13.25; 
 		 C[2] = 87.3; 
 		 C[3] = 2.84; 
 		 C[4] = 44.855; 
 		 C[5] = 47.795; 
 		 C[6] = 49.755; 
 		 C[7] = 30.155; 
 		 C[8] = 0.402; 
	 } 
	 if (c == 85)   { 
 		 C[1] = 9.05; 
 		 C[2] = 73.5; 
 		 C[3] = 4.52; 
 		 C[4] = 15.945; 
 		 C[5] = 11.045; 
 		 C[6] = 43.875; 
 		 C[7] = 1.245; 
 		 C[8] = 0.33; 
	 } 
	 if (c == 86)   { 
 		 C[1] = 6.25; 
 		 C[2] = 81.3; 
 		 C[3] = 3.24; 
 		 C[4] = 31.135; 
 		 C[5] = 14.965; 
 		 C[6] = 1.735; 
 		 C[7] = 40.445; 
 		 C[8] = 0.37; 
	 } 
	 if (c == 87)   { 
 		 C[1] = 5.05; 
 		 C[2] = 41.7; 
 		 C[3] = 4.12; 
 		 C[4] = 20.355; 
 		 C[5] = 12.025; 
 		 C[6] = 49.265; 
 		 C[7] = 25.255; 
 		 C[8] = 0.178; 
	 } 
	 if (c == 88)   { 
 		 C[1] = 13.55; 
 		 C[2] = 65.7; 
 		 C[3] = 8.36; 
 		 C[4] = 4.675; 
 		 C[5] = 27.215; 
 		 C[6] = 37.505; 
 		 C[7] = 31.135; 
 		 C[8] = 0.374; 
	 } 
	 if (c == 89)   { 
 		 C[1] = 11.55; 
 		 C[2] = 84.9; 
 		 C[3] = 9.48; 
 		 C[4] = 39.465; 
 		 C[5] = 6.635; 
 		 C[6] = 26.235; 
 		 C[7] = 19.375; 
 		 C[8] = 0.17; 
	 } 
	 if (c == 90)   { 
 		 C[1] = 6.45; 
 		 C[2] = 84.3; 
 		 C[3] = 7.64; 
 		 C[4] = 7.615; 
 		 C[5] = 44.855; 
 		 C[6] = 18.885; 
 		 C[7] = 45.835; 
 		 C[8] = 0.174; 
	 } 
	 if (c == 91)   { 
 		 C[1] = 6.05; 
 		 C[2] = 30.3; 
 		 C[3] = 5.64; 
 		 C[4] = 42.895; 
 		 C[5] = 5.655; 
 		 C[6] = 45.345; 
 		 C[7] = 39.955; 
 		 C[8] = 0.266; 
	 } 
	 if (c == 92)   { 
 		 C[1] = 14.35; 
 		 C[2] = 88.5; 
 		 C[3] = 7.8; 
 		 C[4] = 33.095; 
 		 C[5] = 32.605; 
 		 C[6] = 15.945; 
 		 C[7] = 28.685; 
 		 C[8] = 0.478; 
	 } 
	 if (c == 93)   { 
 		 C[1] = 14.15; 
 		 C[2] = 37.5; 
 		 C[3] = 8.52; 
 		 C[4] = 46.815; 
 		 C[5] = 24.765; 
 		 C[6] = 5.655; 
 		 C[7] = 37.015; 
 		 C[8] = 0.318; 
	 } 
	 if (c == 94)   { 
 		 C[1] = 9.25; 
 		 C[2] = 71.1; 
 		 C[3] = 5.32; 
 		 C[4] = 28.195; 
 		 C[5] = 18.885; 
 		 C[6] = 7.615; 
 		 C[7] = 37.995; 
 		 C[8] = 0.338; 
	 } 
	 if (c == 95)   { 
 		 C[1] = 7.75; 
 		 C[2] = 69.3; 
 		 C[3] = 9.64; 
 		 C[4] = 35.055; 
 		 C[5] = 13.005; 
 		 C[6] = 42.405; 
 		 C[7] = 43.385; 
 		 C[8] = 0.122; 
	 } 
	 if (c == 96)   { 
 		 C[1] = 5.95; 
 		 C[2] = 39.3; 
 		 C[3] = 8.68; 
 		 C[4] = 8.595; 
 		 C[5] = 29.665; 
 		 C[6] = 36.035; 
 		 C[7] = 30.645; 
 		 C[8] = 0.202; 
	 } 
	 if (c == 97)   { 
 		 C[1] = 14.65; 
 		 C[2] = 68.7; 
 		 C[3] = 9.88; 
 		 C[4] = 46.325; 
 		 C[5] = 12.515; 
 		 C[6] = 39.955; 
 		 C[7] = 5.165; 
 		 C[8] = 0.142; 
	 } 
	 if (c == 98)   { 
 		 C[1] = 12.75; 
 		 C[2] = 62.1; 
 		 C[3] = 9.8; 
 		 C[4] = 24.275; 
 		 C[5] = 30.645; 
 		 C[6] = 11.535; 
 		 C[7] = 29.665; 
 		 C[8] = 0.238; 
	 } 
	 if (c == 99)   { 
 		 C[1] = 12.15; 
 		 C[2] = 55.5; 
 		 C[3] = 4.44; 
 		 C[4] = 35.545; 
 		 C[5] = 47.305; 
 		 C[6] = 48.775; 
 		 C[7] = 12.025; 
 		 C[8] = 0.118; 
	 } 
 } /*  End of update_variables_lhs  */ 


static int nfigsloaded; 
int load_data_figures(double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, int jj, double y[], double C[], double AUX[], double x, FILE * myfile) { 
	 int j = nfigsloaded; 
	 int i; int errorover=0; 
		  FIGS[0][j] = x; 
		  if (! finite(AUX[1])) { printf("Overflow error on variable trucks_in \n"); errorover++; } 
		  FIGS[1][j] = AUX[1];
		  if (figdata[1].max < FIGS[1][j]) { figdata[1].max =  FIGS[1][j]; } 
		  if (figdata[1].min > FIGS[1][j]) { figdata[1].min =  FIGS[1][j]; } 
		  if (! finite(AUX[2])) { printf("Overflow error on variable queue_line_a \n"); errorover++; } 
		  FIGS[2][j] = AUX[2];
		  if (figdata[2].max < FIGS[2][j]) { figdata[2].max =  FIGS[2][j]; } 
		  if (figdata[2].min > FIGS[2][j]) { figdata[2].min =  FIGS[2][j]; } 
		  if (! finite(AUX[3])) { printf("Overflow error on variable queue_line_b \n"); errorover++; } 
		  FIGS[3][j] = AUX[3];
		  if (figdata[3].max < FIGS[3][j]) { figdata[3].max =  FIGS[3][j]; } 
		  if (figdata[3].min > FIGS[3][j]) { figdata[3].min =  FIGS[3][j]; } 
		  if (! finite(AUX[4])) { printf("Overflow error on variable flow_trucks \n"); errorover++; } 
		  FIGS[4][j] = AUX[4];
		  if (figdata[4].max < FIGS[4][j]) { figdata[4].max =  FIGS[4][j]; } 
		  if (figdata[4].min > FIGS[4][j]) { figdata[4].min =  FIGS[4][j]; } 
	 nfigsloaded ++;
	 if (nfigsloaded >= MAX_FIGS ) { int xx, yy; 
		 for (xx=0;xx<MAX_FIGS;xx++) { 
		 fprintf(myfile,"%f\t",FIGS[0][xx]); 
		  for (yy=1;yy<=4;yy++) 
			  fprintf(myfile, "%.6e\t",FIGS[yy][xx]); 
		 fprintf(myfile,"\n"); 
} nfigsloaded =0; } 
 return errorover; 
} 


int main (void) 
{ 
	 time_t init_time=time(NULL), end_time; double mydifftime;
	 int i,j, resode=0;  int myc=0; 
	 double h, x=1.0, *y, *dydx, *yout, *C, *AUX; 
	 FILE *stream; 
	 char filename[5]; 
	 int nfigs = 4;
	 eqfig figdata[6];
	 double FIGS[6][MAX_FIGS]; 
	 int maxloop =0, nloop =0; int figres=0; 
	 y    =  dvector(1,N);  
	 dydx =  dvector(1,N); 
	 yout =  dvector(1,N);  
	 C   =   dvector(1,NINC); 
	 AUX =   dvector(1,N); 
	 nfigsloaded = 0; 
	 for(i=1;i<=nfigs;i++) { figdata[i].max = -HUGE_VAL; figdata[i].min=HUGE_VAL; }; 
	  strcpy(figdata[1].name,"trucks_in");
	  strcpy(figdata[2].name,"queue_line_a");
	  strcpy(figdata[3].name,"queue_line_b");
	  strcpy(figdata[4].name,"flow_trucks");
	 while (myc < 100) { 
	 load_constants(C); 
	 update_variables_lhs(myc, C); 
	 load_auxlevel(AUX,y,C); 
	 nfigsloaded = 0; 
	 /* ---rk4--- */ 
	 h = (double)0.1;
	 x = (double)0; 
	 maxloop = (int) (3600/h)+1;
	 nloop = 0; 
	 j=0; 
	 printf("Solving run %d \n", myc+1);
	 sprintf(filename, "%04u", myc+1); 
	 stream = fopen(filename,"w"  ); 
	 while(nloop <= maxloop) {
		 if (nloop % (int)ceil(1/h) == 0) {
		 resode=load_data_figures(FIGS, figdata, nfigs, j, y, C, AUX, x,  stream);
		 if (resode > 0) break; 
			 j++; }  
		 derivs(x,y,dydx,C, AUX); 
		 rk4(y,dydx,N,x,h,yout,derivs, C, AUX); 
		 for(i=0; i<N; i++) 
		 y[i] = yout[i]; 
		 nloop++; 
		 x = x + h;
	 }  
	for(i=0; i<nfigsloaded; i++) { 
		 fprintf(stream,"%f\t",FIGS[0][i]);
		 for(j=1; j<=nfigs; j++) fprintf(stream,"%.5e\t",FIGS[j][i]); 
		 fprintf(stream,"\n"); }; 
		 fclose(stream); 
	 if (resode > 0) { printf( "Error solving run %d \n", myc+1); break;	} 
	 myc++; } 
	 WritePlotScript(1, nfigs, 4320, figdata, 0, myc-1);
	 myreturnsyscall = system("rm output1.pdf; gnuplot plot.scr &"); 
	 free_dvector(yout,1,N); 
	 free_dvector(dydx,1,N); 
	 free_dvector(y,1,N); 
	 free_dvector(C,1,NINC); 
	 free_dvector(AUX,1,N); 
	 end_time = time(NULL); 
	 mydifftime = difftime(end_time, init_time); 
	 printf("##  %f seconds \n",mydifftime); 
	 return 0; 
} 
