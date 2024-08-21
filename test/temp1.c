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
#define NINC 23
#define N 25
#define NFigs 9
	 int myreturnsyscall = 0;

void derivs(double TimeAct, double y[], double dydx[], double C[], double AUX[]) 
{
	nrhs++; 
	AUX[1]=  y[13]*C[8];  		 /* Nacimientos=Adultos*tasadenacimientos; */
	AUX[2]=  y[12]/C[10];  		 /* Madurez=Jovenes/periodomadurez; */
	AUX[3]=  y[12]*C[9];  		 /* MuerteJovenes=Jovenes*tasademuertejovenes; */
	AUX[4]=  y[13]*C[11];  		 /* MuerteAdultos=Adultos*tasademuerteadultos; */
	AUX[5]=  y[13]/C[12];  		 /* Vejez=Adultos/periodmadurez; */
	AUX[6]=  y[14]/C[13];  		 /* MuerteNatural=Ancianos/periodofinal; */
	AUX[7]=  (fmod(TimeAct,C[14])==0.0)?C[18]:0.0;  		 /* MinasNuevas=(fmod(TimeAct,tasadenuevasminas)==0.0)?numerodenuevasminas:0.0; */
	AUX[8]=  (fmod(TimeAct,C[15])==0.0)?1.0:0.0;  		 /* param1=(fmod(TimeAct,tasacierreminas)==0.0)?1.0:0.0; */
	AUX[9]=  (y[16]>0)?1.0:0.0;  		 /* param2=(Minas>0)?1.0:0.0; */
	AUX[10]=  (TimeAct>(C[16]))?C[17]:0.0;  		 /* param3=(TimeAct>(inicioclosemina))?numberofclosemine:0.0; */
	AUX[11]=  AUX[8]*AUX[9]*AUX[10];  		 /* CierreMinas=param1*param2*param3; */
	dydx[12]=  AUX[1]-AUX[2]-AUX[3];  		 /* Jovenes=Nacimientos-Madurez-MuerteJovenes; */
	dydx[13]=  AUX[2]-AUX[4]-AUX[5];  		 /* Adultos=Madurez-MuerteAdultos-Vejez; */
	dydx[14]=  AUX[5]-AUX[6];  		 /* Ancianos=Vejez-MuerteNatural; */
	AUX[15]=  y[12]+y[13]+y[14];  		 /* Poblacion=Jovenes+Adultos+Ancianos; */
	dydx[16]=  AUX[7]-AUX[11];  		 /* Minas=MinasNuevas-CierreMinas; */
	AUX[17]=  AUX[15]*C[21];  		 /* Presupuesto=Poblacion*presupuestoporpersona; */
	AUX[18]=  C[19]*C[20]*1000*y[16]*0.07-AUX[17];  		 /* deltaingresos=cobreproducido*preciocobre*1000*Minas*0.07-Presupuesto; */
	AUX[19]=  (AUX[18]<0)?-AUX[18]:0;  		 /* deltadeudas=(deltaingresos<0)?-deltaingresos:0; */
	dydx[20]=  AUX[18];  		 /* Ingresos=deltaingresos; */
	dydx[21]=  AUX[19];  		 /* Deudas=deltadeudas; */
	AUX[22]=  (y[20]>0)?y[20]/AUX[15]:0;  		 /* PBIpc=(Ingresos>0)?Ingresos/Poblacion:0; */
	AUX[23]=  (y[21]>0)?y[21]/AUX[15]:0;  		 /* DEUpc=(Deudas>0)?Deudas/Poblacion:0; */
}
void load_constants(double C[]) {  
	 /* Filling initial conditions */ 
	 C[1] =  8.80615e+06 ; 		 /* Jovenes=8.80615e+06 ; */
	 C[2] =  2.27175e+07 ; 		 /* Adultos=2.27175e+07 ; */
	 C[3] =  2.99434e+06 ; 		 /* Ancianos=2.99434e+06 ; */
	 C[4] =  63 ; 		 /* Minas=63 ; */
	 C[5] =  16961e+06 ; 		 /* Ingresos=16961e+06 ; */
	 C[6] =  7600e+06 ; 		 /* Presupuesto=7600e+06 ; */
	 C[7] =  0 ; 		 /* Deudas=0 ; */
	 C[8] =  0.031 ; 		 /* tasadenacimientos=0.031 0.0149 0.0186 0.0546 $U ; */
	 C[9] =  0.0095 ; 		 /* tasademuertejovenes=0.0095 0.0073 0.0093 0.0095 $U ; */
	 C[10] =  15 ; 		 /* periodomadurez=15 ; */
	 C[11] =  0.003 ; 		 /* tasademuerteadultos=0.003 0.002 0.00213 0.003 $U ; */
	 C[12] =  50 ; 		 /* periodmadurez=50 ; */
	 C[13] =  15 ; 		 /* periodofinal=15 ; */
	 C[14] =  20 ; 		 /* tasadenuevasminas=20 10 20 50 $U ; */
	 C[15] =  40 ; 		 /* tasacierreminas=40 20 40 60 $U ; */
	 C[16] =  60 ; 		 /* inicioclosemina=60 30 50 70 $U ; */
	 C[17] =  10 ; 		 /* numberofclosemine=10 3 5 10 $U ; */
	 C[18] =  2 ; 		 /* numerodenuevasminas=2 3 5 10 $U ; */
	 C[19] =  400 ; 		 /* cobreproducido=400 100 400 600 $U ; */
	 C[20] =  9840 ; 		 /* preciocobre=9840 9000 10000 12000 $U ; */
	 C[21] =  400 ; 		 /* presupuestoporpersona=400 100 400 700 $U ; */
	 } 
/* End of loading constants */ 
void load_auxlevel(double AUX[], double y[], double C[]) {
	 /* Filling leval and aux variables with initial conditions */ 
	 AUX[1] =  0.0; 	 /* Nacimientos = */ 
	 AUX[2] =  0.0; 	 /* Madurez = */ 
	 AUX[3] =  0.0; 	 /* MuerteJovenes = */ 
	 AUX[4] =  0.0; 	 /* MuerteAdultos = */ 
	 AUX[5] =  0.0; 	 /* Vejez = */ 
	 AUX[6] =  0.0; 	 /* MuerteNatural = */ 
	 AUX[7] =  0.0; 	 /* MinasNuevas = */ 
	 AUX[8] =  0.0; 	 /* param1 = */ 
	 AUX[9] =  0.0; 	 /* param2 = */ 
	 AUX[10] =  0.0; 	 /* param3 = */ 
	 AUX[11] =  0.0; 	 /* CierreMinas = */ 
	 y[12] =  C[1]; 	 /* Jovenes = */ 
	 y[13] =  C[2]; 	 /* Adultos = */ 
	 y[14] =  C[3]; 	 /* Ancianos = */ 
	 AUX[15] =  0.0; 	 /* Poblacion = */ 
	 y[16] =  C[4]; 	 /* Minas = */ 
	 AUX[17] =  C[6]; 	 /* Presupuesto = */ 
	 AUX[18] =  0.0; 	 /* deltaingresos = */ 
	 AUX[19] =  0.0; 	 /* deltadeudas = */ 
	 y[20] =  C[5]; 	 /* Ingresos = */ 
	 y[21] =  C[7]; 	 /* Deudas = */ 
	 AUX[22] =  0.0; 	 /* PBIpc = */ 
	 AUX[23] =  0.0; 	 /* DEUpc = */ 
} 


void update_variables_lhs(int c, double C[]) { 
	 /* Updates the following variables at point given by c */ 
	 
	 if (c == 0)   { 
 		 C[8] = 0.025619; 
 		 C[9] = 0.009346; 
 		 C[11] = 0.00235; 
 		 C[14] = 19.2; 
 		 C[15] = 37.2; 
 		 C[16] = 46.4; 
 		 C[17] = 7.27; 
 		 C[18] = 7.27; 
 		 C[19] = 215; 
 		 C[20] = 10650; 
 		 C[21] = 298; 
	 } 
	 if (c == 1)   { 
 		 C[8] = 0.033559; 
 		 C[9] = 0.008246; 
 		 C[11] = 0.00205; 
 		 C[14] = 42.4; 
 		 C[15] = 38.8; 
 		 C[16] = 35.2; 
 		 C[17] = 6.43; 
 		 C[18] = 5.17; 
 		 C[19] = 555; 
 		 C[20] = 10710; 
 		 C[21] = 478; 
	 } 
	 if (c == 2)   { 
 		 C[8] = 0.040705; 
 		 C[9] = 0.008334; 
 		 C[11] = 0.00223; 
 		 C[14] = 49.6; 
 		 C[15] = 40.4; 
 		 C[16] = 67.2; 
 		 C[17] = 8.53; 
 		 C[18] = 7.55; 
 		 C[19] = 415; 
 		 C[20] = 9150; 
 		 C[21] = 538; 
	 } 
	 if (c == 3)   { 
 		 C[8] = 0.016091; 
 		 C[9] = 0.008994; 
 		 C[11] = 0.00237; 
 		 C[14] = 17.6; 
 		 C[15] = 44.4; 
 		 C[16] = 64; 
 		 C[17] = 8.11; 
 		 C[18] = 7.69; 
 		 C[19] = 225; 
 		 C[20] = 11970; 
 		 C[21] = 514; 
	 } 
	 if (c == 4)   { 
 		 C[8] = 0.021649; 
 		 C[9] = 0.008642; 
 		 C[11] = 0.00239; 
 		 C[14] = 44; 
 		 C[15] = 36.4; 
 		 C[16] = 40; 
 		 C[17] = 4.33; 
 		 C[18] = 3.63; 
 		 C[19] = 155; 
 		 C[20] = 11070; 
 		 C[21] = 586; 
	 } 
	 if (c == 5)   { 
 		 C[8] = 0.018473; 
 		 C[9] = 0.00807; 
 		 C[11] = 0.00201; 
 		 C[14] = 10.4; 
 		 C[15] = 54; 
 		 C[16] = 41.6; 
 		 C[17] = 8.67; 
 		 C[18] = 6.85; 
 		 C[19] = 505; 
 		 C[20] = 11850; 
 		 C[21] = 130; 
	 } 
	 if (c == 6)   { 
 		 C[8] = 0.032765; 
 		 C[9] = 0.008158; 
 		 C[11] = 0.00219; 
 		 C[14] = 45.6; 
 		 C[15] = 58.8; 
 		 C[16] = 59.2; 
 		 C[17] = 6.57; 
 		 C[18] = 5.31; 
 		 C[19] = 315; 
 		 C[20] = 10770; 
 		 C[21] = 166; 
	 } 
	 if (c == 7)   { 
 		 C[8] = 0.050233; 
 		 C[9] = 0.008554; 
 		 C[11] = 0.00221; 
 		 C[14] = 39.2; 
 		 C[15] = 52.4; 
 		 C[16] = 51.2; 
 		 C[17] = 3.77; 
 		 C[18] = 8.95; 
 		 C[19] = 365; 
 		 C[20] = 9030; 
 		 C[21] = 106; 
	 } 
	 if (c == 8)   { 
 		 C[8] = 0.039117; 
 		 C[9] = 0.007938; 
 		 C[11] = 0.00265; 
 		 C[14] = 16; 
 		 C[15] = 21.2; 
 		 C[16] = 39.2; 
 		 C[17] = 3.63; 
 		 C[18] = 4.19; 
 		 C[19] = 105; 
 		 C[20] = 10350; 
 		 C[21] = 598; 
	 } 
	 if (c == 9)   { 
 		 C[8] = 0.038323; 
 		 C[9] = 0.007322; 
 		 C[11] = 0.00275; 
 		 C[14] = 40; 
 		 C[15] = 32.4; 
 		 C[16] = 47.2; 
 		 C[17] = 9.51; 
 		 C[18] = 6.71; 
 		 C[19] = 355; 
 		 C[20] = 9690; 
 		 C[21] = 274; 
	 } 
	 if (c == 10)   { 
 		 C[8] = 0.034353; 
 		 C[9] = 0.007762; 
 		 C[11] = 0.00285; 
 		 C[14] = 37.6; 
 		 C[15] = 22.8; 
 		 C[16] = 64.8; 
 		 C[17] = 4.61; 
 		 C[18] = 7.97; 
 		 C[19] = 485; 
 		 C[20] = 10530; 
 		 C[21] = 250; 
	 } 
	 if (c == 11)   { 
 		 C[8] = 0.045469; 
 		 C[9] = 0.00939; 
 		 C[11] = 0.00257; 
 		 C[14] = 11.2; 
 		 C[15] = 45.2; 
 		 C[16] = 56; 
 		 C[17] = 5.17; 
 		 C[18] = 7.83; 
 		 C[19] = 325; 
 		 C[20] = 11910; 
 		 C[21] = 334; 
	 } 
	 if (c == 12)   { 
 		 C[8] = 0.053409; 
 		 C[9] = 0.007894; 
 		 C[11] = 0.00243; 
 		 C[14] = 36.8; 
 		 C[15] = 56.4; 
 		 C[16] = 36.8; 
 		 C[17] = 9.79; 
 		 C[18] = 9.09; 
 		 C[19] = 145; 
 		 C[20] = 10590; 
 		 C[21] = 190; 
	 } 
	 if (c == 13)   { 
 		 C[8] = 0.026413; 
 		 C[9] = 0.007982; 
 		 C[11] = 0.00229; 
 		 C[14] = 18.4; 
 		 C[15] = 27.6; 
 		 C[16] = 49.6; 
 		 C[17] = 6.29; 
 		 C[18] = 8.81; 
 		 C[19] = 375; 
 		 C[20] = 11430; 
 		 C[21] = 382; 
	 } 
	 if (c == 14)   { 
 		 C[8] = 0.017679; 
 		 C[9] = 0.007542; 
 		 C[11] = 0.00259; 
 		 C[14] = 28.8; 
 		 C[15] = 47.6; 
 		 C[16] = 60; 
 		 C[17] = 8.95; 
 		 C[18] = 8.39; 
 		 C[19] = 535; 
 		 C[20] = 9750; 
 		 C[21] = 238; 
	 } 
	 if (c == 15)   { 
 		 C[8] = 0.044675; 
 		 C[9] = 0.008862; 
 		 C[11] = 0.00233; 
 		 C[14] = 41.6; 
 		 C[15] = 49.2; 
 		 C[16] = 42.4; 
 		 C[17] = 9.37; 
 		 C[18] = 9.23; 
 		 C[19] = 425; 
 		 C[20] = 9330; 
 		 C[21] = 430; 
	 } 
	 if (c == 16)   { 
 		 C[8] = 0.023237; 
 		 C[9] = 0.009214; 
 		 C[11] = 0.00263; 
 		 C[14] = 34.4; 
 		 C[15] = 58; 
 		 C[16] = 52.8; 
 		 C[17] = 7.13; 
 		 C[18] = 7.41; 
 		 C[19] = 575; 
 		 C[20] = 10230; 
 		 C[21] = 262; 
	 } 
	 if (c == 17)   { 
 		 C[8] = 0.051027; 
 		 C[9] = 0.00829; 
 		 C[11] = 0.00225; 
 		 C[14] = 25.6; 
 		 C[15] = 22; 
 		 C[16] = 44; 
 		 C[17] = 3.07; 
 		 C[18] = 9.93; 
 		 C[19] = 465; 
 		 C[20] = 10110; 
 		 C[21] = 142; 
	 } 
	 if (c == 18)   { 
 		 C[8] = 0.020855; 
 		 C[9] = 0.007498; 
 		 C[11] = 0.00231; 
 		 C[14] = 29.6; 
 		 C[15] = 39.6; 
 		 C[16] = 45.6; 
 		 C[17] = 9.23; 
 		 C[18] = 6.43; 
 		 C[19] = 305; 
 		 C[20] = 10290; 
 		 C[21] = 358; 
	 } 
	 if (c == 19)   { 
 		 C[8] = 0.043881; 
 		 C[9] = 0.009038; 
 		 C[11] = 0.00253; 
 		 C[14] = 33.6; 
 		 C[15] = 50; 
 		 C[16] = 44.8; 
 		 C[17] = 7.55; 
 		 C[18] = 4.61; 
 		 C[19] = 175; 
 		 C[20] = 11310; 
 		 C[21] = 118; 
	 } 
	 if (c == 20)   { 
 		 C[8] = 0.035147; 
 		 C[9] = 0.009302; 
 		 C[11] = 0.00281; 
 		 C[14] = 14.4; 
 		 C[15] = 20.4; 
 		 C[16] = 56.8; 
 		 C[17] = 4.05; 
 		 C[18] = 5.03; 
 		 C[19] = 545; 
 		 C[20] = 10830; 
 		 C[21] = 694; 
	 } 
	 if (c == 21)   { 
 		 C[8] = 0.037529; 
 		 C[9] = 0.007674; 
 		 C[11] = 0.00273; 
 		 C[14] = 20; 
 		 C[15] = 34; 
 		 C[16] = 48.8; 
 		 C[17] = 4.19; 
 		 C[18] = 6.99; 
 		 C[19] = 265; 
 		 C[20] = 9630; 
 		 C[21] = 670; 
	 } 
	 if (c == 22)   { 
 		 C[8] = 0.024031; 
 		 C[9] = 0.008774; 
 		 C[11] = 0.00217; 
 		 C[14] = 12.8; 
 		 C[15] = 46; 
 		 C[16] = 40.8; 
 		 C[17] = 3.91; 
 		 C[18] = 3.49; 
 		 C[19] = 115; 
 		 C[20] = 9090; 
 		 C[21] = 178; 
	 } 
	 if (c == 23)   { 
 		 C[8] = 0.022443; 
 		 C[9] = 0.008026; 
 		 C[11] = 0.00209; 
 		 C[14] = 26.4; 
 		 C[15] = 59.6; 
 		 C[16] = 66.4; 
 		 C[17] = 3.21; 
 		 C[18] = 3.07; 
 		 C[19] = 245; 
 		 C[20] = 10950; 
 		 C[21] = 682; 
	 } 
	 if (c == 24)   { 
 		 C[8] = 0.048645; 
 		 C[9] = 0.007806; 
 		 C[11] = 0.00291; 
 		 C[14] = 30.4; 
 		 C[15] = 23.6; 
 		 C[16] = 55.2; 
 		 C[17] = 5.73; 
 		 C[18] = 8.25; 
 		 C[19] = 435; 
 		 C[20] = 11010; 
 		 C[21] = 562; 
	 } 
	 if (c == 25)   { 
 		 C[8] = 0.019267; 
 		 C[9] = 0.008422; 
 		 C[11] = 0.00249; 
 		 C[14] = 46.4; 
 		 C[15] = 30.8; 
 		 C[16] = 53.6; 
 		 C[17] = 9.09; 
 		 C[18] = 9.37; 
 		 C[19] = 185; 
 		 C[20] = 9270; 
 		 C[21] = 394; 
	 } 
	 if (c == 26)   { 
 		 C[8] = 0.024825; 
 		 C[9] = 0.007454; 
 		 C[11] = 0.00227; 
 		 C[14] = 44.8; 
 		 C[15] = 57.2; 
 		 C[16] = 37.6; 
 		 C[17] = 7.69; 
 		 C[18] = 7.13; 
 		 C[19] = 125; 
 		 C[20] = 11610; 
 		 C[21] = 574; 
	 } 
	 if (c == 27)   { 
 		 C[8] = 0.016885; 
 		 C[9] = 0.008202; 
 		 C[11] = 0.00279; 
 		 C[14] = 32.8; 
 		 C[15] = 42.8; 
 		 C[16] = 32; 
 		 C[17] = 4.89; 
 		 C[18] = 3.21; 
 		 C[19] = 385; 
 		 C[20] = 10410; 
 		 C[21] = 202; 
	 } 
	 if (c == 28)   { 
 		 C[8] = 0.028001; 
 		 C[9] = 0.009258; 
 		 C[11] = 0.00295; 
 		 C[14] = 43.2; 
 		 C[15] = 38; 
 		 C[16] = 63.2; 
 		 C[17] = 5.87; 
 		 C[18] = 4.33; 
 		 C[19] = 335; 
 		 C[20] = 9450; 
 		 C[21] = 226; 
	 } 
	 if (c == 29)   { 
 		 C[8] = 0.039911; 
 		 C[9] = 0.00741; 
 		 C[11] = 0.00213; 
 		 C[14] = 23.2; 
 		 C[15] = 54.8; 
 		 C[16] = 54.4; 
 		 C[17] = 6.15; 
 		 C[18] = 9.65; 
 		 C[19] = 255; 
 		 C[20] = 11250; 
 		 C[21] = 502; 
	 } 
	 if (c == 30)   { 
 		 C[8] = 0.030383; 
 		 C[9] = 0.008818; 
 		 C[11] = 0.00277; 
 		 C[14] = 38.4; 
 		 C[15] = 50.8; 
 		 C[16] = 58.4; 
 		 C[17] = 5.45; 
 		 C[18] = 4.47; 
 		 C[19] = 195; 
 		 C[20] = 9390; 
 		 C[21] = 442; 
	 } 
	 if (c == 31)   { 
 		 C[8] = 0.036735; 
 		 C[9] = 0.00851; 
 		 C[11] = 0.00271; 
 		 C[14] = 31.2; 
 		 C[15] = 48.4; 
 		 C[16] = 48; 
 		 C[17] = 7.83; 
 		 C[18] = 8.53; 
 		 C[19] = 205; 
 		 C[20] = 10050; 
 		 C[21] = 646; 
	 } 
	 if (c == 32)   { 
 		 C[8] = 0.051821; 
 		 C[9] = 0.00785; 
 		 C[11] = 0.00251; 
 		 C[14] = 24; 
 		 C[15] = 42; 
 		 C[16] = 33.6; 
 		 C[17] = 3.49; 
 		 C[18] = 5.73; 
 		 C[19] = 165; 
 		 C[20] = 9510; 
 		 C[21] = 634; 
	 } 
	 if (c == 33)   { 
 		 C[8] = 0.047057; 
 		 C[9] = 0.007718; 
 		 C[11] = 0.00299; 
 		 C[14] = 21.6; 
 		 C[15] = 24.4; 
 		 C[16] = 65.6; 
 		 C[17] = 7.41; 
 		 C[18] = 6.29; 
 		 C[19] = 295; 
 		 C[20] = 11190; 
 		 C[21] = 346; 
	 } 
	 if (c == 34)   { 
 		 C[8] = 0.020061; 
 		 C[9] = 0.007586; 
 		 C[11] = 0.00207; 
 		 C[14] = 28; 
 		 C[15] = 51.6; 
 		 C[16] = 50.4; 
 		 C[17] = 9.65; 
 		 C[18] = 5.87; 
 		 C[19] = 235; 
 		 C[20] = 10170; 
 		 C[21] = 286; 
	 } 
	 if (c == 35)   { 
 		 C[8] = 0.028795; 
 		 C[9] = 0.007366; 
 		 C[11] = 0.00247; 
 		 C[14] = 22.4; 
 		 C[15] = 35.6; 
 		 C[16] = 34.4; 
 		 C[17] = 5.03; 
 		 C[18] = 6.15; 
 		 C[19] = 475; 
 		 C[20] = 9990; 
 		 C[21] = 466; 
	 } 
	 if (c == 36)   { 
 		 C[8] = 0.042293; 
 		 C[9] = 0.008114; 
 		 C[11] = 0.00269; 
 		 C[14] = 15.2; 
 		 C[15] = 46.8; 
 		 C[16] = 31.2; 
 		 C[17] = 7.97; 
 		 C[18] = 3.91; 
 		 C[19] = 565; 
 		 C[20] = 9930; 
 		 C[21] = 658; 
	 } 
	 if (c == 37)   { 
 		 C[8] = 0.031177; 
 		 C[9] = 0.00895; 
 		 C[11] = 0.00255; 
 		 C[14] = 12; 
 		 C[15] = 33.2; 
 		 C[16] = 30.4; 
 		 C[17] = 4.75; 
 		 C[18] = 4.75; 
 		 C[19] = 515; 
 		 C[20] = 11790; 
 		 C[21] = 310; 
	 } 
	 if (c == 38)   { 
 		 C[8] = 0.035941; 
 		 C[9] = 0.00917; 
 		 C[11] = 0.00245; 
 		 C[14] = 20.8; 
 		 C[15] = 26.8; 
 		 C[16] = 60.8; 
 		 C[17] = 8.81; 
 		 C[18] = 5.45; 
 		 C[19] = 345; 
 		 C[20] = 11130; 
 		 C[21] = 622; 
	 } 
	 if (c == 39)   { 
 		 C[8] = 0.027207; 
 		 C[9] = 0.008906; 
 		 C[11] = 0.00261; 
 		 C[14] = 32; 
 		 C[15] = 34.8; 
 		 C[16] = 52; 
 		 C[17] = 6.85; 
 		 C[18] = 5.59; 
 		 C[19] = 285; 
 		 C[20] = 11550; 
 		 C[21] = 214; 
	 } 
	 if (c == 40)   { 
 		 C[8] = 0.015297; 
 		 C[9] = 0.008686; 
 		 C[11] = 0.00215; 
 		 C[14] = 47.2; 
 		 C[15] = 28.4; 
 		 C[16] = 57.6; 
 		 C[17] = 8.39; 
 		 C[18] = 6.01; 
 		 C[19] = 525; 
 		 C[20] = 9810; 
 		 C[21] = 418; 
	 } 
	 if (c == 41)   { 
 		 C[8] = 0.031971; 
 		 C[9] = 0.009434; 
 		 C[11] = 0.00241; 
 		 C[14] = 35.2; 
 		 C[15] = 25.2; 
 		 C[16] = 61.6; 
 		 C[17] = 8.25; 
 		 C[18] = 6.57; 
 		 C[19] = 585; 
 		 C[20] = 11730; 
 		 C[21] = 550; 
	 } 
	 if (c == 42)   { 
 		 C[8] = 0.049439; 
 		 C[9] = 0.008378; 
 		 C[11] = 0.00297; 
 		 C[14] = 48.8; 
 		 C[15] = 41.2; 
 		 C[16] = 62.4; 
 		 C[17] = 6.01; 
 		 C[18] = 3.77; 
 		 C[19] = 445; 
 		 C[20] = 10470; 
 		 C[21] = 322; 
	 } 
	 if (c == 43)   { 
 		 C[8] = 0.046263; 
 		 C[9] = 0.008598; 
 		 C[11] = 0.00283; 
 		 C[14] = 24.8; 
 		 C[15] = 53.2; 
 		 C[16] = 32.8; 
 		 C[17] = 9.93; 
 		 C[18] = 8.67; 
 		 C[19] = 275; 
 		 C[20] = 11490; 
 		 C[21] = 490; 
	 } 
	 if (c == 44)   { 
 		 C[8] = 0.052615; 
 		 C[9] = 0.009478; 
 		 C[11] = 0.00293; 
 		 C[14] = 48; 
 		 C[15] = 43.6; 
 		 C[16] = 68.8; 
 		 C[17] = 5.31; 
 		 C[18] = 4.05; 
 		 C[19] = 135; 
 		 C[20] = 10890; 
 		 C[21] = 406; 
	 } 
	 if (c == 45)   { 
 		 C[8] = 0.043087; 
 		 C[9] = 0.008466; 
 		 C[11] = 0.00203; 
 		 C[14] = 13.6; 
 		 C[15] = 31.6; 
 		 C[16] = 68; 
 		 C[17] = 3.35; 
 		 C[18] = 3.35; 
 		 C[19] = 455; 
 		 C[20] = 11370; 
 		 C[21] = 454; 
	 } 
	 if (c == 46)   { 
 		 C[8] = 0.054203; 
 		 C[9] = 0.00763; 
 		 C[11] = 0.00289; 
 		 C[14] = 16.8; 
 		 C[15] = 26; 
 		 C[16] = 36; 
 		 C[17] = 6.99; 
 		 C[18] = 8.11; 
 		 C[19] = 395; 
 		 C[20] = 9570; 
 		 C[21] = 370; 
	 } 
	 if (c == 47)   { 
 		 C[8] = 0.029589; 
 		 C[9] = 0.009126; 
 		 C[11] = 0.00211; 
 		 C[14] = 27.2; 
 		 C[15] = 29.2; 
 		 C[16] = 38.4; 
 		 C[17] = 5.59; 
 		 C[18] = 9.51; 
 		 C[19] = 595; 
 		 C[20] = 9870; 
 		 C[21] = 154; 
	 } 
	 if (c == 48)   { 
 		 C[8] = 0.047851; 
 		 C[9] = 0.009082; 
 		 C[11] = 0.00267; 
 		 C[14] = 36; 
 		 C[15] = 30; 
 		 C[16] = 43.2; 
 		 C[17] = 6.71; 
 		 C[18] = 9.79; 
 		 C[19] = 405; 
 		 C[20] = 11670; 
 		 C[21] = 610; 
	 } 
	 if (c == 49)   { 
 		 C[8] = 0.041499; 
 		 C[9] = 0.00873; 
 		 C[11] = 0.00287; 
 		 C[14] = 40.8; 
 		 C[15] = 55.6; 
 		 C[16] = 69.6; 
 		 C[17] = 4.47; 
 		 C[18] = 4.89; 
 		 C[19] = 495; 
 		 C[20] = 9210; 
 		 C[21] = 526; 
	 } 
 } /*  End of update_variables_lhs  */ 


static int nfigsloaded; 
int load_data_figures(double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, int jj, double y[], double C[], double AUX[], double x, FILE * myfile) { 
	 int j = nfigsloaded; 
	 int i; int errorover=0; 
		  FIGS[0][j] = x; 
		  if (! finite(AUX[15])) { printf("Overflow error on variable Poblacion \n"); errorover++; } 
		  FIGS[1][j] = AUX[15];
		  if (figdata[1].max < FIGS[1][j]) { figdata[1].max =  FIGS[1][j]; } 
		  if (figdata[1].min > FIGS[1][j]) { figdata[1].min =  FIGS[1][j]; } 
		  if (! finite(y[16])) { printf("Overflow error on variable Minas \n"); errorover++; } 
		  FIGS[2][j] = y[16];
		  if (figdata[2].max < FIGS[2][j]) { figdata[2].max =  FIGS[2][j]; } 
		  if (figdata[2].min > FIGS[2][j]) { figdata[2].min =  FIGS[2][j]; } 
		  if (! finite(y[20])) { printf("Overflow error on variable Ingresos \n"); errorover++; } 
		  FIGS[3][j] = y[20];
		  if (figdata[3].max < FIGS[3][j]) { figdata[3].max =  FIGS[3][j]; } 
		  if (figdata[3].min > FIGS[3][j]) { figdata[3].min =  FIGS[3][j]; } 
		  if (! finite(y[21])) { printf("Overflow error on variable Deudas \n"); errorover++; } 
		  FIGS[4][j] = y[21];
		  if (figdata[4].max < FIGS[4][j]) { figdata[4].max =  FIGS[4][j]; } 
		  if (figdata[4].min > FIGS[4][j]) { figdata[4].min =  FIGS[4][j]; } 
		  if (! finite(AUX[17])) { printf("Overflow error on variable Presupuesto \n"); errorover++; } 
		  FIGS[5][j] = AUX[17];
		  if (figdata[5].max < FIGS[5][j]) { figdata[5].max =  FIGS[5][j]; } 
		  if (figdata[5].min > FIGS[5][j]) { figdata[5].min =  FIGS[5][j]; } 
		  if (! finite(AUX[22])) { printf("Overflow error on variable PBIpc \n"); errorover++; } 
		  FIGS[6][j] = AUX[22];
		  if (figdata[6].max < FIGS[6][j]) { figdata[6].max =  FIGS[6][j]; } 
		  if (figdata[6].min > FIGS[6][j]) { figdata[6].min =  FIGS[6][j]; } 
		  if (! finite(AUX[23])) { printf("Overflow error on variable DEUpc \n"); errorover++; } 
		  FIGS[7][j] = AUX[23];
		  if (figdata[7].max < FIGS[7][j]) { figdata[7].max =  FIGS[7][j]; } 
		  if (figdata[7].min > FIGS[7][j]) { figdata[7].min =  FIGS[7][j]; } 
	 nfigsloaded ++;
	 if (nfigsloaded >= MAX_FIGS ) { int xx, yy; 
		 for (xx=0;xx<MAX_FIGS;xx++) { 
		 fprintf(myfile,"%f\t",FIGS[0][xx]); 
		  for (yy=1;yy<=7;yy++) 
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
	 int nfigs = 7;
	 eqfig figdata[9];
	 double FIGS[9][MAX_FIGS]; 
	 int maxloop =0, nloop =0; int figres=0; 
	 y    =  dvector(1,N);  
	 dydx =  dvector(1,N); 
	 yout =  dvector(1,N);  
	 C   =   dvector(1,NINC); 
	 AUX =   dvector(1,N); 
	 nfigsloaded = 0; 
	 for(i=1;i<=nfigs;i++) { figdata[i].max = -HUGE_VAL; figdata[i].min=HUGE_VAL; }; 
	  strcpy(figdata[1].name,"Poblacion");
	  strcpy(figdata[2].name,"Minas");
	  strcpy(figdata[3].name,"Ingresos");
	  strcpy(figdata[4].name,"Deudas");
	  strcpy(figdata[5].name,"Presupuesto");
	  strcpy(figdata[6].name,"PBIpc");
	  strcpy(figdata[7].name,"DEUpc");
	 while (myc < 50) { 
	 load_constants(C); 
	 update_variables_lhs(myc, C); 
	 load_auxlevel(AUX,y,C); 
	 nfigsloaded = 0; 
	 /* ---rk4--- */ 
	 h = (double)1;
	 x = (double)0; 
	 maxloop = (int) (200/h)+1;
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
	 WritePlotScript(1, nfigs, 240, figdata, 0, myc-1);
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
