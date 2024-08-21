/* Based on the lhs code - rosworks.c */

typedef struct 
{ 
   long double outcome[NUMPRCCOUTS]; 
   long double time; 
} timepoint; 
 
typedef struct       /* Structure that holds info for one equation */ 
{ 
   char name[MAXNM]; 
   double top; 
}  equationinfo;  
 
typedef struct 
{ 
     char* name[20]; 
     double min; 
     double peak; 
     double max; 
     char dType; 
} paramT; 
  
typedef struct 
{ 
     long double value; 
     double rank; 
     int tieset; 
} matrixDB; 
 
void  FillHypercube(void); 
void  FillUniformHypercube(double min, double max,int element); 
void  FillTriangularHypercube(double min, double peak, double max,int element); 
void  FillIntervalHypercube(double min, double max,int element);  
void  FillLogHypercube(double min, double max,int element);  
void  RankTieSort(void); 
void  PrintParamSamples(void); 
void  FillInitialConditions(int sim); 
void  UncertaintyAnalysis( void); 
void  CheckInverse( void); 
void  GetPRCC(void ); 
void  FillCijMatrix(void); 
void  CopyVaryPar(void); 

void  AdamsBash(long double *DidIncrement, long double *NextIncrement,long double *WorkValArray);  
void  RungeKutta(long double *DidIncrement,long double *NextIncrement,long double *WorkValArray, long double *derivative); 

void  MakeLog(void); 
void  MakeGraphData(void); 
void  WriteTopValues(void); 
void  Error(char* msg); 
void  ErrorSim2(char* msg); 
int   RandomInteger(int low, int high);    
void  StringClean(int typeclean); 
void  CalculatedOut(int numkept,long double t); 

void  Rosenbrock(long double *ScaleArray, long double *DidIncrement,long double *NextIncrement,long double *WorkValArray, 
      long double* derivArray,long double *timederiv,long double eps); 
void  RKAdapt( long double *ScaleArray, long double *DidIncrement, long double *NextIncrement, long double *WorkValArray, 
      long double *derivArray,long double eps); 


void  FillJacobian(long double t); 
void  BSStep(long double *ScaleArray, long double *DidIncrement,long 
double *NextIncrement, long double *WorkValArray, long 
      double *derivArray, long double *timederiv,long double eps); 
void  ODEStepper(long double TimeEnd); 

void  PrintMatrix(FILE* output,long double* matrix,int n,int m); 
void  LUDecomp(long double *matrix00, int n, int *indx, long double *d); 
void  LUBackSub(long double* matrix00, int n, int *indx, long double *b); 
void  InvertMatrix(long double *matrix,long double *inverse, int dim); 
void  MatMult(long double *matrix1,long double *matrix2,long double *product,int n,int m,int p); 
void  MatCopy(long double *target, long double *source,int n, int m); 
static void *GetBlock(size_t nbytes); 
void  OutputRunParams(int sim); 
void  GetPreviousParams(int sim); 
void  FiniteDiffJacob (int n, long double *indvar,long double *fvec, long double *df00); 
void  ModifiedMid( long double *yin, int numeq,long double start,long double bigstepsize, 
						 int numsubsteps,long double *yout,long double *deriv); 
void  PolyExtrapolate( int numpts,long double xest,long double *yest,long double *yextrap,long double *ydiff,int numeq); 
void  SemiImpModMid ( long double *yin,long double *timederiv,int numeq, long double start, long double bigstepsize,  
      int numsubsteps, long double *yout,long double *deriv); 
void  GetTimeDeriv(int n, long double *indvar,long double *fvec, long double *df00); 
void  PrintLabelSample(int element); 
void SaveOutComeData(int ptssaved,int sim); 
void PlaceNextOutCome(int currentstep); 
void CheckStrongEffect(void);         /* Subroutine that checks to see if a strong effect (PRCC==1 or PRCC==-1) parameter exists */ 
  
/* input output files */ 
FILE  *infile[2];               /* Input file with the parameters and ranges specified*/ 
FILE  *check[6];             /* Output files for results */ 
FILE  *warn1;                /* Dump error messages, progress messages   FILE= Warning_1 */ 
FILE  *prcc[NUMPRCCOUTS];   /* Files for printing results of time-dependent PRCCs */ 

/*Variable declarations*/ 
long double   TimeMatrix[STEPSTOKEEP][NUMEQ+1];      /*Array for  writing solutions ever SaveFreq time units, equations row1=initial conditions*/ 
matrixDB      LHSMatrix[MAXPARS][NUMSIM];        /* 3-D Matrix -parameter and outcome values & ranks for each run*/ 
long double   RankMatrix[MAXVARPARS][NUMSIM]; /* Matrix of ranks for parameters that vary  */ 
long double   CijMatrix[MAXVARPARS][MAXVARPARS];  /* Calculated from RankMatrix, Cij^1 used to calculate PRCCs */ 
long double   CijInverse[MAXVARPARS][MAXVARPARS];  /* Used in inversion of Cij in the calculation of the PRCCs */ 
paramT        pArray[MAXPARS];                    /* Array of parameter distribution information K parameters */ 
long double   CalcTime[STEPSTOKEEP][NUMCALCOUT] ; /* Array of calculated outcomes */ 
long double   Jacobian[NUMEQ][NUMEQ];       /* The Jacobian matrix used in the Rosenbrock stiff solver */ 
timepoint     OutComeData[NUMSIM][POINTSTOKEEP];    



int          RKCalc; 
long double  TimeAct;          /* The current timepoint for the integration*/ 
long double  TimeIncrement;    /* The current time increment*/ 
int OutNum;                    /* The number of the outcome of interest */  
equationinfo  Eqns[NUMEQ]; 
equationinfo  CalcOut[NUMCALCOUT]; 

long double OutDay;            /* The outcome time chosen */ 
char        OutName[MAXLINE];  /* The name of the outcome variable used for PRCC*/ 
long double CurrOutTime;       /* The time for the current PRCC outcome variable - should be close to OutDay for Local */ 
int         LastStep; 
int         MinSaveNo=STEPSTOKEEP;    

long double PRCCOffTime =ENDDAY;/* Indicator for shutting off doing the PRCCs -if an outcometime is not reached for one of thesimulations */ 
int         DoLocal=0;          /* Indicator for calculating the PRCC - if the outcome time falls in the current interval ofintegration */     
int         DoGlobal=0;         /* Indicator for doing the globAL PRCC */ 
int         DonePRCC=0;         /* Indicator to stop running PRCC analyses, all requested points have been done */ 
int         Method;             /* method used to solve differential equation */      
long double StartTime;          /* Beginning point for the integration */ 
long double EndDay;             /* Ending point for the integration */ 
int         NumParDer;          /* The number of nonindependent parameters (derived from other parameters ) */ 
int         NumPar;             /* The number of parameters */ 
int         NumVarPar ;         /* The number of parameters,the number of parameters that vary*/  
int         NumSim;             /* The number of simulations*/ 
int         sim;                /* counter for current simulation number */      
int         Call=1;             /* How many PRCCs have been calculated on the current run*/ 
int         FirstCall=1; 
int         PrintPRCCResults;   /* Do we print results to PRCC output files now 1=Yes, 0=No*/ 
,MAXPARS*NUMSIM*si
time_t        now;                                /* Used to track the time */ 
time_t        start_time,current_time; 
double        run_duration; 
struct tm     *ptr; 
char          *c; 

int           NewI[NUMSIM];      
int            Invert=1;                          /* Used for log files to flag whether inversion of matrix was successful  */ 
int            ErrFlag,FirstErr=0;      
char          *parfile,*eqnfile; 
int           Warn=0; 
char          junk[80];     
char          line[MAXLINE];    
char          errmessage[MAXLINE]; 
int           NumEq = NUMEQ; 
int           NumCalcOut ; 
int           LabelPar=LABELNUM; 
long double   SaveFreq=SAVEFREQ; 
int           currentoutcome=0;     
int           SimDone=0; 
int           StrongEffect,StrongPar,Parity=1; 

 
/* Variables used in new polynomial extrapolation, modified midpoint, Burlisch-Stoer,            */ 
/* semiimplicit modified midpoint, , stiff Burlisch Stoer, and finite differencing subroutines   */ 
/* that need to be defined more clearly                                                          */ 
 
long double  d[NUMEQ][KMAXX]; 
long double  x[KMAXX]; 
 
/* array to keep track of max values of equation variables, used to 
   defined yranges of graphs */ 
long double TopValues[NUMEQ]; 
 
/*SYSTEM VARS HERE*/ 

#include "C_VarDeclars"
 
  
main(int argc, char **argv) 
{ 
     int i, j,element; 
     time(&now); 
     ptr=localtime(&now); 
     c=asctime(ptr); 
     puts(c); 
     printf("Running LHS code version 1.0: Last update 11/11/00\n\n\n\n"); 
 
      
/*Read in parameters,select outcome variable and number of runs */ 
     if(argc<2) 
     { 
        printf("Please enter the name of parameter file to be used\n"); 
        exit(1); 
     } 
     infile[0] = fopen(argv[1], "r"); 
     parfile=argv[1]; 
     if(infile[0] == NULL) Error("Cannot find parameter file2 "); 
 
     if(argc<3) 
     { 
        printf("Please enter the name of the equation file to be used\n"); 
        exit(1); 
     } 
     infile[1] = fopen(argv[2], "r"); 
     eqnfile = argv[2]; 
     if(infile[1] == NULL) Error("Cannot find equation file2 "); 
     
     if (argc <4)  
     { 
        NumSim=NUMSIM+1; 
        do    
        { 
           printf( "How many runs would you like to  do? Current max is %d \n",NUMSIM); 
           scanf("%d[^\n]", &NumSim); 
           gets(junk); 
        } 
       while (NumSim>NUMSIM); 
     } 
     else  NumSim=atoi(argv[3]); 
      
     if (argc<5) /* NEEDS FIXING FOR BATCH MODE */ 
     {    
        OutNum=999; 
        do   
        { 
            printf("\n\n\nWhich outcome variable would you like?\n" ); 
            system("more ./.outcomenames"); 
            scanf("%d[^\n]", &OutNum); 
            gets(junk); 
        } 
        while ( (OutNum>=NUMEQ+1+NUMCALCOUT) || (OutNum<0));           
     } 
     else OutNum=atoi(argv[4]); 
           
     if (argc<6)     
     { 
        OutDay=ENDDAY+1.0; 
        do  
        {      
           printf( "What time would you like to see? Current max time is %f \n",(double)(ENDDAY)); 
           scanf("%lf[^\n]", &OutDay); 
           gets(junk); 
        } 
        while (OutDay>ENDDAY); 
     } 
     else OutDay=atof(argv[5]); 
     printf( "\n\n\n\nThe outcome day was %lf\n\n\n\n\n",OutDay); 
      
           
     if ( argc<7)  
     { 
        Method=0; 
        do  
        { 
           printf( "What method would you like to use to solve the equations ? \n\n"); 
           printf("\t1\t=\tFourth order Runge-Kutta\n"); 
           printf("\t2\t=\tAdams-Bashforth-Moulton 4th order Predictor Corrector\n"); 
           printf("\t\t=\tRosenbrock (stiff) Solver\n"); 
           printf("\t4\t=\tRunge-Kutta with Adaptive Stepsize\n"); 
           printf("\t5\t=\tBurlisch-Stoer with Adaptive Stepsize\n"); 
           printf("\t6\t=\tStiff Burlisch-Stoer with Adaptive Stepsize\n"); 
           printf("\t7\t=\tRosenbrock - Jacobian precalculated \n"); 
           scanf("%d[^\n]", &Method); 
           gets(junk); 
        } 
        while((Method<1) || (Method>7)); 
     } 
     else Method=atoi(argv[6]); 
           
     if (argc<8) StartTime = 0.0; 
     else StartTime= atof(argv[7]); 
     if (argc<9) EndDay=ENDDAY; 
     else EndDay=atof(argv[8]); 
       
 
     printf("Start time: %lf \t End Time: %lf ENDDAY: %lf\n",StartTime,EndDay,ENDDAY); 
     if(((OutDay+TINY>StartTime+TINY)&& (OutDay+TINY<=EndDay+TINY)) &&(PRCCMODE==0))   /* If we're doing local PRCCs, and the*/ 
     {                                                                                 /* PRCC outcome time falls into the interval*/  
         DoLocal=1;                                                                    /* of integration (might not if we're doing */ 
     } 
     if ((PRCCMODE==2) && (NumSim>1)) DoGlobal=1; 
/*Set up output files */ 
      
     if ((DoLocal==1) || (DoGlobal==1))check[1]=fopen("Check_Ranks","w"); 
     if (StartTime<TINY) check[2]=fopen("Check_Param", "w"); 
     if((DoLocal==1) || (DoGlobal==1))check[3]=fopen("Check_Inverse","w"); 
     if( (DoLocal==1) || (DoGlobal==1)) check[4]=fopen("Check_Uncert","a"); 
     if ((DoLocal==1)|| (DoGlobal==1)) check[5]=fopen("Check_PRCC","a"); 
     check[0]=fopen("Check_","w"); 
     warn1=fopen("Log_File","a"); 
      
      
     if(StartTime<TINY) fprintf(check[2], "\nCheck_Param: Parameter Selections for each of the N runs:\nDate: %s ", c); 
     if ((DoLocal==1) || (DoGlobal==1)) 
     { 
        fprintf(check[1], "\nCheckRanks: Values and ranks for parameters that vary \nDate: %s ",c); 
        fprintf(check[3], "Check_Inverse : Check the inversion of Cij matrix\nDate: %s ", c); 
        fprintf(check[4], "\nCheck_Uncert: Results for the uncertainty analysis: \nDate: %s ", c); 
        fprintf(check[5], "\nCheck_PRCC: Results for the PRCC analysis: \nDate: %s ", c); 
        fprintf(check[1],"_________________________________________________________________________\n\n"); 
        for(i=3;i<6;i++)fprintf(check[i],"_________________________________________________________________________\n"); 
        fprintf(warn1,"_________________________________________________________________________\n\n"); 
        fprintf(warn1, "Error messages for run on : %s ", c); 
        for (i=0;i<6;i++) fflush(check[i]); 
     }  
 
/* Initialize matrices here */ 
 
     memset(OutComeData,0,NUMSIM*POINTSTOKEEP*sizeof(timepoint)); 
     memset(LHSMatrix,0,MAXPARS*NUMSIM*sizeof(matrixDB)); 
     memset(RankMatrix,0,MAXVARPARS*NUMSIM*sizeof(long double)); 
     memset(CijMatrix,0,MAXVARPARS*MAXVARPARS*sizeof(long double)); 
     memset(CijInverse,0,MAXVARPARS*MAXVARPARS*sizeof(long double)); 
     memset(Jacobian,0,(NUMEQ)*(NUMEQ)*sizeof(long double)); 
     memset(NewI,0,NUMSIM*sizeof(int)); 
     memset(d,0,(NUMEQ)*KMAXX*sizeof(long double));	 
     memset(x,0,(KMAXX)*sizeof(long double)); 
 
/* Begin primary analysis   */ 
 
     ReadParameters();  
     if (NumVarPar<1)                     /* if there are no varying parameters, shut off the PRCC analysis */ 
     { 
        DoGlobal=0; 
        DoLocal=0; 
        DonePRCC=1; 
     } 
     ReadEquations();     
     if (StartTime<=TINY)  
     { 
        FillHypercube(); 
        printf("\n\nNumber of independent parameters read in :%d\n",NumPar); 
        PrintParamSamples(); 
     } 
     else for(i=0;i<NumSim;i++) GetPreviousParams(i); 
     for (i=0;i<NumEq;i++) 
       TopValues[i] = 0.0; 
 
     for(sim=0; sim < NumSim; sim++) 
     {     
        SimDone=0;                                                  /* Set the flag for terminating the simulation off */ 
        start_time=time(0);                                         /* Get the starting time for the simulation */ 
        memset(TimeMatrix,0,(NUMEQ+1)*STEPSTOKEEP*sizeof(long double));/* Initialize the time matrix  */ 
        FillInitialConditions(sim);                                /* Fill in the initial conditions */ 
        ODEStepper(EndDay);                                        /* Integrate the equations */ 
        MakeGraphData();                                           /* Write out the data for the graph files */ 
        if(EndDay<ENDDAY) OutputRunParams(sim);                    /* If we're doing piecewise continuous integration, and were not*/ 
                                                                   /* to the endpoint yet, save the current parameter samplings    */ 
                                                                   /* and endpoint values of each of the equations to an external file */ 
     } 
       
     WriteTopValues(); 
      
     if (DoGlobal==1) CurrOutTime=0.0; 
     if ( ((DoLocal==1) || (DoGlobal==1)) && (NumSim>1))   /* if we want to do prcc analysis and all outcome times are reached */ 
     { 
       do 
        { 
              if ((DoLocal==1) || ((DoGlobal==1) && (Call==1)))PrintPRCCResults=1; 
              PlaceNextOutCome(Call);  
              Invert=1; 
              Parity=1; 
              StrongEffect=0; 
              StrongPar=NumVarPar; 
              UncertaintyAnalysis(); 
              RankTieSort(); 
              CopyVaryPar();  
              FillCijMatrix(); 
              InvertMatrix((long double *) CijMatrix[0],(long double *)CijInverse[0],NumVarPar+1); 
              if (Invert==0) CheckStrongEffect(); 
              if (Invert==1) CheckInverse(); 
              GetPRCC(); 
              if (DoGlobal==1)  
              { 
                 if (CurrOutTime>=(EndDay-SaveFreq)) DonePRCC=1; 
                 if ((DonePRCC==1) && (CurrOutTime<EndDay)) 
                 {      
                    sprintf(errmessage, "\n\n\nNOTE: Some simulations were terminated before integration endpoint, calculations terminate early"); 
                    fprintf(check[1],errmessage); 
                    fprintf(check[3],errmessage); 
                    fprintf(check[4],errmessage); 
                    fprintf(check[5],errmessage); 
                    fprintf(warn1, errmessage); 
                 } 
              } 
              if (DoLocal==1) DonePRCC=1;   
              Call++;  
        } 
        while ((DonePRCC==0)); 
     } 
      
     sprintf(errmessage,"Run completed successfully"); 
     if ((DoLocal==1) || (DoGlobal==1)) MakeLog(); 
     return 0; 
} 
 
 
/***************************************************************************************************************/ 
/* This subroutine reads in the values of the parameters into the array pArray from the file                   */ 
/*													       */	 
/* Values read in are:                                                                                         */ 
/*                                                                                                             */ 
/*         pArray[i].name ( the name of the ith parameter)                                                     */ 
/*         pArray[i].min  ( the lower end of the distribution for sampling the ith parameter                   */ 
/*         pArray[i].peak ( the value in the range of the distribution where the peak of the triangle occurs)  */ 
/*         pArray[i].max  ( the upper end of the distribution for sampling the ith parameter                   */ 
/*         pArray[i].dType ( the distribution type for the ith parameter (T=triangular, U=Uniform )            */ 
/*         NumPar          (the number of parameters read in)                                                  */ 
/***************************************************************************************************************/ 
 
void ReadParameters(void) 
{ 
 
   int lastchar; 
   NumPar=0; 
   NumVarPar=0; 
 
   while (fgets(line,MAXLINE,infile[0])!=NULL)  
   { 
      StringClean(0); 
      if(strlen(line)==0) break; 
      if (strncmp(line,"Initial",7)==0) continue; 
      if(strncmp(line,"Parameters",10)==0)continue; 
      else  
      { 
         lastchar=(strlen(line)-2);  
         if(line[lastchar]=='D') 
	   pArray[NumPar].dType = 'D'; 
         else  
         { 
            if ((line[lastchar]=='U')||(line[lastchar]=='T') ||(line[lastchar]=='C') || (line[lastchar]=='L')) 
              sscanf(line, "%s %lf %lf %lf %c\n",pArray[NumPar].name,&pArray[NumPar].min,&pArray[NumPar].peak,&pArray[NumPar].max,&pArray[NumPar].dType); 
	      if (pArray[NumPar].max<pArray[NumPar].min) 
              { 
                 sprintf(errmessage, "Minimum for %s exceeds the maximum in parameter file", pArray[NumPar].name); 
                 Error(errmessage); 
              } 
              if (pArray[NumPar].max!=pArray[NumPar].min) NumVarPar++; 
              NumPar++; 
         }  
      }  
   }  
   printf("Number of independent parameters: %d\n",NumPar); 
   fclose(infile[0]); 
}   
 
 
 
/* The next three subroutines take the upper and lower limits of the          */ 
/* distributions, divide them into equiprobable regions, sample the midpoints */ 
/* of the regions,  randomize them and put them in the LHC matrix             */ 
 
void FillHypercube(void) 
{ 
     int      i,element,loc; 
     double   temp; 
      
    srand((int) time(NULL)); 
    for(element = 0; element< NumPar; element++) 
    { 
        printf("Reading in Parameter %d from input file...\n", element+1); 
        switch(pArray[element].dType)     /* which type of row to fill in */ 
        { 
            case 'U': FillUniformHypercube( pArray[element].min,pArray[element].max,element); break; 
            case 'C': FillIntervalHypercube(pArray[element].min,pArray[element].max,element);break; 
            case 'T': FillTriangularHypercube(pArray[element].min,pArray[element].peak,pArray[element].max,element);break; 
            case 'L': FillLogHypercube(pArray[element].min,pArray[element].max,element);break; 
            case 'D': break; 
            default: 
            { 
	        fprintf(warn1,"%c",pArray[element].dType);  
                fflush(warn1); 
		Error("Parameters distn type must be T or U.");  
		break; 
            } 
	 } 
         for(i=0; i<NumSim; i++) 
         { 
            loc = RandomInteger(i, NumSim-1); 
            SWAP(LHSMatrix[element][i].value, LHSMatrix[element][loc].value); 
         } 
    if (element==LabelPar) PrintLabelSample(element);       
    } 
} 
 
void FillTriangularHypercube(double min, double peak, double max,int element) 
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
    LHSMatrix[element][i].value = (ximin + ximax)/2; 
    ximin = ximax; 
  } 
} 
 
 
 
 /* This subroutine fills in one row of the LHSMatrix with sampled values       */ 
 /* for the current parameter, which has a uniform distribution. Itakes the     */ 
 /* lower and upper limits of the distribution and determines what the          */ 
 /* midpoint of each of region would be if the distribution were divided into   */ 
 /* NumSim equiprobable regions. These midpoints are the sampled values that    */ 
 /* are initially put into the LHS matrix in order from smallest to largest     */ 
 
 void FillUniformHypercube( double distmin, double distmax,int element) 
 { 
      double incr,min; 
      int    k; 
 
      incr = (distmax-distmin)/NumSim; 
      min=distmin; 
      min += incr/2; 
      for(k = 0; k< NumSim; k++) 
      { 
	 LHSMatrix[element][k].value = min; 
	 min += incr; 
      } 
 } 
 
 
 /* This subroutine sampes evenly across the entire interval instead of just */ 
 /* at the midpoints of equiprobable regions                                 */ 
 
 void FillIntervalHypercube( double distmin, double distmax,int element) 
 { 
      double incr,min; 
      int    k; 
 
      incr = (distmax-distmin)/(NumSim-1);   /* width of each subinterval is max-min / number of simulations */       
      min=distmin; 
      for(k = 0; k< NumSim; k++)   
      { 
	 LHSMatrix[element][k].value = min;   /* store the last midpoint value in the LHSMatrix */      
	 min += incr;			     /* next midpoint = current midpoint + size of the interval */	 
      } 
 } 
 
 
 void FillLogHypercube(double distmin, double distmax,int element) 
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
	 LHSMatrix[element][k].value = exp(min*log(10));    
	 min += increment;   
     } 
 } 
 
 
 
 
 /* Equation Solvers: The next eight subroutines fill the TimeArray which is the matrix of all the     */ 
 /* solutions ranging over time - first row contains the initial conditions, second contains          */ 
 /* values of equations after 1 time increment obtained by solving (e.g. using 4th order Runge-Kutta) */  
 
 /* FillInitialConditions fills in the initial conditions, and when a stiff solver is being used, */ 
 /* the constant values of the Jacobian that don't need to be solved. */ 
 /* It then calls the driver subroutine (ODEStepper) which starts and stops the integration*/ 
 
 void FillInitialConditions(int sim)          
 { 
    long double derivArray[NUMEQ] ; 
    long double TimeDerivArray[NUMEQ]; 
    int i; 
 
    printf("\nRunning simulation # %d....\n", sim + 1);    
 
 /* ASSIGN EQUATION-SPECIFIC VARIABLES HERE */ 

#include "C_VarAssigns"
 
 /*PUT INITIAL CONDITIONS AND CONSTANT JACOBIAN TERMS HERE */ 
#include "C_InitialCond"

 
    TimeMatrix[0][NUMEQ]=StartTime; 
    FirstErr=0; 
    TimeAct=StartTime; 
    for (i=0;i<NUMEQ;i++) TimeDerivArray[i] = 0.0;   /* This will need to be changed to accomodate time dependencies */ 
    GetDerivatives((long double *) TimeMatrix[0], (long double *) derivArray, TimeAct);     
    if ((Method==3) || ((Method==6) || (Method==7)))GetTimeDeriv(1,(long double *) TimeMatrix[0],  
    (long double *) derivArray, (long double *) TimeDerivArray);  
    TimeIncrement =(long double) INCREMENT; 
    if (Method==2) LastStep=STEPSTOKEEP+5; 
    else LastStep=STEPSTOKEEP+1; 
    if ((Method==3) || (Method==6))  FiniteDiffJacob(NumEq,(long double *) TimeMatrix[0],(long double *) derivArray,  
     (long double *)Jacobian[0]); 
    if (Method==7) FillJacobian(TimeAct); 
 } 
 
 
 /* This subroutine performs one step of a fourth order Runge-Kutta solver,with fixed stepsized. It is called by ODEStepper */ 
 /* Comment 1: Compute K1/h = f(t,w) */ 
 /* Comment 2: Compute K2/h = f(t+h/2,w+K1/2) store this in dyt */ 
 /* Comment 3: Compute K3/h = f(t+h/2,w+K2/2) - store this in dym*/ 
 /* Comment 4: Comput w+K3 , store this in yt */ 
 /* Comment 5: Compute K2+K3/h - Store this in dym */ 
 /* Comment 6: Compute K4/h = f(t+h,w+K3) - Store this in dyt */ 
 /* Comment 7: Compute the wnew = wold + (k1+2(k2+k3)+k4)/6  */ 
 /*              If w<0, set w = 0                             */ 
 
 void RungeKutta(long double *DidIncrement, long double *NextIncrement,long double *WorkValArray, long double *derivative) 
 { 
    long double  dym[NUMEQ], dyt[NUMEQ] ,yt[NUMEQ]; 
    int          k; 
    long double  EstimateArray[NUMEQ]; 
 
    memset(dym,0,(NUMEQ)*sizeof(long double)); 
    memset(dyt,0,(NUMEQ)*sizeof(long double)); 
    memset(yt,0,(NUMEQ)*sizeof(long double)); 
 
    for(k=0; k<NUMEQ;k++)yt[k]=(*(WorkValArray+k))+(TimeIncrement*0.500000000000)*(*(derivative+k)); 
 /* 1*/ 
 
    GetDerivatives((long double *)yt, (long double*)dyt,TimeAct+TimeIncrement*0.5);  /* 2 */ 
 
    for(k=0;k<NUMEQ;k++)yt[k]=(*(WorkValArray+k))+TimeIncrement*0.5*dyt[k]; /* 3 */ 
    GetDerivatives((long double *)yt, (long double *)dym,TimeAct+TimeIncrement*0.5); 
 
    for(k=0; k<NUMEQ; k++) 
    {         
       yt[k]= (*(WorkValArray+k))+TimeIncrement*dym[k];          /* 4 */ 
       dym[k] += dyt[k];                                         /* 5 */ 
    } 
 
    GetDerivatives((long double *)yt, (long double *)dyt,TimeAct+TimeIncrement); /* 6*/ 
    for(k=0;k<NUMEQ;k++)     /* 7 */ 
    { 
       (*(WorkValArray+k))+=(TimeIncrement/6.0)*((*(derivative+k))+dyt[k]+2.0*dym[k]); 
       if (*(WorkValArray+k)<-TINY) 
       {       
	  (*(WorkValArray+k))=0.0; 
	  if (Warn==0) fprintf(warn1, "\tWarning: some equation values were set to zero \n \twhen solver produced negative outcome   values at time %lf\n", TimeAct); 
	  Warn=1; 
       } 
    } 
    GetDerivatives((long double *) WorkValArray,(long double *) derivative,TimeAct); 
    TimeAct+=TimeIncrement; 
    *DidIncrement= TimeIncrement;  
    *NextIncrement= TimeIncrement; 
    RKCalc=1; 
 } 
 
 /* This subroutine does one iteration of an Adam-Bashforth-Moulton fourth-order                                 */ 
 /* predictor corrector. Note that  in order to optimize code calculation, the predictor has been rewritten      */ 
 /* in terms of previous calculations after the first call. Also, note that after the predictor is calculated,   */ 
 /* derivatives are moved back by one before calculating the corrector - to save recalculation. derivs[0][j] is  */ 
 /* the value of jth derivative at iteration current-4                                                           */ 
 
 
 
 void AdamsBash(long double *DidIncrement, long double *NextIncrement,long double *WorkValArray) 
 { 
    static long double    derivs[4][NUMEQ];                 
    static long double    LastTimes[4][NUMEQ];   
    int                   k,l,j; 
    long double           t1,t8; 
    static long double    t2[NUMEQ],t2prev[NUMEQ]; 
    long double           pred[NUMEQ]; 
    long double           EstimateArray[NUMEQ]; 
 
 
     t1 = TimeIncrement/24.0; 
     t8 = TimeIncrement/3.0; 
     if (RKCalc==1) for (j=0;j<4;j++) for(k=0;k<NUMEQ;k++)LastTimes[j][k]=TimeMatrix[j][k]; 
     for (j=0;j<4;j++) GetDerivatives((long double *) LastTimes[j],(long double *) derivs[j],TimeAct); 
     if (RKCalc == 1) 
     { 
	for(k=0; k<NUMEQ; k++) 
	{ 
	   t2[k]=5*derivs[2][k]; 
	   pred[k]=LastTimes[3][k]+t1*(55*derivs[3][k]-59*derivs[2][k]+37*derivs[1][k] - 9*derivs[0][k]); 
	} 
    } 
    else    
    { 
       for(k=0; k<NUMEQ; k++) 
       { 
	  t2prev[k] = t2[k]; 
	  t2[k]=5*derivs[2][k]; 
	  pred[k]= LastTimes[2][k]+t8*(8*derivs[3][k]-t2[k] +t2prev[k]-derivs[1][k] - derivs[0][k]); 
       } 
    } 
    for (k=0; k<NUMEQ;k++) for (l=1;l<4;l++)  
    { 
       derivs[l-1][k] =derivs[l][k];    
       LastTimes[l-1][k]= LastTimes[l][k]; 
    } 
 
    GetDerivatives ( (long double *) pred, (long double *) derivs[3],TimeAct);   
    for (k=0; k< NUMEQ; k++)  
    { 
       EstimateArray[k] =LastTimes[3][k]+t1*(9*derivs[3][k]+19*derivs[2][k]- t2[k] + derivs[0][k]); 
       if (EstimateArray[k]<-TINY)  
       { 
	  EstimateArray[k]=0.0; 
	  if (Warn==0) fprintf(warn1, "\tWarning: some equation values were set to zero \n \twhen solver produced negative outcome values at time %lf\n", TimeAct); 
	  Warn=1; 
       } 
       *(WorkValArray+k) = EstimateArray[k]; 
       LastTimes[3][k]=(*(WorkValArray+k)); 
    } 
 
    GetDerivatives((long double *)WorkValArray,(long double *) derivs[3],TimeAct); 
    RKCalc=0; 
    TimeAct+=TimeIncrement; 
    *DidIncrement= TimeIncrement; 
    *NextIncrement=TimeIncrement;    
 } 
 
 /* This is an adaptive stepsize ordinary differential equation driver. It         */ 
 /* starts and stops the integration, handles the storage of results and controls  */ 
 /* the stepsize of the integration subject to error conditions,                   */ 
 /* Note that not every step is kept, instead every SAVEFREQ amount of time, the   */ 
 /* currently calculated  value is written to the TimeMatrix which has dimension   */ 
 /* STEPSTOKEEP*(NUMEQ+1) therefore, we keep up to STEPSTOKEEP+1 steps, 
 TempValArray    */ 
 /* actually holds the current solution values                                     */ 
 
 
 void ODEStepper(long double TimeEnd) 
 { 
    int            i,numkept=0,step;        /* numkept = # points have been saved to the TimeMatrix ,  
					       step = # iterations have been done by the method, */ 
    long double    ScaleArray[NUMEQ];   							 
    long double 	  DidIncrement=TimeIncrement;    /*Increment done in last step */ 
    long double    NextIncrement=TimeIncrement;   /*Estimated increment for next step */ 
    long double 	  TempValArray[NUMEQ];             /* Holds solution values for current step, may be stored or not */ 
    int            NumBadSteps = 0, NumGoodSteps=0;  
    long double    LastDiff, Diff=0.0;  
    int            SaveLocalOutCome = 0;               /* Indicator for writing current time point to TimeMatrix because it is the desired outcome time point */ 
    long double derivArray[NUMEQ] ; 
    long double TimeDerivArray[NUMEQ];   
    int First=0; 
    int OutComesSaved=0; 
 
   TimeAct=StartTime; 
   GetDerivatives((long double *) TimeMatrix[0], (long double *) derivArray, TimeAct); 
    if ((Method==3) || ((Method==6) || (Method==7)))GetTimeDeriv(1,(long double *) TimeMatrix[0], (long double *) derivArray, (long double *) TimeDerivArray); 
    TimeIncrement =(long double) INCREMENT; 
    if (Method==2) LastStep=STEPSTOKEEP+5; 
    else LastStep=STEPSTOKEEP+1; 
    if ((Method==3) || (Method==6))  FiniteDiffJacob(NumEq,(long double *) TimeMatrix[0],(long double *) derivArray, 
    (long double *)Jacobian[0]); 
   CalculatedOut(0,StartTime); 
   memset(TempValArray,0,(NUMEQ)*sizeof(long double)); 
   for(i=0;i<(NUMEQ);i++) TempValArray[i]=TimeMatrix[0][i]; 
   for(step=0;step<MAXSTEPS;step++) 
   {  
      GetDerivatives((long double *) TempValArray,(long double *) derivArray,TimeAct); 
      for (i=0;i<NUMEQ;i++) ScaleArray[i] =fabs(TempValArray[i])+fabs(*(derivArray+i))*TimeIncrement + TINY2; 
      if (DoLocal==1)   /* pick out the point closest to the desired local outcome to save to OutcomeData array*/ 
      { 
	 LastDiff=Diff; 
	 Diff=OutDay-TimeAct; 
 
	 if (First==0) if(Diff+TINY<=0.00+TINY) if(LastDiff+TINY>=0.0+TINY)  
	 { 
	    SaveLocalOutCome=1;  
	    First=1; 
	 }  
      }   
      if ((Method==2) && (step<4))  
      { 
	 numkept++; 
	 TimeMatrix[numkept][NUMEQ]=TimeAct; 
	 if (KEEPDISPLAY==1) printf("Current calculation time: %lf Step:%d\n",TimeAct,step); 
	 for (i=0;i<NUMEQ;i++) TimeMatrix[numkept][i]=(*(TempValArray+i)); 
	 CalculatedOut(numkept,TimeAct); 
      } 
      if (  (STEPSTOKEEP >0) && (numkept <(STEPSTOKEEP-1)) ) 
      if ((  (TimeAct >= (numkept+1)*SaveFreq+StartTime) || (SaveLocalOutCome==1)  ))   /* if we've gone SaveFreq time units since last save OR we want to save the current point to the outcome array */ 
      { 
	 numkept++; 
	 TimeMatrix[numkept][NUMEQ]=TimeAct; 
	 if (KEEPDISPLAY==1) printf("Current calculation time: %lf Step: %d\n",TimeAct,step); 
	 for (i=0;i<NUMEQ;i++) TimeMatrix[numkept][i]= (*(TempValArray+i)); 
	 CalculatedOut(numkept,TimeAct);  
	 if ((DoGlobal==1) || ((SaveLocalOutCome==1) && (OutComesSaved==0)) ) SaveOutComeData(numkept,sim); 
	 if (SaveLocalOutCome==1) OutComesSaved++; 
	 SaveLocalOutCome=0;     /* Don't save any more local outcomes */ 
      } 
      current_time = time(0); 
      run_duration = difftime(current_time,start_time); 
      if (run_duration > MAXSIMMIN*60.0) 
      { 
	    sprintf(errmessage,"Maximum time for simulation %d exceeded at integration point : %lf \n",sim,TimeAct); 
	    ErrorSim2(errmessage); 
      }    
      if (( TimeAct + TimeIncrement - TimeEnd) > 0.0) TimeIncrement=TimeEnd -TimeAct; /* if stepsize overshoots the desired end point for integration, correct this*/ 
      if (SimDone==0) switch(Method) 
      { 
	 case 1:RungeKutta(&DidIncrement,&NextIncrement,(long double *) TempValArray,(long double *) derivArray);break; 
	 case 2:   
	 { 
	    if (step<4)  RungeKutta(&DidIncrement,&NextIncrement,(long double *)TempValArray,(long double *) derivArray); 
	    else AdamsBash(&DidIncrement,&NextIncrement,(long double *) TempValArray); 
	    break; 
	 }  
	 case 3:  Rosenbrock((long double *) ScaleArray,&DidIncrement,&NextIncrement,(long double *) TempValArray,(long double *)derivArray,(long double *) TimeDerivArray,EPS);break; 
	 case 4:  RKAdapt((long double *) ScaleArray,&DidIncrement,&NextIncrement,(long double *) TempValArray,(long double *) derivArray,EPS);break; 
	 case 5:  BSStep((long double *) ScaleArray,&DidIncrement, &NextIncrement,(long double *) TempValArray, (long double *) derivArray, (long double *)TimeDerivArray,EPS);break; 
	 case 6:  BSStep((long double *) ScaleArray,&DidIncrement, &NextIncrement,(long double *) TempValArray, (long double *) derivArray, (long double *)TimeDerivArray,EPS);break; 
	 case 7:  Rosenbrock((long double *) ScaleArray, &DidIncrement, &NextIncrement, (long double *) TempValArray, (long double *) derivArray,(long double *) TimeDerivArray,EPS);break; 
	 default: Error("Error in adaptive stepsize routine, unknown solver method chosen");break; 
      } 
      else return; 
      if (DidIncrement==TimeIncrement) NumGoodSteps++;     /*step succeeded at current time step*/ 
      else NumBadSteps++;				  /* needed  a new step size	*/ 
 
      if (TimeAct>=TimeEnd)   		                  /* if we're past the end, save the last step and finish*/ 
      {                                                    /* if we're doing local PRCC, save the outcome */ 
	 numkept++; 
	 TimeMatrix[numkept][NUMEQ]=TimeAct; 
	 for(i=0;i<NUMEQ;i++) TimeMatrix[numkept][i]= (*(TempValArray+i)); 
	 if ((DoLocal==1) && (OutComesSaved<1)) SaveOutComeData(numkept,sim); 
	 CalculatedOut(numkept,TimeAct); 
	 LastStep=numkept+1; 
	 MinSaveNo=numkept; 
	 return; 
      } 
      if ((NextIncrement)<=MININCREMENT) ErrorSim2("Step size in solver too small"); 
      if (DoGlobal==1) if (NextIncrement>0.8*SaveFreq)NextIncrement=0.8*SaveFreq; 
      if (NextIncrement>HMAX) NextIncrement=HMAX; 
      if ((Method>3) && (Method<8)) TimeIncrement=NextIncrement; 
    } 
    if (TimeAct<EndDay) ErrorSim2("\n\n\n\tMaximum Number of steps exceeded before end point \n\tof integration reached. To remedy, increase parameter MAXSTEPSD in file \n\tsetrosen.c and recompile executable\n\n"); 
 
 } 
 
 /* This subroutine performs one step of a Rosenbrock stiff solver. It will repeat the step, updating    */ 
 /* the current stepsize until error bounds are met or the maximum number of trys is exceeded. When      */ 
 /* error bounds are met, the current stepsize is returned as the stepsize that was done (DidIncrement)  */ 
 /* the next stepsize is estimated (NextIncrement), and control returns to the driver routine ODEStepper.*/ 
 /* TimeIncrement is preserved as a global variable due to coding constraints for other olvers.          */ 
 /* CurrentIncrement is the local version of TimeIncrement that is being tried out befowre the step is   */ 
 /* accepted.                                                                                            */ 
 
 void Rosenbrock(long double *ScaleArray, long double *DidIncrement,long double *NextIncrement,long double *WorkValArray,long double 
 *derivArray,long double *timederiv,long double eps)  
 { 
 
    long double 	g1[NUMEQ], 
		 g2[NUMEQ], 
		 g3[NUMEQ], 
		 g4[NUMEQ]; 
 
    long double  IncInverse,		/* Time increment variables*/ 
		 TimeSav,                /* Time at beginning of current iteration*/ 
		 CurrentIncrement  ;     /*  Current trial time increment*/ 
 
    long double 	Err[NUMEQ],          /*Error condition variables*/ 
		 errmax; 
 
    int          i,j,try; 
    long double  temp ;			/* A temporary storage variable*/ 
 
    long double  index[NUMEQ];     /* Used to track the permutation of the row variables*/  
    long double  d;		    /* d tracks the parity of the permutation */											 
    long double  LastTime[NUMEQ];    /* Storage for the last solution and the last derivative solution */ 
    long double  LastDeriv[NUMEQ];    
 
 
    long double   AMatrix[NUMEQ][NUMEQ], 
		  AInverse[NUMEQ][NUMEQ], 
		  ProdMatrix[NUMEQ]; 
 
    TimeSav = TimeAct;					 
    CurrentIncrement=TimeIncrement; 
    IncInverse = (1/CurrentIncrement) ; 
    for (i=0;i<NumEq; i++)  
    { 
       LastTime[i] =(*(WorkValArray+i)); 
       LastDeriv[i]= (*(derivArray+i)); 
    } 
 
    GetTimeDeriv(1,(long double *) LastTime, (long double *) LastDeriv,(long double *) timederiv); 
    if (Method==3) FiniteDiffJacob(NumEq,(long double *) LastTime,(long double *) LastDeriv, (long double *) Jacobian[0]); 
    if (Method==7) FillJacobian(TimeAct);    
 
    for (try=0; try<=MAXTRY; try++)		/* begin current step, and try to meet error bounds for up to MAXTRY attempts*/ 
    { 
       for(i=0;i<NumEq;i++)  
       {				/* calculate I - h*dfdy */ 
	  for(j=0;j<NumEq;j++) AMatrix[i][j] = (-1.0*Jacobian[i][j]); 
	  AMatrix[i][i] += ((1.0)/(GAMMA*CurrentIncrement)); 
       } 
       LUDecomp((long double *) AMatrix[0], NumEq, (int *) index, &d); 
       for (i=0;i<NumEq;i++) g1[i] = LastDeriv[i]+timederiv[i]*CurrentIncrement*C1X; 
       LUBackSub((long double*) AMatrix[0], NumEq, (int *) index, (long double *) g1); 
 
       for (i=0;i<(NumEq);i++) *(WorkValArray+i) = LastTime[i]+A21*g1[i];/*intermediate value and update time*/ 
       TimeAct = TimeSav+ A2X*CurrentIncrement; 
       GetDerivatives((long double *) WorkValArray,(long double *) derivArray,TimeAct); 
       for (i=0;i<NumEq;i++) g2[i]=(*(derivArray+i))+timederiv[i]*CurrentIncrement*C2X +C21*g1[i]/CurrentIncrement; 
       LUBackSub((long double*) AMatrix[0], NumEq,(int *) index,(long double *) g2); 
 
 
       for (i=0;i<NumEq; i++) *(WorkValArray+i) = LastTime[i]+A31*g1[i]+A32*g2[i]; 
       TimeAct= TimeSav+ A3X*CurrentIncrement; 
       GetDerivatives((long double *) WorkValArray,(long double *) derivArray, TimeAct); 
       for (i=0;i<NumEq;i++) g3[i] = (*(derivArray+i))+ CurrentIncrement *C3X* timederiv[i] +(C31*g1[i]+C32*g2[i])/CurrentIncrement; 
 
       LUBackSub((long double*) AMatrix[0], NumEq, (int *)index,(long double *) g3); 
       for (i=0;i<NumEq;i++) g4[i] = (*(derivArray+i))+CurrentIncrement*C4X*timederiv[i]+(C41*g1[i]+C42*g2[i]+C43*g3[i])/CurrentIncrement; 
       LUBackSub((long double *) AMatrix[0], NumEq, (int *)index,(long double *) g4); 
       for (i=0;i<NumEq;i++)   
       { 
	  *(WorkValArray+i)= LastTime[i] + B1*g1[i]+B2*g2[i]+B3*g3[i]+B4*g4[i]; 
	  Err[i] = E1*g1[i]+E2*g2[i]+E3*g3[i]+E4*g4[i]; 
	  if (*(WorkValArray+i)<-TINY) 
	  { 
	      (*(WorkValArray+i))=0.00000; 
	      if (Warn==0) fprintf(warn1, "\tWarning: some equation values were set to zero \n \twhen solver produced negative outcome values\n"); 
	       Warn=1; 
	  } 
       } 
       TimeAct = TimeSav + CurrentIncrement; 
       if (TimeAct==TimeSav)  
       { 
	  fprintf(warn1,"Insignificant stepsize in solver for simulation %d at time point %lf\n",sim,TimeAct); 
	  ErrorSim2( "Insignificant stepsize in solver");  
       }  
 
       errmax=0.0;	/* determine error of the worst offender equation */ 
       for (i=0;i<NumEq; i++)    
       { 
	  temp = fabs(Err[i]/(*(ScaleArray+i))) ; 
	  errmax =  ((errmax > temp)  ? errmax: temp) ; 
       } 
       errmax= errmax/eps; 
 
       if(errmax <=1.0) 		/*if the worst error is acceptable accept the current step, go to next with a larger trial increment if error is small*/ 
       { 
	  *DidIncrement = CurrentIncrement;  
	  *NextIncrement= ((errmax > ERRCON) ?SAFETY*CurrentIncrement*exp(PGROW*log(errmax)):GROW*CurrentIncrement); 
	   memset(LastTime,0,NumEq*sizeof(long double)); 
	   memset(LastDeriv,0,NumEq*sizeof(long double)); 
	   memset(Jacobian,0,NumEq*NumEq*sizeof(long double)); 
	   memset(g1,0,NumEq*sizeof(long double)); 
	   memset(g2,0,NumEq*sizeof(long double)); 
	   memset(g3,0, NumEq*sizeof(long double)); 
	   memset(g4,0, NumEq*sizeof(long double)); 
	   memset(index,0,NumEq*sizeof(int)); 
	   memset(AMatrix,0, NumEq*NumEq*sizeof(long double)); 
	   memset(AInverse,0,NumEq*NumEq*sizeof(long double));   
	   memset(ProdMatrix,0,NumEq*sizeof(long double)); 
	   memset(Err,0,NumEq*sizeof(long double)); 
	  return; 
       }  
       else  /* the truncation error is too large, reduce the stepsize*/ 
       { 
	  *NextIncrement=SAFETY*CurrentIncrement*exp(PSHRNK*log(errmax)); 
	  if (CurrentIncrement+TINY >= 0.0+TINY) CurrentIncrement =  ((*NextIncrement > (SHRNK*CurrentIncrement)) ? *NextIncrement:SHRNK*CurrentIncrement); 
	  else if (CurrentIncrement+TINY<0.0+TINY) CurrentIncrement  = ( (*NextIncrement < (SHRNK*CurrentIncrement)) ?*NextIncrement:SHRNK*CurrentIncrement) ; 
       }  
    }      
    Error("Exceeded Maximum number of trial in subroutine Rosenbrock\n"); 
 }  
 
 
 /* This entire subroutine performs one step of an adaptive stepsize Runge-Kutta Solver */ 
 
 
 void RKAdapt( long double *ScaleArray, long double *DidIncrement, long double *NextIncrement, long double *WorkValArray, long 
 double *derivArray,long double eps) 
 { 
    int         i,j; 
    long double errmax, CurrentIncrement, TempIncrement,TimeNew; 
 
    long double ak2[NUMEQ]; 
    long double ak3[NUMEQ]; 
    long double ak4[NUMEQ]; 
    long double ak5[NUMEQ]; 
    long double ak6[NUMEQ]; 
    long double ytemp[NUMEQ]; 
    long double yerr[NUMEQ]; 
    long double temp; 
    long double EstimateArray[NUMEQ]; 
    long double dc1,dc3,dc4,dc5,dc6; 
    long double   c1 = (37.0/378.0),c3=(250.0/621.0),c4=(125.0/594.0),c6=(512.0/1771.0); 
 
 
    dc1 = c1 - (2825.0/27648.0); 
    dc3 = c3 - (18575.0/48384.0); 
    dc4 = c4-(13525.0/55296.0); 
    dc5 = (-277.00/14336.0); 
    dc6 = c6-0.25; 
 
 /* Array to hold the predicted values at the current trial step, if the error is in bounds the step succeeds and the values for the current working value array (WorkValArray are given the values that were estimated (held in EstimateArray) */ 
    CurrentIncrement=TimeIncrement; 
    for (;;) 
    { 
	for(i=0;i<NumEq;i++) ytemp[i] =(*(WorkValArray+i))+0.2*CurrentIncrement*(*(derivArray+i)); 
	GetDerivatives ( (long double*) ytemp, (long double*) ak2,TimeAct+0.2*CurrentIncrement); 
	for (i=0;i<NumEq;i++) ytemp[i] = (*(WorkValArray+i))+CurrentIncrement*(0.075*(*(derivArray+i))+0.225*ak2[i]); 
	GetDerivatives((long double*) ytemp, (long double *) ak3,TimeAct+0.3*CurrentIncrement); 
	for (i=0;i<NumEq;i++) ytemp[i] = *(WorkValArray+i) +CurrentIncrement*(0.3*(*(derivArray+i))+(-0.9)*ak2[i]+1.2*ak3[i]); 
	GetDerivatives((long double * ) ytemp, (long double *)ak4,TimeAct+0.6*CurrentIncrement); 
	for (i=0;i<NumEq;i++) ytemp[i] = *(WorkValArray+i)+CurrentIncrement*((-11.0/54.0)*(*(derivArray+i))+2.5*ak2[i]+(-70.0/27.0)*ak3[i]+(35.0/27.0)*ak4[i]); 
	GetDerivatives ((long double *) ytemp, (long double *)ak5,TimeAct+CurrentIncrement); 
	for (i=0;i<NumEq;i++)  
	      ytemp[i] =(*(WorkValArray+i))+CurrentIncrement*((1631.0/55296.0)*(*(derivArray+i)) +(175.0/512.0)*ak2[i] + (575.0/13824.0)*ak3[i] + (44275.0/110592.0)*ak4[i]+ (253.0/4096.0)*ak5[i]); 
	GetDerivatives((long double *) ytemp, (long double *)ak6,TimeAct+0.875*CurrentIncrement); 
	for (i=0;i<NumEq;i++) EstimateArray[i]  = (*(WorkValArray+i))+CurrentIncrement*(c1*(*(derivArray+i)) + c3*ak3[i] + c4*ak4[i] +c6*ak6[i]); 
	for (i=0;i<NumEq;i++) yerr[i] =CurrentIncrement*(dc1*(*(derivArray+i))+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]); 
 
	errmax = 0.0; 
	for (i=0;i<NumEq;i++) 
	{  
	    temp = fabs(yerr[i]/(*(ScaleArray+i))); 
	    errmax = ((errmax > temp) ? errmax: temp); 
	} 
	errmax = errmax/eps; 
	if (errmax < 1.0) break;                                                                       /* Step succeeded,  compute size of next step */ 
														      /* Otherwise, update the time increment and recompute the step */ 
	TempIncrement  = SAFETY2 *CurrentIncrement*exp(PSHRNK2*log(errmax));    
	CurrentIncrement = ((0.1*CurrentIncrement<TempIncrement) ? TempIncrement: 0.1*CurrentIncrement); 
	TimeNew=TimeAct+CurrentIncrement; 
	if (CurrentIncrement<MININCREMENT)  
	{ 
	       if (FirstErr==0) printf("The value of the current increment is : %lg\n", CurrentIncrement); 
	       sprintf(errmessage,"Stepsize underflow for simulation %d:Possibly due to a parameter that cause an equation to approach infinity",sim); 
	       if (FirstErr==0) ErrorSim2(errmessage); 
	       FirstErr = 1; 
	}        
    }   
    if (errmax>ERRCON2) *NextIncrement =SAFETY2*CurrentIncrement*exp(PGROW2*log(errmax));  
    else *NextIncrement= 5.0*CurrentIncrement; 
    *DidIncrement = CurrentIncrement; 
    TimeAct+=(*DidIncrement); 
    for (i=0;i<NumEq;i++)  
    { 
	if (EstimateArray[i] < -TINY)  
	{ 
	    (*(WorkValArray+i)) = 0.00; 
	    if (Warn==0) fprintf(warn1, "\tWarning: some equation values were set to zero \n \twhen solver produced negative outcome values at time %lf\n", TimeAct); 
	    Warn=1; 
	} 
	else (*(WorkValArray+i)) = EstimateArray[i];     
    }      
 }    
 
 
 
 
 /* This subroutine provides the values of the derivatives for the solver subroutines. It does so by     */ 
 /* reading in the values of the input array which carries approximated 
 solutions of the NUMEQ and then */ 
 /* calculating the value of the  derivative at the point using these solutions. The actual code is      */ 
 /* written by the setup program and is inserted after the comment line                                  */ 
 
 void GetDerivatives(long double *inarray, long double *outarray, long double t) 
 { 
    int    j; 
    long double simp[MAXPARS];    
 
 /* ASSIGN, CALCULATE SIMPLIFICATIONS AND CALCULATE DIFFERENTIAL EQUATIONS HERE */ 
 
#include "C_Equations"
 
 } 
 
 
 
 /* This subroutine fills in the nonconstant values of the Jacobian matrix when they     */ 
 /* are used for the stiff solvers. The actual code is written by the setup program and  */ 
 /* is inserted after the comment line.                                                  */ 
 
 
 void FillJacobian(long double t) 
 { 
 /*PUT JACOBIAN CALCULATIONS HERE*/ 

#include "C_Jacob"

 } 
 
 
 /*This subroutine calculates outcomes of interest that are a (nonidentity)function of the equation outputs */ 
 
 void CalculatedOut(int numkept,long double t) 
 { 
     int i; 
 
 /*PUT CALCULATED OUTCOMES HERE */
#include "C_OutComeCalcs"
 
 } 
 
 
 /* This subroutine takes the LHSMatrix values that vary and ranks them for the RankMatrix.  It accounts for    */ 
 /* ties in the ranking, by giving tied values the average of the ranks that they would have received if they   */ 
 /* had differed slightly (see Mendenhall, section on rank correlation coefficients for detail)                 */ 
 /* Note as well that the number of values to be sorted for each run depends on whether we are ranking the      */ 
 /* parameter selections  and the outcome of interest for PRCC analysis and printing to the file Check_Ranks    */ 
 /* (as in PRCCMODE=0) , or just rankingthe parameter selections to print to the file Check_Ranks (when         */ 
 /* PRCCMODE =1)                                                                                                */ 
 
 void RankTieSort(void) 
 {	 
     int          count,ties=0,first; 
     int          i,j,k; 
     double  ranksum, newrank; 
     int          NumToSort; 		 
 
     if (DonePRCC==0) NumToSort=NumPar+1;  /* We're doing PRCC analysis */ 
     else NumToSort = NumPar;  /* We're not doing PRCCanalysis*/ 
     for(k=0;k<NumToSort;k++)  
     { 
	if((pArray[k].min != pArray[k].max) | (k==NumPar))         /* if we have a varying parameter or its the outcome */ 
	{ 
	   for(i=0; i<NumSim; i++) LHSMatrix[k][i].rank = 1.0;     /* initialize kth row ranks to 1.0 */ 
	   ties=0; 
	   for(i=1;i<NumSim;i++)                                   
	   { 
	      first=1; 
	      for(j=0;j<i;j++)  
	      { 
		 if(LHSMatrix[k][i].value >= LHSMatrix[k][j].value) LHSMatrix[k][i].rank+=1.0; 
		 else LHSMatrix[k][j].rank+=1.0; 
		 if(LHSMatrix[k][i].value +TINY==LHSMatrix[k][j].value+TINY)  
		 { 
		    if(first==1) ties++; 
		    LHSMatrix[k][i].tieset=ties; 
		    LHSMatrix[k][j].tieset=ties; 
		    first=0; 
		 }     
	      }     
	   } 
	   for(i=1;i<=ties;i++)  
	   { 
	      ranksum=0.0;	    
	      count=0; 
	      for(j=0;j<NumSim;j++) if(LHSMatrix[k][j].tieset==i)  
	      { 
		 count++; 
		 ranksum+=LHSMatrix[k][j].rank; 
	      }  
	      newrank= (ranksum/count); 
	      for(j=0;j<NumSim;j++) if(LHSMatrix[k][j].tieset==i) LHSMatrix[k][j].rank=newrank; 
	   }    
	}    
     }     
  } 
 
 /* This subroutine takes the LHSMatrix values that vary and ranks them for the RankMatrix.  It accounts for  
 ties in the ranking, by giving tied values the average of the ranks that they would have received if they  
 had differed slightly (see Mendenhall, section on rank correlation coefficients for detail) 
 
  void RankTieSort(void) 
  { 
     int count,ties,first; 
     int i,j,k; 
     long double ranksum, newrank; 
 
     for(k=0;k<(K+1);k++)  
     { 
	if((pArray[k].min != pArray[k].max) | (k==K)) 
	{ 
	   for(i=0; i<N; i++) LHSMatrix[k][i].rank = 1; 
	   for(i=1;i<N;i++)  
	   { 
	      first=1; 
	      for(j=0;j<i;j++)  
	      { 
		 if(LHSMatrix[k][i].value >= LHSMatrix[k][j].value) LHSMatrix[k][i].rank+=1; 
		 else LHSMatrix[k][j].rank+=1; 
		 if(LHSMatrix[k][i].value ==LHSMatrix[k][j].value)  
		 { 
		    if(first==1) ties++;               
		    LHSMatrix[k][i].tieset=ties; 
		    LHSMatrix[k][j].tieset=ties; 
		    first=0; 
		 }     
	      }     
	   } 
	   for(i=1;i<=ties;i++)  
	   { 
	      ranksum=0.0;           
	      count=0; 
	      for(j=0;j<N;j++) if(LHSMatrix[k][j].tieset==i)  
	      { 
		 count++; 
		 ranksum+=LHSMatrix[k][j].rank; 
	      }  
	      newrank= (ranksum/count); 
	      for(j=0;j<N;j++) if(LHSMatrix[k][j].tieset==i) LHSMatrix[k][j].rank=newrank; 
	   }    
	}    
     }     
  } 
 */ 
 
 /*This subroutine outputs the sample of K parameters for each of the N runs to the matrix Check_Param*/ 
 
 void PrintParamSamples(void) 
 { 
      int k, l; 
 
      fprintf(check[2], "\n"); 
      if( NumSim<=9) 
      { 
	 for(k=0; k<NumPar; k++) 
	 { 
	    fprintf(check[2],"\n\n"); 
	    fprintf(check[2],"%s\t",pArray[k].name);  
	    for(l=0; l<NumSim; l++)  
	    { 
		if(LHSMatrix[k][l].value > 0.00001) fprintf(check[2],"%lf ",LHSMatrix[k][l].value); 
		else fprintf(check[2],"%lg  ",LHSMatrix[k][l].value); 
	    } 
 
	 } 
      } 
      if (NumSim>9)  
      { 
	 for(k=0; k<NumSim; k++) 
	 { 
	    fprintf(check[2],"\nRun: %d \n",k); 
	    for(l=0; l<NumPar; l++) 
	    { 
	       if (LHSMatrix[l][k].value >0.00001)fprintf(check[2],"%s:%lf\t",pArray[l].name,LHSMatrix[l][k].value); 
	       else fprintf(check[2],"%s:%lg\t",pArray[l].name,LHSMatrix[l][k].value); 
	       if(l!=0 && l%4==0) fprintf(check[2] ,"\n"); 
	    } 
	    fprintf(check[2], "\n"); 
	 } 
     }      
     fprintf(check[2],"\n"); 
 } 
 
 
 
 /* Uncertainty Analysis :  This subroutine calculates the mean, median, and variance for the PRCC outcome */ 
 /* variable. Results are output to the file Check_Uncert                                                  */ 
 
 void UncertaintyAnalysis(void) 
 { 
      long  double total, mean, median, var,numer, max,min; 
      int   i; 
 
      total = 0.0; 
      numer = 0.0; 
      max = MINDOUBLE; 
      min= MAXDOUBLE; 
      median = LHSMatrix[NumPar][(NumSim+1)/2].value; 
      for(i=0;i<NumSim;i++)  
      { 
	   total+=LHSMatrix[NumPar][i].value; 
	   numer += (long double) (LHSMatrix[NumPar][i].value *LHSMatrix[NumPar][i].value); 
	   if (LHSMatrix[NumPar][i].value >max) max=LHSMatrix[NumPar][i].value; 
	   if (LHSMatrix[NumPar][i].value < min) min=LHSMatrix[NumPar][i].value; 
      } 
      mean = (long double) total/NumSim; 
      var = (numer -NumSim*mean*mean)/(NumSim-1); 
 
      if ((DoGlobal==1) && (KEEPTDUNCERT==1) ) 
      { 
	 if (Call==1)fprintf(check[4],"Time\t\tMean\t\tMedian\t\tVariance\t\tMaximum\t\tMinimum\n"); 
	 fprintf(check[4],"%lf \t %lg \t %lg \t %lg \t %lg \t %lg\n",CurrOutTime,mean,median,var,max,min);  
      }  
      if (DoLocal==1) fprintf(check[4],"\tMean:\t\t%lg\n\tMedian:\t\t%lg\n\tVariance:\t%lg\n\tMaximum:\t%lg\n\tMinimum:\t%lg\n\n",mean,median,var,max,min );  
      if ((DoGlobal==1) && (KEEPTDUNCERT==0)) fprintf(check[4],"Data for time-dependent uncertainty Analysis not kept\n"); 
 } 
 
 
 /*Begin Sensitivity Analysis */ 
 
 
 /*This subroutine copies the ranks for the parameters that are nonconstant in the LHS matrix into       */ 
 /* the matrix RankMatrix (dimension NumVarPar+1 * N)  where NumVarPar=the number of parameters that are nonconstant*/ 
 /* This stuff needs to be moved to the parameter reads  */ 
 
 
 void CopyVaryPar(void) 
 { 
      int  i, j, m,n; 
 
      NumVarPar=0; 
      for(i=0; i<NumPar; i++) 
      { 
	   if(pArray[i].min != pArray[i].max) 
	   { 
	      for(j=0; j<NumSim; j++) RankMatrix[NumVarPar][j] = LHSMatrix[i][j].rank; 
	      NewI[NumVarPar]=i; 
	      NumVarPar++;               
 
	      if (Call==1)  
	      { 
		  if (NumVarPar==1) fprintf(check[5],"\n\nParameters that vary:\n \n"); 
		  if((pArray[i].min <100.0) && (pArray[i].min>0))fprintf(check[5],"%s\t\tmin:%lg\tmax:%lg\n",pArray[i].name,pArray[i].min,pArray[i].max); 
		  else if (pArray[i].min>=100)fprintf(check[5],"%s\t\tmin:%lg\tmax:%lg\n",pArray[i].name,pArray[i].min,pArray[i].max); 
		  else if (pArray[i].min==0) fprintf(check[5],"%s\t\tmin:%5.4lf\tmax:%lg\n",pArray[i].name,pArray[i].min,pArray[i].max); 
		  fprintf(check[1]," \n Varying parameter: %s \n",pArray[i].name); 
		  for(j=0; j<NumSim; j++) fprintf(check[1], " %f \t %lg\n",RankMatrix[NumVarPar-1][j],LHSMatrix[i][j].value); 
		  if (Call==1) printf("%s copied for PRCC\n",pArray[i].name); 
	      } 
	   } 
      } 
 
      if (PRCCOffTime>=CurrOutTime) 
      { 
	 for(m=0; m<NumSim; m++) RankMatrix[NumVarPar][m] = LHSMatrix[NumPar][m].rank; 
	 if ((DoLocal==1) || ((DoGlobal==1)&& (KEEPTDRANKS==1)))  
	 { 
	    fprintf(check[1],"\n \n Outcome Values:\n"); 
	    for (m=0;m<NumSim;m++) fprintf(check[1], " %f \t%lg\n",RankMatrix[NumVarPar][m],LHSMatrix[NumPar][m].value); 
	    fprintf(check[1]," \n \nRank matrix: \n"); 
	    PrintMatrix(check[1],(double *)RankMatrix[0],NumVarPar+1,NumSim);  
	 }  
     } 
     for (i=1;i<6;i++) fflush(check[i]); 
 
 } 
 
 
 
 
 
 /*Using the matrix of parameters that vary, this subroutine calculates the Cij values for the which are  */ 
 /* used to calculate the PRCCs. We obtain the ijth entry of the cijMatrix via the formula:               */ 
 
 /* Cij = ((\sum_t(r_it-\mu) (r_jt-\mu)) \over (\sum_t((r_it-\mu)^2) \sum_s ((r_js - \mu)^2)  ) ^(1/2)     */ 
 /*    = ((\sum_t (r_it r_jt) - \mu(\sum_t r_it) - \mu (\sum_t r_jt) + N \mu^2)                            */ 
 /*       \over (\sum_t(r_it^2 - 2 \mu \sum r_it + N \mu^2) )^2)^1/2                                       */ 
 /*    = ((\sum_t r_it r_jt) - 2\mu (\sum_t r_it) + N \mu^2)                                               */ 
 /*        \over ( \sum_t(r_it^2) - 2 \mu \sum r_it + N \mu^2) ) ^2 )^1/2                                  */ 
 /*    = ((\sum_t r_it r_jt) - 2 \mu (N\mu) + N \mu^2) \over (\sum_t(r_it^2 - 2 \mu N \mu  + N \mu^2) )    */ 
 /*    = ((\sum_t r_it r_jt) - N \mu^2) \over (\sum_t (r_it^2 - N \mu^2)                                   */ 
 /*                                                                                                        */      
 /* since                                                                                                  */ 
 /*          \mu           = {1+2+... N \over N} = {N(N+1) \over 2N}                                       */ 
 /*                        = {N+1 \over 2}                                                                 */ 
 /* and                                                                                                    */ 
 /*          \sum_t r_it^2 = 1^2+2^2+... N^2                                                               */ 
 /*                        = {N(N+1)(2N+1) \over 6 } )                                                     */ 
 /*                                                                                                        */ 
 /* Cij = {( (\sum_t r_it r_jt) - N \mu^2) \over (N(N+1)(2N+1) /6 - N (N+1)^2/4 )}                         */ 
 /*    = {( (\sum_t r_it r_jt) - N \mu^2) \over (N(N+1)/2) ( 2N+1/3 - N+1/2) }                             */ 
 /*    = {( (\sum_t r_it r_jt) - N \mu^2) \over (N(N+1)/2) (N-1/6) }                                       */ 
 /*    = { (12 (\sum_t rit rjt ) - N \mu^2) \over  (N(N+1)(N-1)) }                                         */ 
 
 void FillCijMatrix(void) 
 { 
      long double      mu, sum[3], denominator; 
      int              i,j,t,ll; 
 
      mu = (long double) (1.0+NumSim)/2.0; 
      for(i = 0; i< (NumVarPar+1); i++) 
      { 
	 for(j=0; j< (NumVarPar+1); j++) 
	 { 
	   for(ll=0;ll<3;ll++) sum[ll] = 0.0; 
	    for(t = 0; t< NumSim; t++) 
	    { 
		sum [0] += (RankMatrix[i][t])*(RankMatrix[j][t]); 
		sum[1] += (RankMatrix[i][t])*(RankMatrix[i][t]); 
		sum [2] += (RankMatrix[j][t])*(RankMatrix[j][t]); 
	    } 
	    for(ll=0;ll<3;ll++) sum[ll]-=NumSim*mu*mu; 
	    denominator = sqrt(sum[1]*sum[2]);                
	    if (denominator != 0.0) CijMatrix[i][j] = sum[0]/denominator; 
	    else 
	    { 
	       printf("Error: Zero denominator at outcome time :%lf in Fill Cij\n"); 
	       Invert=0; 
	    }  
	 } 
      }      
      if ((DoLocal==1) || ((DoGlobal==1) && (KEEPTDRANKS==1))) 
      { 
	 fprintf(check[1],"Cij Matrix:\n"); 
	 PrintMatrix(check[1],(long double *)CijMatrix[0],NumVarPar+1,NumVarPar+1); 
      } 
 } 
 
 
 
 /* This subroutine prints out the Cij matrix, the inverted matrix and checks the inversion  */ 
 /* by multiplying the CijMatrix by CijInverse and then checking to see if the off-diagonal  */ 
 /* entries of the product matrix (CheckMatrix) are reasonably close to 0 and if the diagnal */ 
 /* entries are close to 1.0. If any entry does not meet this criterion, Invert is set to 0  */ 
 /* flagging that the matrix inversion was not successful.                                   */ 
 /* If inversion is successful, CijMatrix and CijInverse are written out to Check_Inverse    */ 
 /* If it is not, we write an error message to Check_Inverse and check for a strong effect   */ 
 /* parameter.                                                                               */ 
 
 void CheckInverse (void) 
 { 
   int          i, j, k; 
   long double CheckMatrix[MAXVARPARS][MAXVARPARS]; 
 
      memset(CheckMatrix,0,(MAXVARPARS)*(MAXVARPARS)*sizeof(long double)); 
 
      MatMult((long double *) CijMatrix[0],(long double *)CijInverse[0],(long double *)CheckMatrix[0],NumVarPar+1,NumVarPar+1,NumVarPar+1); 
      Invert=1; 
      for(i=0; i<(NumVarPar+1);i++) 
      { 
	 for(j=0;j<(NumVarPar+1);j++) 
	 {    
	    if (i!=j) if ((-.000001>CheckMatrix[i][j]) || (CheckMatrix[i][j]>0.000001))Invert=0; 
	    if (i==j) if ((0.9999>CheckMatrix[i][j]) || (CheckMatrix[i][j]>1.000001)) Invert=0; 
	 } 
      } 
      if (Invert==1) printf("Matrix inversion successful\n"); 
      else CheckStrongEffect(); 
      if ((DoLocal==1) || ((DoGlobal==1) && (KEEPTDINVERT==1))) 
      { 
	 fprintf(check[3], "\n \nCij Matrix:\n  "); 
	 PrintMatrix(check[3],(long double *)CijMatrix[0],NumVarPar+1,NumVarPar+1); 
	 if (Invert==1)  
	 { 
	    fprintf(check[3],"\n \nInverted Cij Matrix: \n "); 
	     PrintMatrix(check[3],(long double *)CijInverse[0],NumVarPar+1,NumVarPar+1); 
	     fprintf(check[3], "\n \nCheck that product produces the identity:\n\n ") ; 
	     PrintMatrix(check[3],(long double *)CheckMatrix[0],NumVarPar+1,NumVarPar+1); 
	 } 
	 else fprintf(check[3], "\n\nSingular Matrix for outcome %s at time %lf\n",OutName,CurrOutTime); 
     }  
 } 
 
 
 /* Using the inverted CijMatrix, this procedure calculates the PRCC valuees of each parameter with  */  
 /* the chosen outcome variable and an approximate t-test of the value. It then outputs the results  */ 
 /* to a variety of different PRCC files,depending on the type of analysis that is being done        */ 
 /* The fundamental types of analyses are either a point PRCC (when DoLocal=1) or a time-dependent   */ 
 /* PRCC (do Global==1) that is done for each point of the PRCC outcome variable that is output to   */ 
 /* the data files (ie every SaveFreq time units). When a local PRCC is done, we output the values   */ 
 /* to the file Check_PRCC. When a time-dependent PRCC is done, we may elect to keep the data in     */ 
 /* Check_PRCC (KeepTDPRCC=1) or not (KeepTDPRCC=0). A separate set of data files is also created    */ 
 /* when the time-dependent PRCC is done; these will be denoted by PRCC_all_* where * is the name of */ 
 /* the PRCC outcome variable. This data file is then used by gnuplot and the setup program to create*/ 
 /* the file PRCC_all_*.ps or PRCC_all_*.eps (depending on current graph option .                    */ 
 /* This subroutine also accounts for the possibility that the CijMatrix was singular and could not  */ 
 /* be inverted. When this happens, (either during the matrix inversion, or during the subroutine    */ 
 /* CheckInverse(), CheckStrongEffect() is called to screen the rank matrix (and thus the CijMatrix )*/ 
 /* for a row that is the same or exactly opposite to the row corresponding to the outcome variable's*/ 
 /* ranks. When this occurs, Invert is set to 0 , flagging that the matrix inversion was not         */ 
 /* successful, and the the indicator variable StrongEffect is set to 1; meaning a parameter was     */ 
 /* found to have either correlation of 1 or -1 with the ranks of the outcome variable. StrongPar    */ 
 /* denotes the number of the parameter that was found to have a strong effect on the outcome. So    */ 
 /* then in this case, we must set the PRCC for the parameter pArray[StrongPar] equal to 1.0 or -1.0 */ 
 /* and set all other PRCC values to be nonexistent (-). This is necessary so that the graphs will   */ 
 /* not be truncated by gnuplot.                                                                     */ 
 
 /*                                     COMMENT LINES                                                    */ 
 /* C1: We loop over all of the rows of the CijInverse; calculating the PRCC for the ith                 */ 
 /*     varying parameter                                                                                */ 
 /* C2: If the matrix inversion was successful, and we didn't have a strong effect parameter             */ 
 /*     the PRCC formula is:                                                                             */ 
 /*       PRCC[i] = -CijInverse[i][NumVarPar]/ sqrt( CijInverse[NumVarPar][NumVarPar]*CijInverse[i][i] ) */ 
 /* C3: Otherwise, if the inversion was not successful, and we have a strong effect variable             */ 
 /*     AND this is that variable, set PRCC[i]=1.0 or -1.0 as appropriate. If this is not the strong     */ 
 /*     effect parameter, set it equal to missing.                                                       */ 
 /* C4: If we're doing a pointwise PRCC, or we're doing a time-dependent PRCC and we want to save the    */ 
 /*     PRCC values in a table format, (KeepTDPRCC=1), then append the values to the file Check_PRCC     */ 
 /* C6: If we're doing time-dependent PRCC, write the data out to the file PRCC_all_* for the outcome    */ 
 /*     variable *                                                                                       */ 
 /* C7: If we're on the first varying parameter, open the file PRCC_all_*  write the time and the value  */ 
 /*     of the PRCC for this parameter, don't start a new line                                           */ 
 /* C8: If we're not on the first and we're not on the last , just write the PRCC for the current par    */ 
 /* C9: If we're on the last varying parameter, write the PRCC and end the row                           */ 
 
 void GetPRCC(void) 
 { 
    long double PRCCNum, PRCCDenom,PRCC,TTest ; 
    int i;  
 
    if ((DoGlobal==1) && (Call==1) )fprintf(check[5],"\nIf PRCC=1.0 or-1.0,t-test is effectively infinite\n"); 
 
   if ((Invert==1) || (StrongEffect==1))  for(i=0; i< NumVarPar; i++)                  /* 1 */ 
   { 
      PRCC=0.0; 
      if ((StrongEffect !=1) && (Invert==1))            /*2*/ 
      { 
	 PRCCNum = -1*(CijInverse[i][NumVarPar]); 
	 if((CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar])< 0.000) 
	 {    
	    if((CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar])<-0.01)Error("PRCC denominator <0 in GetPRCC-no square root"); 
	    else 
	    { 
	       PRCCDenom=-1*(CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
	       fprintf(warn1, "%g \t %g \n", PRCCNum,CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
	       printf("PRCC Denom close to zero and negative, opposite sign given, check warning file"); 
	    } 
	 } 
 
	 PRCCDenom= sqrt(CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
	 (PRCCDenom!=0.0000) ?  PRCC = PRCCNum/(PRCCDenom):   Error("PRCCDenom is zero in GetPRCC"); 
	 (-1.0<=PRCC<=1.0) ? TTest = PRCC*(sqrt((NumSim-2)/(1-PRCC*PRCC))):Error("Error in Calculating T-test Value"); 
      } 
      if ((StrongEffect==1) && (i==StrongPar)) PRCC=(long double) Parity*1.0;   /* 3 */ 
 
      if ((DoLocal==1) ||((DoGlobal==1) && (KEEPTDPRCC==1)))          /* 4 */ 
      { 
	 if ((i==0) && (Call==1)) fprintf(check[5],"\nT-test based on:%d degrees of freedom\n\n",NumSim-2); 
	 if ((i==0) && (DoGlobal==1)) fprintf(check[5],"\n\nOutcome time:%lf\n",CurrOutTime); 
	 if (((StrongPar==i) || (StrongEffect!=1)) && (PRCC>=0.00)) fprintf(check[5], "%s\t %lf\t%lf\n",pArray[NewI[i]].name,PRCC,TTest); 
	 if (((StrongPar==i) || (StrongEffect!=1)) && (PRCC<0.00)) fprintf(check[5],"%s\t%lf\t%lf\n",pArray[NewI[i]].name, PRCC,TTest); 
      } 
 
      if (DoGlobal==1) /* 6 */ 
      { 
	 if (i==0)            /* 7 */ 
	 { 
	    if (FirstCall==1)  
	    { 
	       sprintf(line,"PRCC_all_%s",OutName); 
	       prcc[0]= fopen(line,"w"); 
	       FirstCall=0; 
	    } 
	    if (StrongEffect==0) fprintf(prcc[0],"%lf  %lf  ",CurrOutTime,PRCC); 
	    if (NumVarPar >1) 
	    { 
	      if ((StrongEffect==1) && (i==StrongPar)) fprintf(prcc[0],"%lf %lf  ",CurrOutTime,PRCC); 
	      if ((StrongEffect==1) && (i!=StrongPar)) fprintf(prcc[0],"%lf  -   ",CurrOutTime); 
	    } 
	    else  
	    { 
	      if ((StrongEffect==1) && (i==StrongPar)) fprintf(prcc[0],"%lf %lf  \n",CurrOutTime,PRCC); 
	      if ((StrongEffect==1) && (i!=StrongPar)) fprintf(prcc[0],"%lf  -   \n",CurrOutTime); 
	    } 
	 } 
	 if ((i>0) && (i<NumVarPar-1))                 /* 8 */ 
	 { 
	    if (StrongEffect==0) fprintf(prcc[0],"%lf   ",PRCC); 
	    if ((StrongEffect==1) && (i==StrongPar)) fprintf(prcc[0],"%lf  ",PRCC); 
	    if ((StrongEffect==1) && (i!=StrongPar)) fprintf(prcc[0]," -  "); 
	 } 
	 if ((i==NumVarPar-1) && (i!=0))                           /* 9 */ 
	 { 
	    if (StrongEffect==0) fprintf(prcc[0],"%lf\n",PRCC); 
	    if ((StrongEffect==1) && (i==StrongPar)) fprintf(prcc[0],"%lf\n",PRCC);    
	    if ((StrongEffect==1) && (i!=StrongPar)) fprintf(prcc[0]," - \n"); 
	 } 
      } 
    }   
 /*   if (DoLocal==1)  
    { 
       fprintf(check[5], "\n\n"); 
       if (StrongEffect==1) fprintf(check[5],"\tNote: Values with PRCC=1.0 or -1.0 have extreme\n\teffect on outcome,t-test is effectively infinite\n\n\n"); 
    } */ 
    else 
    { 
       if (FirstCall==1)  
       { 
	  sprintf(line, "PRCC_all_%s",OutName); 
	  prcc[0]=fopen(line,"w"); 
	  FirstCall=0; 
       } 
       fprintf(prcc[0],"%lf  ",CurrOutTime); 
       for (i=0;i<NumVarPar-1;i++) fprintf(prcc[0]," -"); 
       fprintf(prcc[0]," - \n"); 
    } 
 
 } 
 
 
 /* This subroutine creates the log file (Log_File) that is appended each time the program is run. It contains information    */ 
 /* about the system inputs, equations, parameter values, time of run, solution method , type of PRCC done, type of equations */ 
 /* the version of the code used for the run, and outcome of the single matrix inversion if a local PRCC is done              */ 
 
 void MakeLog(void) 
 {  
      int      i; 
 
      fprintf(warn1,"%s\n",errmessage); 
      fprintf(warn1, "_________________________________________________________________________\n\n"); 
      fprintf(warn1, "Information for current run:\n\n"); 
      fprintf(warn1, "\tDate: \t\t\t\t%s", c); 
      switch(Method) 
      { 
	 case 1:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Runge-Kutta"); break;  
	 case 2:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Adams-Bashforth");break; 
	 case 3:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Adaptive Stiff Rosenbrock");break; 
	 case 4:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Adaptive Stepsize RungeKutta");break; 
	 case 5:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Adaptive Stepsize Burlisch-Stoer");break; 
	 case 6:  fprintf(warn1, "\tSolver Method:\t\t\t%s\n", "Stiff Burlisch-Stoer");break; 
	 case 7:  fprintf(warn1,"\tSolver Method:\t\t\t%s\n","Rosenbrock - Jacobian precalculated"); break; 
	 default: fprintf(warn1, "\tSolver Method:\t\t\t%s\n","Unknown\n");break; 
      }      
      fprintf(warn1, "\tName of parameter file: \t%s \n", parfile); 
      if ((Method==1) || (Method==2))fprintf(warn1, "\tTime increment:\t\t\t%f\n",(double) INCREMENT); 
      if ((Method>3) && (Method<8)) fprintf(warn1,"\tTolerance:\t\t\t%lg\n",EPS); 
      fprintf(warn1, "\tNumber of runs: \t\t%d\n", NumSim); 
      if (DoGlobal==0) fprintf(warn1, "\tInversion successful:\t\t%s\n",(Invert==1) ? "Yes": "No");  
      fprintf(warn1,"\tPRCC Type:\t\t\t%s\n",((DoGlobal==1)?"Time-dependent":"Single Point")); 
      fprintf(warn1,"\tPRCC Variable:\t\t\t%s\n",OutName); 
      fprintf(warn1,"\tEquation type:\t\t\t%s\n",((StartTime<TINY)&& (EndDay>ENDDAY-2.0)) ? "Continuous":"Discontinuous"); 
      fprintf(warn1, "\tCode last updated: \t\t11/11/00\n\n\n"); 
      fprintf(warn1, "Parameters: \n\n"); 
      for(i=0; i<NumPar; i++)  
      { 
	 if(pArray[i].min <=999.0) fprintf(warn1,"\t%s\t\t%lg\t0.0\t%lg\t%c\n",pArray[i].name,pArray[i].min,pArray[i].max,pArray[i].dType); 
	 else fprintf(warn1,"\t%s\t\t%lg\t0.0\t%lg\t%c\n",pArray[i].name, pArray[i].min,pArray[i].max,pArray[i].dType); 
      } 
      fprintf(warn1,"\n\n\n"); 
      printf("\n\nValues for parameters for each run written to file Check_Param \n"); 
      printf("Parameter ranks written to file Check_Ranks\n"); 
      if(Invert==1) printf("PRCC values written to file Check_PRCC\n"); 
      if (Invert==1) printf("Matrix inversion results written to file Check_Inverse \n"); 
      printf("Values for Uncertainty Analysis appended to file Check_Uncert\n"); 
      printf("Error messages and run log appended to file Log_File\n\n"); 
 
 }	 
 
 /* This subroutine writes equations solutions to various output files.  Two types of file are produced:                */ 
 /*    (1)   A file called Simulation* where * denotes the simulation number. It contains data for                      */ 
 /*          all equation solutions and all outcome variables that are calculated from these solutions                  */ 
 /*          and is in the format:                                                                                      */         
 /*            time0   Eqn0(time0) ... Eqn_NumEq-1(time0) CalcOut_0 (time0) ... CalcOut_NumCalcOut-1(time0)              */ 
 /*            time0   Eqn0(time1) ... Eqn_NumEq-1(time1) CalcOut_0 (time1) ... CalcOut_NumCalcOut-1(time1)              */ 
 /*                                                                                                                      */ 
 /*                                                                                                                      */ 
 /*                                                                                                                      */ 
 /*            timeEnd   Eqn0(timeEnd) ... Eqn_NumEq-1(timeEnd) CalcOut_0 (timeEnd) ... CalcOut_NumCalcOut-1(timeEnd)    */ 
 
 /*     (2)  For each of the equations and each of the calculated outcomes, a file of the same name is produced with     */ 
 /*          data written in the format:                                                                                 */ 
 /*                                                                                                                      */ 
 /*                     time0     Eqn(time0)               (start simulation 0 data)                                     */ 
 /*                     .                                                                                                */ 
 /*                     .                                                                                                */ 
 /*                     timeEnd   Eqn(timeEnd)             (end simulation 0 data)                                       */ 
 /*                     time0     Eqn(time0)               (start simulation 1 data)                                     */ 
 /*                     .                                                                                                */ 
 /*                     timeEnd   Eqn(timeEnd)             (end simulation 1 data )                                      */ 
 /*                     ....                                                                                             */ 
 /*          The first file type is used for making the graphs in gnuplot. The second is used if we exceed the number    */ 
 /*          of simulations that can be represented in color in gnuplot.  Also, these files have been kept in case       */ 
 /*          because they can be more easily manipulated for other types of data representation, and because             */ 
 /*          it is easier to look at data for individual equations in this way instead of tracking the columns in  a     */ 
 /*          huge output file.                                                                                           */ 
 
 
 
 
 
 void MakeGraphData(void) 
 { 
   int i,j; 
   FILE *data; 
   long double max, min, maxtime, mintime; 
   long double MaxOutValues[NumCalcOut], MinOutValues[NumCalcOut]; 
   long double MaxOutTimes[NumCalcOut], MinOutTimes[NumCalcOut]; 
   long double MaxEqValues[NumEq], MinEqValues[NumEq]; 
   long double MaxEqTimes[NumEq], MinEqTimes[NumEq]; 
 
   long double epsilon = 0.0000000000001; 
   sprintf(line,"Simulation%d",sim); 
   if (StartTime==0) data = fopen(line,"w"); 
   else data=fopen(line,"a"); 
   for(j=0;j<LastStep;j++) 
     { 
       fprintf(data, "%lg",TimeMatrix[j][NUMEQ]); 
       for(i=0;i<NUMEQ;i++) fprintf(data,"    %lg",TimeMatrix[j][i]); 
       for(i=0;i<NUMCALCOUT;i++) fprintf(data,"    %lg",CalcTime[j][i]); 
       fprintf(data,"\n"); 
     } 
   fflush(data); 
   fclose(data); 
 
   printf("%d\n", NumEq); 
 
   for (i=0;i<NumEq;i++) 
     { 
       sprintf(line,"%s",Eqns[i].name); 
       if ((StartTime==0) && (sim==0)) data=fopen(line,"w"); 
       else data=fopen(line,"a"); 
 
       max = TimeMatrix[0][i];   /* set max and min to initial values */ 
       min = TimeMatrix[0][i]; 
       maxtime = StartTime; 
       mintime = StartTime; 
 
       for(j=0;j<LastStep;j++) 
	 { 
	   fprintf(data,"%lg %lg\n",TimeMatrix[j][NUMEQ],TimeMatrix[j][i]); 
	   if (TimeMatrix[j][i] > max) { /* if value is larger than old max, reset max */ 
	     max = TimeMatrix[j][i]; 
	     maxtime = TimeMatrix[j][NUMEQ]; 
	   } 
	   if (TimeMatrix[j][i] < min) { 
	     min = TimeMatrix[j][i]; 
	     mintime = TimeMatrix[j][NUMEQ]; 
	   } 
	 } 
       fflush(data); 
       fclose(data); 
       MaxEqValues[i] = max; 
       MinEqValues[i] = min; 
       MaxEqTimes[i] = maxtime; 
       MinEqTimes[i] = mintime; 
        
       if (max > TopValues[i]) 
	 TopValues[i] = max; 
     } 
 
   for (i=0;i<NumCalcOut;i++) { 
     sprintf(line,"%s",CalcOut[i].name); 
     if ((StartTime==0) && (sim==0)) data=fopen(line,"w"); 
     else data=fopen(line,"a"); 
 
     max = CalcTime[0][i];   /* set max and min to initial values */ 
     min = CalcTime[0][i]; 
     maxtime = StartTime; 
     mintime = StartTime; 
 
     for(j=0;j<LastStep;j++) 
       { 
	 fprintf(data,"%lg  %lg\n",TimeMatrix[j][NUMEQ],CalcTime[j][i]); 
	 if (CalcTime[j][i] > max)  {/* if current value is larger than old max, reset max */ 
	   max = CalcTime[j][i]; 
	   maxtime = TimeMatrix[j][NUMEQ]; 
	 } 
	 if (CalcTime[j][i] < min) { 
	   min = CalcTime[j][i]; 
	   mintime = TimeMatrix[j][NUMEQ]; 
	 } 
       } 
     fflush(data); 
     fclose(data); 
     MaxOutValues[i] = max; 
     MinOutValues[i] = min; 
     MaxOutTimes[i] = maxtime; 
     MinOutTimes[i] = mintime; 
   } 
 
   if ((StartTime==0) && (sim==0)) { 
     data = fopen("Check_MaxMins", "w"); 
     fprintf(data, "Check_MaxMins: Max and Min values for each of the N runs:\n"); 
   } 
   else 
     data = fopen("Check_MaxMins", "a"); 
    
   fprintf(data, "Run %d:\n", sim); 
   for (i=0;i<NumEq; i++)  
     fprintf(data, "\t%s: %lg (@t=%lg)\t%lg (@t=%lg)\n", Eqns[i].name, MaxEqValues[i], MaxEqTimes[i], MinEqValues[i],  MinEqTimes[i]); 
   for (i=0;i<NumCalcOut;i++)  
     fprintf(data, "\t%s: %lg (@t=%lg)\t%lg (@t=%lg)\n", CalcOut[i].name, MaxOutValues[i], MaxOutTimes[i], MinOutValues[i], MinOutTimes[i]); 
   fprintf(data, "\n"); 
   fflush(data); 
   fclose(data); 
 
 } 
 
/* Write values in array TopValues--containing max value of each 
   equation variable, to be used to set yranges of graphs--to file 
   Check_TopValues */ 
void WriteTopValues(void) 
{ 
  FILE *data; 
  int i; 
 
  data = fopen("Check_TopValues", "w"); 
  for (i=0;i<NumEq; i++) { 
    /*    printf("NumEq[%d]: %lg\n", i, TopValues[i]); */ 
    fprintf(data, "%lg\n", TopValues[i]); 
  } 
   fflush(data); 
   fclose(data); 
} 
 
 
/*Utility Subroutines */ 
 
void Error(char* msg) 
{ 
  printf("ERROR:%s\n",msg); 
  exit(1); 
} 
/* This subroutine provides an alternate route for error conditions when 
we do not wish to abort the entire analysis, but rather just want to 
terminate the current simulation (as would be the case if we have stepsize 
underflow, a time overrun; etc . Instead of terminating the entire program 
it sends an error message to the files Check_Ranks, Check_PRCC, Log_File, 
Check_Uncert, and Check_Param, sets the PRCC to terminate at the current 
time point if time-dependent PRCCs are being done and sets an indicator to  
terminate the current simulation (SimDone=1)  It also resets the time to  
the endpoint of the integration, in order to avoid continuous looping 
without reaching the check point for SimDone =1 */ 
 
void ErrorSim2(char* msg) 
{ 
  printf("Simulation error:%s\n",msg); 
  fprintf(warn1,"Simulation error:%s\n",msg); 
  if (TimeAct<PRCCOffTime) PRCCOffTime=TimeAct;  
  SimDone=1;           /* Terminate the current simulation */ 
} 
 
 
int RandomInteger(int low, int high)      
{ 
  int k; 
  double d; 
  d = (double) rand()/((double) RAND_MAX + 1); 
  k = (int) (d*(high - low + 1)); 
  return(low + k); 
} 
 
 
/* This subroutine cleans up the input data and writes it to a             */ 
/* standardized format.  If typeclean==0, it takes the string, deletes all */ 
/* non-printing characters, but keeps blanks, converts multiple blanks to  */ 
/* a single blank. If typeclean==1 it deletes all non-printing characters, */ 
/* and deletes all blanks.                                                 */ 
          
/*                          Comment Lines 
*/ 
/* #1:   Loop through entire string character by character '\0' at end   */ 
/* #2:   If current character is a graphing character, write it to line2 */ 
/*       Flag that we've found a nonblank character.(first=1             */ 
/*       Move indices for each line forward by 1                         */ 
          
/*       If current character is nonvisible:                             */ 
/* #3    If we're not deleting all blanks, and we havent found           */ 
/*       the next blank, see if this character is.If it is put it in the */ 
/*       string Skip all consecutive blanks until we find a visible      */ 
/*       character again.                                                */ 
/* #4    Otherwise, we are deleting all nonvisible characters (including */ 
/*       blanks) , so move to next character                             */ 
/* #5    Replace original string (line) with edited copy (line2)         */ 
     
 
 
void StringClean(int typeclean) 
{ 
   int   i=0,j=0,k,first=0; 
   char  line2[MAXLINE]; 
 
   for (k=0;k<MAXLINE;k++) line2[k]=0; 
   while(line[j]!='\0')       /* #1 */ 
   { 
     if (isgraph(line[j]))    /* #2 */ 
      { 
         line2[i]=line[j]; 
         first=1; 
         i++; 
         j++; 
      } 
      else  if (typeclean!=1)   /* #3 */ 
      { 
         if (first==1) 
         { 
            if (isspace(line[j])) line2[i]=' '; 
            i++; 
            first=0; 
         } 
         j++; 
      } 
     else j++;                  /* #4 */ 
   } 
   for (k=0;k<MAXLINE;k++) line[k]=0; 
   strcpy(line,line2); 
} 
 
 
 
 
/* This subroutine takes a a matrix (matrix00) and inverts it using an LU Decomposition followed by    */ 
/* a back substitution. The inverted matrix is stored in inverse. During the inversion, matrix00 is    */ 
/* destroyed. We make a backup copy up front and then restore it at the end of the subroutine          */ 
/* Note that for this subroutine, we use a pointer to the first row of the matrix and then use         */ 
/* pointer arithmetic instead of passing the array directly. We do this to make the process faster, to */  
/* avoid issues between various C compilers with respect to passing of two dimensional arrays.         */ 
/* Furthermore the matrix used for the backup copy is declared as a two dimensional array instead of   */  
/* being allocated dynamically. This is to avoid repetitive dynamic allocation which may result in     */ 
/* memory leaks when a great number of allocations and deallocations are done.                         */ 
 
 
 
void InvertMatrix(long double *matrix,long double *inverse,int dim) 
{ 
   long double parity=1.0; 
   int *index; 
   long double sum=0.0,temp; 
   int i,j,k; 
   long double *col; 
   long double matrixcopy[MAXVARPARS][MAXVARPARS];     /* avoiding repetitive dynamic allocation, change to use code elsewhere */  
   long double determinant=0.0; 
 
   index = (int *) GetBlock(dim*sizeof(int)); 
   col = (long double *) GetBlock(dim*sizeof(long double)); 
 
   MatCopy((long double *) matrixcopy[0], (long double *)matrix,dim,dim);  /* make a copy of the matrix to be inverted */ 
   memset(index,0,dim*sizeof(int)); 
 
   LUDecomp((long double *) matrix,dim, index, &parity);        /* Invert matrix00, destroying it, inverse written to inverse */ 
   for (j=0;j<dim;j++) parity*=(*(matrix+dim*j+j)); 
   if (parity==0.0) Invert=0; 
   for(j=0;j<dim;j++)  
   { 
      for(i=0;i<dim;i++) col[i]=0.0; 
      col[j]=1.0; 
      LUBackSub((long double*) matrix,dim,(int *) index,(long double*)col); 
      for(i=0;i<dim;i++) (*(inverse+dim*i+j))=col[i]; 
   } 
     
   MatCopy((long double*) matrix,(long double *) matrixcopy[0],dim,dim);  /* Copy the backup of the matrix00 into the destroyed 
matrix */ 
   free(col); 
   free(index); 
   memset(matrixcopy,0,(MAXVARPARS)*(MAXVARPARS)*sizeof(long double));      /* reinitialize matrixcopy*/   
} 
 
 
/*   This subroutine prints out a copy of the argument matrix , of dimensions n by m to the argument file.   */ 
/*   given by *output. Again the matrix is passed using a pointer to avoid array issues and to speed calculation*/ 
 
 
void PrintMatrix(FILE *output,long double *matrix, int n,int m) 
{ 
   int i,j; 
 
   for (i=0;i<n;i++)  
   { 
      fprintf(output,"\n"); 
      for (j=0;j<m;j++)  
      {   
      if ((0.00000<=(*(matrix+m*i+j)))&&(*(matrix+m*i+j)<10.0))fprintf(output," %f  ",*(matrix+m*i+j)); 
      else fprintf(output,"%f  ",*(matrix+m*i+j)); 
      } 
   } 
   fprintf(output,"\n\n\n"); 
} 
 
 
 
void MatCopy( long double *target, long double *source, int n, int m) 
{ 
   int i,j; 
   
   for(i=0;i<n;i++) for(j=0;j<m;j++) *target++=*source++; 
} 
 
/* matrix 1 is nxm matrix 2 is mxp product is nxp */ 
 
void MatMult( long double *matrix1, long double *matrix2, long double *product, int n, int m,int p) 
{ 
   int i,j,k; 
   long double sum; 
   
   for (i=0;i<n;i++)  
   {  
      for(j=0;j<p;j++) 
      { 
         sum = 0.0 ; 
         for(k=0;k<m;k++) sum += (*(matrix1+n*i+k))*(*(matrix2+m*k+j)); 
         *(product+n*i+j) = sum; 
      } 
   }    
} 
 
 
 
 
/* d is output as +1 or -1 depending on whether the number of  
row interchanges was negative or positive This implementation  
comes almost verbatim from Numerical Recipes.*/ 
 
void LUDecomp(long double *matrix00, int n, int *indx, long double *d) 
{ 
int i,j,k,imax; 
long double big, dum, sum, temp; 
long double *ImpScaling;   
 
   ImpScaling = (long double* ) GetBlock(n*sizeof(long double));  
   *d = 1.0; 
 
/* To begin, find the largest element in each row,store in the vector */ 
/*   ImpScaling this will be used for scaling later. If the biggest */ 
/*  element in any row is zero, then the matrix is singular         */ 
    
    for (i=0;i<n;i++) 
    {  
      big = 0.0; 
      for (j=0;j<n;j++)  
      { 
         temp=fabs(*(matrix00+n*i+j)); 
         if (temp > big) big = temp; 
      } 
      if ((big ==0.0) && (DoLocal==1) ) 
      {           
          MakeLog(); 
          Error ("Singular matrix in routine LUDecomp"); 
      } 
      if ((big==0.0) && (DoGlobal==1)) 
      { 
         Invert=0; 
         return; 
      } 
 
 
      ImpScaling[i] = (1.0/big);  
  
   } 
 
/* Use Crout's method to solve for the alphas and betas  */ 
 
   for (j=0;j<n;j++)  
   { 
      for (i=0;i<j;i++)        /* first solve betaij= aij-sum{k=1,j-1} (alphaik)(betakj) */ 
      { 
          sum = (*(matrix00+n*i+j)); 
          for (k=0;k<i;k++) sum -=(*(matrix00+n*i+k))* (*(matrix00+n*k+j)); 
          (*(matrix00+n*i+j)) = sum; 
      } 
      big = 0.0; 
      for (i=j;i<n;i++)     
      { 
          sum = (*(matrix00+n*i+j)) ; 
          for (k=0;k<j;k++) sum -=(*(matrix00+n*i+k))*(*(matrix00+n*k+j)); 
          (*(matrix00+n*i+j)) = sum; 
          if (sum>= 0.0) dum = ImpScaling[i]*(sum); 
          else dum=(-1)*ImpScaling[i]*sum;  
          if (dum >= big ) 
          { 
             big = dum; 
             imax = i; 
          }  
      }  
      if ( j!=imax)     /* need to interchange rows*/ 
      { 
         for( k=0;k<n;k++) 
         {   
            dum = (*(matrix00+n*imax+k)); 
            (*(matrix00+n*imax+k)) = (*(matrix00+n*j+k)); 
            (*(matrix00+n*j+k)) = dum; 
         } 
         *d = -(*d); 
         ImpScaling[imax] = ImpScaling[j]; 
      } 
      *(indx+j) = imax; 
      /* if ((*(matrix00+n*j+j)) != 0.0 ); else  (*(matrix00+n*j+j)) = TINY;  pivotvalue is numerically zero */ 
      if (j!=(n-1)) 
      { 
         dum = 1.0/(*(matrix00+n*j+j)); 
         for(i=j+1;i<n;i++)(*(matrix00+n*i+j)) *=dum; 
      } 
   } 
   free(ImpScaling); 
} 
 
 
/* This subroutine does the backsubstitution procedure. Note that this implementation of the LUBackSub differs from Numerical */ 
/* Recipes in that the A matrix is not specifically passed in but is passed in  */ 
/* as a pointer to the first element - this is to avoid two-dimensional array issues */ 
/* in gcc  */ 
 
void LUBackSub(long double* matrix00, int n, int *indx, long double *b) 
{ 
   int i,k= 0,ip,j; 
   long double sum; 
 
   for (i=0;i<n;i++)  
   { 
      ip = (*(indx+i)); 
      sum = (*(b+ip)); 
      *(b+ip) = (*(b+i)); 
      if (k>-1) for (j=k;j<=i-1;j++) sum -=(*(matrix00+n*i+j))*(*(b+j)); 
      else if (sum) k = i;   /*if a nonzero element was encountered movek forward */ 
       *(b+i) = sum; 
   } 
   for ( i=n-1;i>=0;i--)  
   { 
      sum = (*(b+i)); 
      for (j=i+1;j< n;j++) sum -= (*(matrix00+n*i+j))*(*(b+j)); 
      *(b+i) = sum/(*(matrix00+n*i+i)); 
   } 
}  
 
 
 
static void *GetBlock(size_t nbytes) 
{ 
        void *result; 
 
        result = (void *)malloc(nbytes); 
        if(result == NULL) {printf("No memory available"); exit(1);} 
        return(result); 
}   
 
 
void OutputRunParams(int sim) 
{ 
    FILE *parsample; 
    int i;  
    sprintf(line,"ParSample%d",sim); 
    parsample = fopen(line,"w"); 
    for(i=0;i<NUMEQ;i++)fprintf(parsample,"%lg\n",TimeMatrix[LastStep-1][i]); 
    for(i=0;i<NUMCALCOUT;i++)fprintf(parsample,"%lg\n",CalcTime[LastStep-1][i]); 
    for(i=NUMEQ;i<MAXPARS;i++)fprintf(parsample,"%lg\n",LHSMatrix[i][sim].value); 
    fclose(parsample); 
} 
 
 
 
 
void GetPreviousParams(int sim) 
{ 
   int i; 
   FILE *paramsample; 
    
 
   sprintf(line,"ParSample%d",sim); 
   paramsample = fopen(line,"r"); 
   for(i=0;i<NUMEQ;i++) 
   { 
       fgets(line,MAXLINE,paramsample); 
       sscanf(line,"%lg\n",&LHSMatrix[i][sim].value); 
   } 
   for(i=0;i<NUMCALCOUT;i++) 
   { 
       fgets(line,MAXLINE,paramsample); 
       sscanf(line,"%lg\n",&CalcTime[0][i]);  
   } 
   for(i=NUMEQ;i<MAXPARS;i++) 
   { 
      fgets(line,MAXLINE,paramsample); 
      sscanf(line,"%lg\n",&LHSMatrix[i][sim].value); 
   } 
   fclose(paramsample); 
   PrintLabelSample(LabelPar); 
} 
 
 
 
/* This subroutine computes a forward difference approximation to the Jacobian matrix (df/dy) where f(y,t) = dy/dt        */ 
/* The forward finite difference approximation is based on the formula: f'(x_0) = f(x_0+h) - f(x_0) / h )                 */  
/* This method is inherently unstable, with a relation ship of the truncation error scaling with h^3 and the              */ 
/* roundoff error scaling with 1/h. Thus decreasing truncation error increases roundoff error and vice versa. To avoid    */ 
/* this issue we use a trick where we set the independent variable equal to a temp variable, and then set the independent */ 
/* variable equal to the temp value + h and then h to the value of the independent variable minus the temp varialbe. This */ 
/* makes the internal respresentation of the independent variable and the next point x_0+h differ by a constant           */ 
/* See the LHS help file for further explanation of the method                                                            */ 
 
 
void FiniteDiffJacob (int n, long double *indvar,long double *fvec, long double *df00) 
{ 
   int i, j; 
   long double h, temp, *fpdelta; 
 
 
   fpdelta = GetBlock(n*sizeof(long double)); 
 
   for(j=0;j<n;j++) 
   { 
      temp=(*(indvar+j));     /*preserve the original value of the independent variable */ 
      h = EPSJ*fabs(temp); 
      if (h==0.0) h=EPSJ; 
      *(indvar+j)=temp+h;       /* Trick to reduce finite precision error */ 
      h = (*(indvar+j)) - temp; 
      GetDerivatives((long double *) indvar, (long double *)fpdelta, TimeAct);   /* get f(y+h) */ 
      *(indvar+j) = temp;               				         /* restore original y value */ 
      for (i=0;i<n;i++) (*(df00+n*i+j)) = ((*(fpdelta+i)) -(*(fvec+i)))/h; 
   } 
   free(fpdelta); 
} 
 
 
 
/* This subroutine calculates df/dt where dy/dt=f(t,y)  and t occurs explicitly on the right-hand */ 
/* side of the equation */ 
 
void GetTimeDeriv(int n, long double *indvar,long double *fvec, long double *df00) 
{  
   int i, j; 
   long double h, temp, *fpdelta,newtime; 
          
       
   fpdelta = GetBlock(n*sizeof(long double)); 
 
   for(j=0;j<n;j++) 
   { 
      temp=TimeAct;     /*preserve the original value of the independent variable */ 
      h = EPSJ*fabs(temp); 
      if (h==0.0) h=EPSJ; 
      newtime= temp+h;       /* Trick to reduce finite precision error */ 
      h = newtime - temp; 
      GetDerivatives((long double *) indvar, (long double *)fpdelta, newtime);   /* get f(y+h) */ 
      for (i=0;i<NUMEQ;i++) (*(df00+n*i+j)) = ((*(fpdelta+i))-(*(fvec+i)))/h; 
   } 
   free(fpdelta); 
} 
 
 
 
 
 
 
/* 
         ModifiedMid((long double *) ysav,(long double*)timederiv,NUMEQ,TimeAct,CurrentIncrement,nseq[k],(long double *)yseq, 
         (long double *) derivArray); 
    
    
This subroutine performs one step of the modified midpoint method, which advances a vector of 
dependent variables y(t) from a point t to a point t+H by a sequence of n substeps each of size 
h = H/n See Numerical recipes in C ppg 722 for further details */ 
 
void ModifiedMid( long double *yin, int numeq,long double start,long double bigstepsize, int numsubsteps,long double *yout,long double *deriv) 
{ 
   int i,j,n; 
   long double timelocal,swap, h2, h; 
   long double ym[NUMEQ], yn[NUMEQ]; 
 
 
   h = bigstepsize/numsubsteps; 
   h2 = 2.0*h; 
   timelocal = start+h;  
   for (j=0;j<NUMEQ;j++)        /* first step*/ 
   { 
      ym[j] = (*(yin+j)); 
      yn[j] = (*(yin+j))+h*(*(deriv+j)); 
   } 
       
 
   GetDerivatives((long double * )yn,(long double *) yout,timelocal);  /* Make sure this is compatible */ 
   
 
   for (n=1;n<numsubsteps;n++)    /* steps 2 to numsubsteps - the array indices ARE right*/ 
   { 
      for(j=0;j<numeq;j++) 
      { 
         swap=ym[j]+h2*(*(yout+j)); 
         ym[j]=yn[j]; 
         yn[j]=swap; 
      } 
      timelocal+=h; 
      GetDerivatives((long double *) yn,(long double *) yout,timelocal); 
   } 
   for(i=0;i<numeq;i++) (*(yout+i))= 0.5*(ym[i]+yn[i]+h*(*(yout+i)));    /* last step */ 
 
} 
 
 
 
/* This subroutine uses polynomial extrapolation to to evaluate  
numeq equations at t = 0 by fitting a polynomial to a sequence of  
estimates with progressively smller values t = xest, and correpsonding  
function vectors yest[1,... numeq]. This call is number numpts in the  
sequence of calls. Extrapolated function values are output as  
yextrap[1,... numeq] and their estimated error is output as ydiff[1,...numeq]. 
 
PolyExtrapolate(k,xest,(long double *) yseq,(long double *)WorkValArray,(long double *) yerr,NUMEQ); 
 
 
*/ 
 
 
void PolyExtrapolate( int numpts,long double xest,long double *yest,long double *yextrap,long double *ydiff,int numeq) 
{ 
   int l,j; 
   long double q, f2, f1, delta; 
   long double c[NUMEQ];  
 
   x[numpts-1]=xest;    /* Save the current value of the independent variable */ 
   for(j=0;j<NUMEQ;j++) (*(ydiff+j))= (*(yextrap+j)) = (*(yest+j));  /*initialize error and extrapolated value estimates to the last estimate */ 
   if (numpts==1) for (j=0;j<NUMEQ;j++) d[j][0]=(*(yest+j));    /* store first estimate (the ds) in the first column */ 
   else 
   { 
      for (j=0;j<NUMEQ;j++) c[j]=(*(yest+j)); 
      for (l=1;l<numpts;l++) 
      { 
         delta = 1.0/(x[numpts-l]-xest-1);  /* difference between the current independent variable value and the the current-lthvalue*/ 
         f1 = xest*delta; 
         f2 = x[numpts-l-1]*delta; 
         for (j=0;j<NUMEQ;j++) 
         { 
            q=d[j][l-1]; 
            d[j][l-1]=(*(ydiff+j)); 
            delta=c[j]-q; 
            (*(ydiff+j)) = f1*delta; 
            c[j] = f2*delta; 
            *(yextrap+j) += (*(ydiff+j)); 
 
         } 
      } 
      for (j=0;j<NUMEQ;j++) d[j][numpts-1]=(*(ydiff+j)); 
    } 
 
} 
  
 
/*This subroutine performs one step of a Burlisch - Stoer method (either a nonstiff method (Method ==5)  or a stiff  
method (Method==6) ) for integration of differential equations. It is called by 
ODEStepper, which controls the stepsize and starts and stops the integration. BSStep takes values at beginning of the interval 
(stored in WorkValArray), together with the derivative values at that point, (stored in derivArray) and advances the 
estimate by one step. It does so by estimating the value of the functions at the time point t+H by a sequence of n  
substeps of size h via a modified midpoint method. The interval is first spanned by a few points, and then is repeated  
with finer and finer coverage of the interval, these results are then used to extrapolate (via polynomial exptrapolation)  
to the answer that would have been gotten using an infinite number of substeps (i.e if h had been 0). 
 
NOTE:  For the stiff solver, KMAXX = 7, IMAXX = 8, For the nonstiff solver, KMAXX = 8, IMAXX = 9   - this is written to  
temp.c by the setup program as #define KMAXX  7 (or 8)  
                               #define IMAXX KMAXX+1    
 
*/ 
   
void BSStep(long double *ScaleArray, long double *DidIncrement, long double *NextIncrement, long double *WorkValArray, 
long double *derivArray, long double *timederiv,long double eps) 
{ 
   long double  eps1, errmax, fact, red , scale , workdone, workmin, xest; 
   long double  err[KMAXX]; 
   long double  yerr[NUMEQ]; 
   long double  ysav[NUMEQ]; 
   long double  yseq[NUMEQ]; 
   int          reduct,exitflag=0; 
   int          i,j,k,l,m; 
   long double CurrentIncrement; 
 
   static int          first =1, kmax, colopt,nvold=-1; 
   static long double  epsold = -1.0, TimeNew;    
   static long double  a[IMAXX+1]; 
   static long double  alf[KMAXX+1][KMAXX+1]; 
   static int          nseq[IMAXX+1]; 
   long double test; 
 
   if (Method==5) for (i=0;i<IMAXX+1;i++) nseq[i]=2*i;  /* for nonstiff solver nseq = {0,2,4,6 ,8 ,10,12,14,16,18} */ 
   if (Method==6)                                    /* for stiff solver    nseq = {0,2,6,10,14,22,34,50,70}    */ 
   {  
      nseq[0] = 0; 
      nseq[1] = 2; 
      nseq[2] = 6; 
      nseq[3] = 10; 
      nseq[4] = 14; 
      nseq[5] = 22; 
      nseq[6] = 34; 
      nseq[7] = 50; 
      nseq[8] = 70; 
   }    
 
   if ((eps != epsold) || ((Method==6) && (NumEq!=nvold)))  
   { 
      *NextIncrement = TimeNew = -1.0e29;  /* New tolerance, so reinitialize to impossible values */ 
      eps1 = SAFE1*eps; 
      a[1] = nseq[1]+1;    /* Check this - this is difficult, nseq is 0 based in NR as well a, a is never initialized for the 0 element*/ 
      for(k=1;k<=KMAXX;k++) a[k+1] = a[k] + nseq[k+1];  /* a runs from 0 to IMAXX+1-1 = IMAXX = KMAXX+1, so loop runs, k=1, KMAXX, initialize a[2] to a[KMAXX+1] */ 
      for(j=2;j<=KMAXX; j++) 
      { 
         for ( k=1;k<j;k++) alf[k][j] = pow(eps1,((a[k+1]-a[j+1]) /((a[j+1] - a[1]+1.0) * (2*k+1))));   /* again alf is initialize from row 1, column 2 ? */ 
      } 
      epsold = eps; 
      if (Method==6)    /* if we use the stiff method, add cost of evaluating the wJacobian to the overall work calculation */ 
      { 
         nvold = NumEq; 
         a[1] +=NumEq;   
         for (k=1;k<=KMAXX; k++) a[k+1] = a[k] + nseq[k+1]; 
      } 
      for(colopt = 2;colopt<KMAXX; colopt++) if (a[colopt+1]> a[colopt]*alf[colopt-1][colopt]) break;    /* Check eqn 16.4.13 */ 
      kmax = colopt;              /* the maximum number of columns to try is equal to the estimatated  optimal column - these are done based on 0 based arrays used as 1 based arrays , may need to scale kmax*/ 
   } 
    
   CurrentIncrement = TimeIncrement;  
   reduct = 0; 
 
   for (j=0;j<NUMEQ;j++) ysav[j] =  (*(WorkValArray+j));  /* both are zero based here */ 
   if (Method == 6)  
   { 
      FiniteDiffJacob(NumEq,(long double *) WorkValArray,(long double *) derivArray,(long double *)Jacobian[0]); 
      GetTimeDeriv(1, (long double *) WorkValArray, (long double *) derivArray, (long double *) timederiv); 
   }      
   if ((TimeAct !=TimeNew) || (CurrentIncrement!=(*(NextIncrement))))  
   { 
      first=1; 
      colopt=kmax; 
   } 
 
   for (;;) 
   { 
      for(k=1;k<=kmax;k++) 
      { 
         TimeNew = TimeAct+ CurrentIncrement; 
         if (TimeNew+TINY == (TimeAct)+TINY) ErrorSim2("Step size underflow in Burlisch Stoer stepper algorithm"); 
         if (Method ==5) ModifiedMid((long double *)ysav,NUMEQ,TimeAct,CurrentIncrement,nseq[k],(long double *)yseq,(long 
double *) derivArray); 
         if (Method ==6) SemiImpModMid ( (long double *)ysav,(long double*) timederiv,NUMEQ,TimeAct,CurrentIncrement,nseq[k], (long double 
*)yseq, (long double *) derivArray); 
 
 
         xest = (CurrentIncrement/nseq[k]) *(CurrentIncrement/nseq[k]);  /* This is why use the 1 based loop, both nseq (this and NR have a 0 based nseq array  */ 
         PolyExtrapolate(k,xest,(long double *) yseq,(long double *) WorkValArray,(long double *) yerr,NUMEQ); 
         for (i=0;i<NUMEQ;i++) if (*(WorkValArray+i)<-TINY) 
         { 
             (*(WorkValArray+i))=0.00000; 
             if (Warn==0) fprintf(warn1, "\tWarning: some equation values were set to zero \n \twhen solver produced negative outcome values\n"); 
              Warn=1;  
         } 
         if(k!=1) 
         { 
            errmax = TINY; 
            for (i=0;i<NUMEQ;i++) if (errmax <fabs(yerr[i]/ScaleArray[i]))errmax= fabs(yerr[i]/ScaleArray[i]); 
            errmax = errmax/eps; 
            m = k-1;  /* this is the same - make translation where its used below */ 
            err[m-1] = pow(errmax/SAFE1,1.0/(2*m+1)); 
         } 
         if ((k!=1) && ( k>=colopt-1 || first )) 
         { 
            if (errmax <1.0) 
            { 
               exitflag = 1; 
               break; 
            } 
            if ( (k==kmax )||( k==colopt+1)) 
            {    
               red = SAFE2/err[m-1]; 
               break;   
            } 
            else if ( (k==colopt) && (alf[colopt-1][colopt] <err[m-1])) 
            { 
               red = 1.0/err[m-1];   
               break; 
            } 
            else if ((colopt == kmax) && (alf[m][kmax-1] <err[m-1]))  /* alf is used as a one based array*/ 
            { 
               red = alf[m][kmax-1]* SAFE2/err[m-1]; 
               break; 
            } 
            else if (alf[m][colopt] < err[m-1]) 
            { 
               red = alf[m][colopt-1]/err[m-1]; 
               break; 
            } 
             
         }    
      } 
      if (exitflag) break; 
      if (red>REDMIN) red=REDMIN;   
      if (red<REDMAX) red =REDMAX; 
      CurrentIncrement *= red; 
      reduct = 1; 
   } 
   TimeAct = TimeNew; 
   *DidIncrement = CurrentIncrement; 
   first =0; 
   workmin = 1.0e35; 
   for (l = 0;l<m;l++) 
   { 
      ((err[l]>SCALMX) ? (fact = err[l]): (fact =SCALMX)); 
      workdone = fact * a[l+1]; 
      if (workdone < workmin) 
      {    
         scale = fact ; 
         workmin = workdone; 
         colopt = l + 1; 
      } 
   } 
   *NextIncrement = CurrentIncrement/scale; 
   if ((colopt >= k) && ((colopt !=kmax) && (!reduct))) 
   { 
      test = scale/alf[colopt -1][colopt]; 
      ((test> SCALMX)? (fact=test): (fact=SCALMX)); 
      if ( a[colopt+1]*fact <= workmin) 
      {      
         *NextIncrement = CurrentIncrement/fact; 
         colopt ++; 
      } 
   } 
   memset(yseq,0,(NUMEQ)*sizeof(long double)); 
   memset(ysav,0,(NUMEQ)*sizeof(long double)); 
   memset(yerr,0,(NUMEQ)*sizeof(long double)); 
   memset(x,0,KMAXX*sizeof(long double)); 
   memset(err,0,KMAXX*sizeof(long double)); 
   memset(d,0,(NUMEQ)*KMAXX*sizeof(long double)); 
} 
 
 
/* This subroutine  performs one step of  a semi-implicit modified midpoint method which is used  
by the stiff version of the Burlisch - Stoer method */ 
 
void SemiImpModMid ( long double *yin,long double *timederiv,int numeq,long double start, long double bigstepsize,  
int numsubsteps, long double *yout,long double *deriv) 
{ 
 
   int i,j,nn; 
   int index[NUMEQ]; 
   long double d,h,x; 
   long double a[NUMEQ][NUMEQ]; 
   long double del[NUMEQ]; 
   long double ytemp[NUMEQ]; 
 
   h= bigstepsize/numsubsteps; 
   for (i=0;i<NUMEQ;i++) 
   { 
      for (j=0;j<NUMEQ;j++) a[i][j] = -h*(Jacobian[i][j]); 
      ++a[i][i]; 
   } 
   LUDecomp((long double *) a[0],NUMEQ,(int *)index,&d); 
   for (i=0;i<NUMEQ;i++) *(yout+i)= h*(*(deriv+i)+h*(*(timederiv+i))); 
   LUBackSub((long double *) a[0],NUMEQ,(int *) index, (long double *)yout); 
   for (i=0;i<NUMEQ;i++) *(ytemp+i) = (*(yin+i)) + (*(del+i) =(*(yout+i))); 
   x = start +h; 
   GetDerivatives((long double *) ytemp, (long double *) yout,x); 
   for (nn=1;nn<numsubsteps;nn++) 
   { 
      for(i=0;i<NUMEQ;i++) *(yout+i) = h*(*(yout+i))-(*(del+i)); 
      LUBackSub((long double *) a[0], NUMEQ,(int *) index, (long double *)yout); 
      for (i=0;i<NUMEQ;i++) (*(ytemp+i)) += (*(del+i) += 2.0*(*(yout+i))); 
      x +=h; 
      GetDerivatives((long double *) ytemp, (long double *) yout, x); 
   } 
   for (i=0;i<NUMEQ;i++) (*(yout+i)) = h*(*(yout+i))- (*(del+i)); 
   LUBackSub((long double *) a[0], NUMEQ, (int *) index, (long double *)yout); 
   for (i=0;i<NUMEQ;i++) (*(yout+i)) += (*(ytemp+i)); 
   memset(ytemp,0,(NUMEQ)*sizeof(long double)); 
   memset(del,0,(NUMEQ)*sizeof(long double)); 
   memset(a,0,(NUMEQ)*(NUMEQ)*sizeof(long double)); 
   memset(index,0,(NUMEQ)*sizeof(int)); 
} 
 
 
void PrintLabelSample(int element) 
{ 
   int i; 
   FILE *labelsample; 
 
 
   labelsample = fopen("LabelSample","w"); 
   for (i=0;i<NumSim;i++)  
   { 
       fprintf(labelsample,"%lg\n",LHSMatrix[element][i].value); 
   } 
   fclose(labelsample); 
} 
 
 
void SaveOutComeData(int ptssaved,int sim) 
{ 
   int i; 
   int j=0; 
      
   if (DoGlobal==1) 
   { 
       if (OutNum<NumEq) OutComeData[sim][ptssaved].outcome[0]=TimeMatrix[ptssaved][OutNum]; 
       else OutComeData[sim][ptssaved].outcome[0]=CalcTime[ptssaved][OutNum-NumEq]; 
       OutComeData[sim][ptssaved].time=TimeAct; 
   }   
   if (DoLocal==1) 
   { 
        if (OutNum<NumEq) OutComeData[sim][1].outcome[0]=TimeMatrix[ptssaved][OutNum]; 
        else OutComeData[sim][1].outcome[0] = CalcTime[ptssaved][OutNum-NumEq]; 
        OutComeData[sim][1].time = TimeAct;      
        fprintf(warn1,"The outcome variable saved for simulation %d was %lf at time %lf\n",sim,OutComeData[sim][1].outcome[0],OutComeData[sim][1].time);   
        if((DoLocal==1) && (sim==0))  
        { 
            fprintf(check[5], "\n \nThe outcome variable was %s at time %lf\n ", OutName,OutDay); 
            fprintf(check[4], "\nThe outcome variable was %s at time %lf \n \n", OutName,OutDay); 
            fprintf(check[1], "\n \nThe outcome variable was %s at time %lf \n ",OutName, OutDay); 
        } 
   } 
         
} 
 
/* C3: If the difference between the max and the min outcome times that was used is greater than the maximum  
difference allowed, AND this difference is not greater than the current time (which would mean that some outcomes 
were still initialized to zero because the integration terminated early ) , and we're far enough along in the  
integration to worrry about this, then print out the actual outcome values and times that were used with a warning message */ 
/* C4: If the difference betweent he max and the min outcome times that was used is greater than the current time -1 and 
we're far enough along in the integration to care, we probably have time points that were not reached due the early 
termination of one of the simulations. Therefore,  we should just terminate the 
PRCC calculations at this time point as well (none will exist beyond this point because the outcome value was not  
calculated beyond this point for the simulation that terminated early) */ 
 
              
void PlaceNextOutCome(Call) 
{ 
      int k; 
      double min, max; 
    
      CurrOutTime=0.0; 
      for(k=0;k<NumSim;k++)   /* copy the outcomes into the lhsmatrix and pick out the minimum outcome time that was used */ 
      { 
          if (k==0) min=max=OutComeData[k][Call].time; 
          else 
          { 
             if (OutComeData[k][Call].time<min) min=OutComeData[k][Call].time; 
             if (OutComeData[k][Call].time>max) max=OutComeData[k][Call].time; 
          } 
          LHSMatrix[NumPar][k].value=OutComeData[k][Call].outcome[0]; 
      } 
      CurrOutTime=max;                  /* Set the current outcome time to the minimum that was used */ 
      printf("\n\nPRCC for outcome %s at time %lf being calculated\n",OutName,CurrOutTime); 
      if ( ((max-min>MAXPRCCDIFF) && (max-min<CurrOutTime-1.0)) && (CurrOutTime>1.0))   /* 3 */ 
      { 
        printf("Warning: Maximum difference of %lf at PRCC outcome time %lf\n",max-min,CurrOutTime); 
        printf("Reduce maximum stepsize to make sampled outcome time points closer\n"); 
        for (k=0;k<NumSim;k++) printf("Outcomes Placed: %lf %lf\n",LHSMatrix[NumPar][k].value,OutComeData[k][Call].time); 
      } 
      if ((max-min>=CurrOutTime-1.0) && (CurrOutTime>1.0)) DonePRCC=1; 
      if (DoGlobal==1)  
      { 
         if (Call==1) fprintf(check[5], "\n \nThe outcome variable was %s\n",OutName,CurrOutTime); 
         if (Call==1) fprintf(check[4], "\nThe outcome variable was %s\n\n",OutName,CurrOutTime); 
         if ((Call==1) && (KEEPTDRANKS==1)) fprintf(check[1], "\n \nThe outcome variable was %s at time %lf\n",OutName, CurrOutTime); 
         if ((Call==1) && (KEEPTDINVERT==1)) fprintf(check[3],"\n\nThe outcome variable was %s at time %lf\n",OutName, CurrOutTime); 
      } 
      if (CurrOutTime>=PRCCOffTime) DonePRCC=1; 
} 
      
 
 
 
/*  This subroutine is called if the matrix inversion revealed a singular CijMatrix. In this case we need to                 */ 
/*  see if the matrix was singular because one RankMatrix rows was the same or opposite the bottom row (row NumVarPar)       */ 
/*  (opposite meaning each element is equal to NumSim+1-the corresponding element in the bottom row.                         */ 
/*  If this is the case, we need to signal that the matrix was singular because one of the parameters                        */ 
/*  strongly correlated with the outcome values and set the PRCC for that parameter equal to either 1.0 or -1.0              */ 
 
/*                                       LINE COMMENTS                                                                   */ 
/* #1:   If we already checked the parity and it's equal to 1, rows are different                                        */ 
/* #2:   If we haven't checked parity or we know parity is -1, check if 
they are opposite                                */ 
/* #3:   If they're opposite, if it's the first row element, set parity to -1, otherwise do nothing (same is still = 1). */ 
/* #4:   Otherwise, if they're also not opposite, the rows are different (set same to 0)                                 */ 
/* #5:   Go to the next row element                                                                                      */   
/* #6:   If we still haven't found an element that differs, and there's another element in the row, compare the next two.*/ 
 
 
void CheckStrongEffect(void) 
{ 
     
   int i,j,same; 
 
   printf("Singular matrix at time %lf.\nChecking for strong effect parameter....\n",CurrOutTime); 
   StrongEffect=0; 
   StrongPar=NumVarPar; 
   for (j=0;j<NumVarPar;j++) 
   { 
      same=1;  
      i=0; 
      do 
      {    
         if (RankMatrix[j][i]!=RankMatrix[NumVarPar][i]) 
         { 
            if ((i>0) && (Parity==1)) same=0;   /* #1 */ 
            if ((i==0) || (Parity==-1))         /* #2 */ 
            { 
               if (RankMatrix[j][i]==NumSim+1-RankMatrix[NumVarPar][i]) if (i==0) Parity=-1;    /* #3 */     
               else same=0;                                                  /* #4 */ 
            } 
         }   
         i++;                                /* #5 */ 
      } 
      while ((same==1) && (i<NumSim));       /* #6 */    
      if (same==1) StrongPar=j; 
   } 
   if (StrongPar<NumVarPar) StrongEffect=1; 
   if (StrongEffect==1) printf("Strong effect parameter:%s\n",pArray[NewI[StrongPar]].name); 
   else printf("No strong effect parameter found\n"); 
} 
 
 
 
void ReadEquations(void) 
{ 
   int typeline=0;    /* Type of line read in from equation file: EqnName  CutoffValue (0) 
   or Equation Calculation (1) */ 
   int i,j,k; 
   int firstblank,firstequal;    /* used to mark first occurence of blank in current line */ 
   int Section=0; 
   FILE* feqns; 
   char stuff[MAXLINE]; 
           
   feqns =fopen(eqnfile,"r"); 
   if(feqns == NULL) Error("Cannot find equation file"); 
   NumEq=0; 
   NumCalcOut=0; 
   while (fgets(line,MAXLINE,feqns) !=NULL) 
   { 
      StringClean(typeline); 
      if(strlen(line)==0) continue; 
      if(strncmp(line,"Equations",9)==0) continue; 
      if (strncmp(line,"Where",5)==0) 
      { 
          Section=1; 
          continue; 
      } 
      if ((strncmp(line,"OutComes",8)==0) || (strncmp(line,"Outcomes",8)==0)) 
      { 
         Section=2; 
         continue; 
      } 
           switch(Section) 
      { 
         case 0 : 
         { 
            if(typeline==0) sscanf(line,"%s %lf\n",Eqns[NumEq].name,&Eqns[NumEq].top); 
            if(typeline==1) sscanf(line,"%s \n",stuff); 
            if(typeline==1) NumEq++; 
            typeline=1-typeline; 
            break; 
         } 
         case 1 : 
         case 2 : 
         { 
            firstequal=-1; 
            j=0; 
            do firstequal++ ; 
            while ((line[firstequal])!='=') ; 
            if (firstequal > MAXNM) Error("Maximum term name length exceeded"); 
            if (Section==2) for(k=0;k<firstequal;k++) if (isspace(line[k])==0) CalcOut[NumCalcOut].name[k] =line[k]; 
            if (Section==2) NumCalcOut++; 
            break; 
         } 
      } 
   } 
   if (OutNum<NumEq) sprintf(OutName,"%s",Eqns[OutNum].name); 
   else sprintf(OutName,"%s",CalcOut[OutNum-NumEq]); 
   printf ("The name of the outcome variable is :%s\n",OutName); 
} 
 
 
/* This is the end of the program! */ 
