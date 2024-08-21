// 
// This is the new prcc 
// Last updated 03/23/03
// 

#include <iostream>
#include <iomanip>   
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cstring> 

#include "nrutil.h"
#include "my_prcc2.h"

#define NUMPRCCOUTS 120 	/* Max number of PRCC outcomes */ 
#define EPSJ 1e-10 
#define TINY 1e-30
#define MAXPRCCDIFF 1.0

using namespace std;

void error(const char * ab, int tok)
{
  printf("Error: %s %d \n", ab, tok);
  exit(1); }

void error(string ab, int linen)
{ 
  cout << "Error: " << ab  << " " << linen << endl;
  exit(1); 
}

void error(string ab)
{ 
  cout << "Error: " << ab  << endl;
  exit(1); }

void warning(string ab)
{ 
  cout << "Warning: " << ab  << endl;
}

typedef struct 
{ 
  double outcome[NUMPRCCOUTS]; 
  double time; 
} timepoint; 
 
 
typedef struct 
{ 
  string name; 
  double min; 
  double peak; 
  double max; 
  char dType; 
} paramT; 
  
typedef struct 
{ 
  double value; 
  double rank; 
  int tieset; 
} matrixDB; 

int Invert=1;   
int StrongEffect =0;
char line[256];    

int myreturn = 0;

void PrintMatrix(double *matrix, int n,int m) 
{ 
  int i,j; 
  for (i=0;i<n;i++)  
    { 
      printf("\n"); 
      for (j=0;j<m;j++)  
	{   
	  if ((0.00000<=(*(matrix+m*i+j)))&&(*(matrix+m*i+j)<10.0)) printf(" %.4e ",*(matrix+m*i+j)); 
	  else printf("%.4e ",*(matrix+m*i+j)); 
	} 
    } 
  printf("\n"); 
} 
 
void MatCopy(double *target,   double *source, int n, int m) 
{ 
  int i,j; 
   
  for(i=0;i<n;i++) for(j=0;j<m;j++) *target++=*source++; 
} 
 
/* matrix 1 is nxm matrix 2 is mxp product is nxp */ 
 
void MatMult(   double *matrix1,   double *matrix2,   double *product, int n, int m,int p) 
{ 
  int i,j,k; 
  double sum; 
   
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
    

static void *GetBlock(size_t nbytes) 
{ 
  void *result; 
 
  result = (void *)malloc(nbytes); 
  if(result == NULL) {printf("No memory available"); exit(1);} 
  return(result); 
}   
 


/* d is output as +1 or -1 depending on whether the number of  
   row interchanges was negative or positive This implementation  
   comes almost verbatim from Numerical Recipes.*/ 
 
void LUDecomp(  double *matrix00, int n, int *indx,   double *d) 
{ 
  int i,j,k,imax; 
  double big, dum, sum, temp; 
  double *ImpScaling;   
 
  ImpScaling = (  double* ) GetBlock(n*sizeof(  double));  
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
      if ((big==0.0))// && (DoGlobal==1)) 
	{ 
	  //cout << " Singular Matrix !" << endl;
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
 
void LUBackSub(  double* matrix00, int n, int *indx,   double *b) 
{ 
  int i,k= 0,ip,j; 
  double sum; 
 
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
 

/* This subroutine takes a a matrix (matrix00) and inverts it using an LU Decomposition followed by    */ 
/* a back substitution. The inverted matrix is stored in inverse. During the inversion, matrix00 is    */ 
/* destroyed. We make a backup copy up front and then restore it at the end of the subroutine          */ 
/* Note that for this subroutine, we use a pointer to the first row of the matrix and then use         */ 
/* pointer arithmetic instead of passing the array directly. We do this to make the process faster, to */  
/* avoid issues between various C compilers with respect to passing of two dimensional arrays.         */ 
/* Furthermore the matrix used for the backup copy is declared as a two dimensional array instead of   */  
/* being allocated dynamically. This is to avoid repetitive dynamic allocation which may result in     */ 
/* memory leaks when a great number of allocations and deallocations are done.                         */ 
 

      


void InvertMatrix(  double *matrix,  double *inverse, int dim) 
{ 
  double parity=1.0; 
  int *index; 
  double sumInvM=0.0,temp; 
  int i,j,k; 
  double *col; 
  double matrixcopy[dim][dim];     /* avoiding repetitive dynamic allocation, change to use code elsewhere */  
  double determinant=0.0; 
  //  cout << "dim = " << dim << endl;
  index = (int *) GetBlock(dim*sizeof(int)); 
  col = (  double *) GetBlock(dim*sizeof(  double)); 
 
  MatCopy((  double *) matrixcopy[0], (  double *)matrix,dim,dim);  /* make a copy of the matrix to be inverted */ 
  memset(index,0,dim*sizeof(int)); 
 
  LUDecomp((  double *) matrix,dim, index, &parity);        /* Invert matrix00, destroying it, inverse written to inverse */ 
  for (j=0;j<dim;j++) parity*=(*(matrix+dim*j+j)); 
  if (parity==0.0) Invert=0; 
  for(j=0;j<dim;j++)  
    { 
      for(i=0;i<dim;i++) col[i]=0.0; 
      col[j]=1.0; 
      LUBackSub((  double*) matrix,dim,(int *) index,(  double*)col); 
      for(i=0;i<dim;i++) (*(inverse+dim*i+j))=col[i]; 
    } 
     
  MatCopy((  double*) matrix,(  double *) matrixcopy[0],dim,dim);  /* Copy the backup of the matrix00 into the destroyed 
								      matrix */ 
  free(col); 
  free(index); 
  memset(matrixcopy,0,(dim)*(dim)*sizeof(  double));      /* reinitialize matrixcopy*/   
} 
 

int main(int argc, char **argv) 
{
  int NCount = 20;
  int mystep = 1;
  int newsteps = 1;
  if (argc != 2) error("You need an integer argument!");
  mystep = atoi(argv[1]);
  if (mystep <= 0) error ("Argument should be an integer more than one.");
  int xr, xc;
  int nvar ,nsim, noutcomes, nsteps;
  int i,j;
  // all the parameters should reside in the input file: lhsmatrix
  ifstream inf("lhsmatrix");  
  if (inf.fail()) error ("File lhsmatrix do not opened "); 
  inf >> nvar >> nsim >> noutcomes >> nsteps; 
  if (inf.fail()) error ("One of the first 4 number in the file lhsmatrix is incorrect ");
  newsteps = (int)(nsteps/mystep);
  // nvar == number of variables/parameters to study
  // nsim == number of simulations/samples
  // noutcomes == number of figures/data obtained
  // nsteps == number of steps to run the simulation
  int arraycount[NCount+1];
  if (noutcomes > NUMPRCCOUTS)
    { 
      cout << "Error! The number of outcomes is greather than a constant=" <<  NUMPRCCOUTS  << endl;
      return 1;
    }
  // create the arrays:
  double TimeMatrix[nsteps+1][nvar+1];
  //  cout << "OK Time Matrix" << endl;
  matrixDB LHSMatrix[nvar+1][nsim+1];        /* 3-D Matrix -parameter and outcome values & ranks for each run*/ 
  // cout << "OK matrixDB" << endl;
  double   RankMatrix[nvar+1][nsim+1]; /* Matrix of ranks for parameters that vary  */ 
  // cout << "OK RankMatrix " << endl;
  double   CijMatrix[nvar+1][nvar+1];  /* Calculated from RankMatrix, Cij^1 used to calculate PRCCs */  
  // cout << "OK Cij Matrix " << endl;
  double   CijInverse[nvar+1][nvar+1];  /* Used in inversion of Cij in the calculation of the PRCCs */  
  // cout << "OK Inverse matrix " << endl;
  paramT   pArray[nvar+1];                    /* Array of parameter distribution information K parameters */ 
  // cout << "OK pArray " << endl;
  double   CalcTime[nsteps+1][nsim+1]; /* Array of calculated outcomes */ 
  // cout << "OK CalcTime " << endl;
  timepoint  OutComeData[nsim+1][newsteps+1];
  //  cout << "OK OutcomeData " << endl;
  int     NewI[nsim+1];  
  FILE  *prcc[noutcomes+1];
  char    OutName[256];
  memset(OutComeData,0,nsim*nsteps+1*sizeof(timepoint)); 
  memset(LHSMatrix,0,nvar*nsim*sizeof(matrixDB)); 
  memset(RankMatrix,0,nvar*nsim*sizeof(  double)); 
  memset(CijMatrix,0,nvar*nvar*sizeof(  double)); 
  memset(CijInverse,0,nvar*nvar*sizeof(  double)); 
  memset(NewI,0,nsim*sizeof(int)); 
  
  string name_outcomes[noutcomes+1];

  //  cout << "\n nvar=" << nvar << "\n nsim=" << nsim << "\n noutcomes=" <<noutcomes << "\n nsteps=" << nsteps << endl;
  for (xr=0; xr < nvar; xr++)
    { 
      for (xc=0; xc < nsim; xc++)
	{
	  inf >> LHSMatrix[xr][xc].value;
	  //  cout << LHSMatrix[xr][xc].value << "\t";
	  if (inf.fail()) error ("Incorrect number in file lhsmatrix, row=",xr+2);
	}
      //     cout << endl;
    }
  inf.close();
  // cout << "End of lhsmatrix !" << endl;
  ifstream lhsdata("lhsdata");
  int NumVarPar=0; 
  int NumPar = nvar;
  int NumSim = nsim;
  int FirsCall = 1;
  if (!lhsdata.is_open()) error("File lhsdata not found!", 1);
  lhsdata >> NumPar;
  // cout << NumPar  << endl;
  for (xr=0; xr < NumPar; xr++)
    {
      lhsdata >> pArray[xr].name >> pArray[xr].min >>  pArray[xr].peak >> pArray[xr].max >> pArray[xr].dType;
      if (lhsdata.fail())
	error("Incorrect value in lhsdata file");
      //cout << xr << endl;		       								   
    }
  lhsdata.close();
  // cout << "End of lhsdata" << endl;
  ifstream outcomenames("lhsoutcome");
  outcomenames >> xr;
  if (xr != noutcomes) error("xr != noutcomes !");
  for (xr=0; xr < noutcomes; xr++)
    {
      outcomenames >> name_outcomes[xr];
      // cout <<  name_outcomes[xr] << endl;
    }
  outcomenames.close();
  double mynumber;
  double lasttime=0.0;
  char filename[32];
  char filename2[32];
  char filename3[32]; 
  int myxr = 0;
  // loading into OutComeData
  for (xc =0; xc < nsim ; xc++)
    { 
      sprintf(filename, "%04u",xc+1);
      ifstream outcdata(filename);
      // cout << filename << ":" <<endl;
      if (!outcdata.is_open()) error("Run/File not found #", xc+1);
      myxr = -1;
      //	while(!outcdata.eof())
      //	 {  
      for (xr=0; xr < nsteps; xr++)
	{
	  outcdata >> mynumber; // time
	  if (outcdata.fail()) error ("Incorrect number in datafile ", xc+1);
	  if ((xr % mystep) == 0) 
	    { myxr++;
	    OutComeData[xc][myxr].time = mynumber; 
	    lasttime = mynumber;
	    }
	  // cout << xr << "=" << OutComeData[xc][xr].time << "\t";
	  for (i=0; i< noutcomes; i++)
	    {  
	      outcdata >> mynumber;
	      if (outcdata.fail()) error ("Incorrect number in datafile ", xc+1);
	      if ((xr % mystep) == 0) 
		OutComeData[xc][myxr].outcome[i] = mynumber;
	    }	 
	}
      nsteps = xr;
      outcdata.close();
    }
  nsteps = newsteps+1;
  // PlaceNextOutCome(Call); 
  // cout << "loading lhsmatrix " << endl;
  //  Place the outcome variable into LHS matrix, at time=Call, and outcome = 0
  // here it will be two large loops for n-outomes and time.
  int Call = 161;
  int CurrOutcome = 18;
  double df = nsim - 2;
  double total, mean, median, var,numer; 
  double stdev, prop, lower, upper;
  double mycurrenttime = 0.0;
  ofstream plotavg("plotavg.gp");
  ofstream plotprcc("plotprcc.gp");
  ofstream plotsignif("plotsignif.gp");
  ofstream plot3d("plot3d.gp");
  
  plotavg << "set term postscript color  10 " << endl;
  plotavg << "set output "<< (char)34 << "stats.ps" << (char)34 << endl;
  // plotavg << "set pointsize 0.500000 " << endl; 
  plotavg << "set xlabel 'Time(Units)' "<< endl;

  plotprcc << "set term postscript color  10 " << endl;
  plotprcc << "set output "<< (char)34 << "prccvals.ps" << (char)34 << endl;
  //  plotprcc << "set pointsize 0.500000 " << endl; 
  plotprcc << "set ytics 0.2 " << endl;
  plotprcc << "set yrange [-1:1] " << endl;
  plotprcc << "set xrange [0:"<<lasttime << "] " << endl;
  plotprcc << "set xlabel 'Time(Units)' "<< endl;

  plotsignif << "set term postscript color  10 " << endl;
  plotsignif << "set output "<< (char)34 << "signif.ps" << (char)34 << endl;
  //  plotprcc << "set pointsize 0.500000 " << endl; 
  //  plotsignif << "set ytics 0.2 " << endl;
  plotsignif << "set yrange [0:0.3] " << endl;
  plotsignif << "set ytics 0.01 " << endl;
  plotsignif << "set xrange [0:"<<lasttime << "] " << endl;
  plotsignif << "set xlabel 'Time(Units)' "  << endl;

  plot3d << "set contour base " << endl;
  plot3d << "set xrange [0:"<<lasttime << "] " << endl;
  plot3d << "set term postscript color  10 " << endl;
  plot3d << "set output "<< (char)34 << "contour.ps" << (char)34 << endl;
  plot3d << "set xlabel 'Time(Units)' "  << endl;
  plot3d << "set ylabel 'Outcome Value [Min-Max]' " << endl;
  plot3d << "set zlabel 'Count'" << endl;
  plot3d << "set hidden3d " << endl;
  plot3d << "set cntrparam levels 8" << endl;
  plot3d << "set ticslevel 1" << endl;

  // string namefile;			 

  if (df < 0.00000) df =0.0;
  for (CurrOutcome = 0; CurrOutcome < noutcomes; CurrOutcome++)
    {
      plotavg << "set title " << (char)34 << "Stats :" << name_outcomes[CurrOutcome] << (char)34 << endl; 
       plot3d << "set title " << (char)34 << "Frequency values of: " << name_outcomes[CurrOutcome] << (char)34 << endl;
      sprintf(filename, "prcc.%04u",CurrOutcome+1);
      //filename  = name_outcomes[CurrOutcome]; 
      ofstream fileout(filename);		


       sprintf(filename2, "contour.%04u",CurrOutcome+1);
       ofstream filetemp(filename2);	

       sprintf(filename3, "density.%04u",CurrOutcome+1);
       ofstream density(filename3);
      for (Call=1; Call < nsteps; Call++) 
	{
	  Invert=1;   
	  StrongEffect =0;
	  mycurrenttime = OutComeData[0][Call].time;
	  int k; 
	  double min, max;   
	  Invert = 1;
	  for(k=0;k<nsim;k++)   
	    /* copy the outcomes into the lhsmatrix and pick out the time "Call" that was used
	       for, each outcome */ 
	    { 
	      LHSMatrix[nvar][k].value=OutComeData[k][Call].outcome[CurrOutcome]; 
	    } 
	  //  for (k=0;k<nsim;k++) 
	  //  printf("Outcomes Placed: %lf %lf\n",LHSMatrix[nvar][k].value, OutComeData[k][Call].outcome[CurrOutcome]); 

	  //  UncertaintyAnalysis(); 
			 
	  total = 0.0; 
	  numer = 0.0; 
	  max = -999999999999.99; 
	  min =  999999999999.99; 
	  median = LHSMatrix[nvar][(nsim+1)/2].value; 
	  for(i=0;i<nsim;i++)  
	    { 
	      total+=LHSMatrix[nvar][i].value; 
	      numer += (double)(LHSMatrix[nvar][i].value *LHSMatrix[nvar][i].value); 
	      if (LHSMatrix[nvar][i].value >max) max=LHSMatrix[nvar][i].value; 
	      if (LHSMatrix[nvar][i].value < min) min=LHSMatrix[nvar][i].value; 
	    } 
	  mean =  total / nsim; 
	  var = (numer-nsim*mean*mean)/(nsim-1); 
	  stdev = 0.0;
	  if (var > 0.0)
	    stdev = sqrt(var)/sqrt((double)nsim); // of the sample
	  lower = mean - (1.96*stdev);
	  upper = mean + (1.96*stdev);
	  if (lower  < 0) lower = mean;
	  if (upper > mean*2) upper = mean;
	  // printf("Time\t\tMean\t\tMedian\t\tVariance\t\tMaximum\t\tMinimum\n"); 
	  // printf(fileout,"\n %d\t %lg \t %lg \t %lg \t %lg \t %lg ",Call,mean,median,var,max,min); 
	  //  cout <<  OutComeData[xc][xr].time;
	  fileout << endl << mycurrenttime << "\t" << mean << "\t" << median  << "\t" << var << "\t" << max << "\t" << min << "\t" << lower << "\t" << upper << "\t";
	  // printf("\tMean:\t\t%lg\n\tMedian:\t\t%lg\n\tVariance:\t%lg\n\tMaximum:\t%lg\n\tMinimum:\t%lg\n\n",mean,median,var,max,min );  
	  // fileout << endl << Call << "\t" << mean << "\t" << lower << "\t" << upper << "\t";
	  
	  // ArrayCount fill
	  double delta = (max - min)/NCount;
	  double ress = 0.0;
	  //  cout << "nvar = " << nvar << " ,d = " << delta << endl;
	   for(i=0;i<NCount;i++)  
		  arraycount[i]=0; 
	  for(i=0;i<nsim;i++)  
		 {
			ress = (LHSMatrix[nvar][i].value - min);
			ress = ress/delta;
			//cout << LHSMatrix[nvar][i].value << "\t" << ress << "\t" << (int)(ress/delta) <<endl;
			arraycount[abs((int)(ress))]++;
		 }
	   
	 
	  for(i=0;i<NCount;i++)  
	     {
		
		filetemp << mycurrenttime << "\t" << i <<"\t" << arraycount[i] << endl;
		 //cout << arraycount[i] << "\t";
		density << arraycount[i] << "\t";
	     }
	   density << endl;
	   filetemp << endl;
	  //void RankTieSort(void) 
	  
	  
	  //		
	  int          count,ties=0,first; 
	  double  ranksum, newrank; 
	  int          NumToSort; 		 
 
	  NumToSort=nvar+1;  /* We're doing PRCC analysis */ 
	  for(k=0;k<NumToSort;k++)  
	    { 
	      if((pArray[k].min != pArray[k].max) || (k==nvar))        
		/* if we have a varying parameter or its the outcome */ 
		{ 
		  for(i=0; i<nsim; i++) LHSMatrix[k][i].rank = 1.0;     /* initialize kth row ranks to 1.0 */ 
		  ties=0; 
		  for(i=1;i<nsim;i++)                                   
		    { 
		      first=1; 
		      for(j=0;j<i;j++)  
			{ 
			  if(LHSMatrix[k][i].value >= LHSMatrix[k][j].value) LHSMatrix[k][i].rank+=1.0; 
			  else LHSMatrix[k][j].rank+=1.0; 
			  if(LHSMatrix[k][i].value + TINY==LHSMatrix[k][j].value+TINY)  
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
		      for(j=0;j<nsim;j++) if(LHSMatrix[k][j].tieset==i)  
			{ 
			  count++; 
			  ranksum+=LHSMatrix[k][j].rank; 
			}  
		      newrank= (ranksum/count); 
		      for(j=0;j<nsim;j++) if(LHSMatrix[k][j].tieset==i) LHSMatrix[k][j].rank=newrank; 
		    }    
		}    
	    }     
  
	  //void CopyVaryPar(void) 
	  NumVarPar=0;
	  int  m,n; 
	  for(i=0; i<NumPar; i++) 
	    { 	 
	      for(j=0; j<NumSim; j++) RankMatrix[NumVarPar][j] = LHSMatrix[i][j].rank; 
	      NewI[NumVarPar]=i; 
	      NumVarPar++;	 
	    } 
     
	  for(m=0; m<NumSim; m++) RankMatrix[NumVarPar][m] = LHSMatrix[NumPar][m].rank; 
	  //  printf("\n %d  Outcome Values:\n ", 1); 
	  // for (m=0;m<NumSim;m++) printf("%d %f \t%lg\n",m, RankMatrix[1][m],LHSMatrix[1][m].value); 

	  //     printf("\n %d  Outcome Values:\n ", NumVarPar); 
	  //  for (m=0;m<NumSim;m++) printf("%d %f \t%lg\n",m, RankMatrix[NumVarPar][m],LHSMatrix[NumPar][m].value); 
	  //  printf(" \n \nRank matrix: \n"); 
	  // PrintMatrix((double *)RankMatrix[0],NumVarPar+1,NumSim);  
   
	  ofstream rnk("lhsranking");
	  //  cout << "Ranking Matrix " << endl;
	  for (i=0; i<NumPar+1; i++)
	    {
	      for  (j=0; j<NumSim;j++)
		rnk << RankMatrix[i][j] << "\t";
	      rnk << endl;
	    }
	  rnk.close();
	  //	for (i=1;i<6;i++) fflush(check[i]); 
		
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
 
	  //oid FillCijMatrix(void) 
	  // 
 
	  double      mu, sum[3], denominator; 
	  int              t,ll; 
	  Invert = 1;
	  mu = (  double) (1.0+NumSim)/2.0; 
	  for(i = 0; i< (NumVarPar+1); i++) 
	    { 
	      for(j=0; j< (NumVarPar+1); j++) 
		{ 
		  for(ll=0;ll<3;ll++) sum[ll] = 0.0; 
		  for(t = 0; t< NumSim; t++) 
		    { 
		      sum[0] += (RankMatrix[i][t])*(RankMatrix[j][t]); 
		      sum[1] += (RankMatrix[i][t])*(RankMatrix[i][t]); 
		      sum[2] += (RankMatrix[j][t])*(RankMatrix[j][t]); 
		    } 
		  for(ll=0;ll<3;ll++) sum[ll]-=NumSim*mu*mu; 
		  denominator = sqrt(sum[1]*sum[2]);                
		  if (denominator != 0.0) CijMatrix[i][j] = sum[0]/denominator; 
		  else 
		    { 
		      //printf("Error: Zero denominator at outcome time :%lf in Fill Cij\n"); 
		      Invert=0; 
		    }  
		} 
	    }      
	  int same=0;
	  int StrongPar=NumVarPar; 
	  int Parity=1;     
	  //  printf("Cij Matrix:\n"); 
	  // PrintMatrix((double *)CijMatrix[0],NumVarPar+1,NumVarPar+1); 
	  //  cout << "Calling InvertMatrix " << endl;
	  InvertMatrix((  double *) CijMatrix[0],(  double *)CijInverse[0],NumVarPar+1); 
	  //  cout << " Invert = " << Invert << endl;
	  // printf("Cij Inverse :\n"); 
	  // PrintMatrix((double *)CijInverse[0],NumVarPar+1,NumVarPar+1); 

	  if (Invert==0)
	    { 
	      // CheckStrongEffect();
	      //  int same; 
	  
	      //		printf("Singular matrix at time %d.\nChecking for strong effect parameter....\n",Call); 
	      StrongEffect=0; 
	  
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
	      //	if (StrongEffect==1) cout << " Strong effect parameter: " << pArray[NewI[StrongPar]].name << endl; 
	      //printf("Strong effect parameter:%s\n",pArray[NewI[StrongPar]].name); 
	      //	else printf("No strong effect parameter found\n");  
	    }
	  if (Invert==1) 
	    {
	  
	      // here inside CheckInverse(); 
	  
	      double CheckMatrix[NumVarPar+1][NumVarPar+1];  
	      memset(CheckMatrix,0,(NumVarPar+1)*(NumVarPar+1)*sizeof(  double)); 
	      MatMult((  double *) CijMatrix[0],(  double *)CijInverse[0],(  double *)CheckMatrix[0],NumVarPar+1,NumVarPar+1,NumVarPar+1); 
	      Invert=1; 
	      for(i=0; i<(NumVarPar+1);i++) 
		{ 
		  for(j=0;j<(NumVarPar+1);j++) 
		    {    
		      if (i!=j) if ((-.000001>CheckMatrix[i][j]) || (CheckMatrix[i][j]>0.000001))Invert=0; 
		      if (i==j) if ((0.9999>CheckMatrix[i][j]) || (CheckMatrix[i][j]>1.000001)) Invert=0; 
		    } 
		} 
	      if (Invert==1) ;//printf("Matrix inversion successful\n"); 
	      else 
		{
		  //CheckStrongEffect(); here
		  //printf("Singular matrix at time %d.\nChecking for strong effect parameter....\n",Call); 
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
	    
		  //	 if (StrongEffect==1)  cout << " Strong effect parameter: " << pArray[NewI[StrongPar]].name << endl;

		  //  printf("Strong effect parameter:%s\n",pArray[NewI[StrongPar]].name); 
		  // else printf("No strong effect parameter found\n");  

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
 
	  // void GetPRCC(void) 
	  // { 
	  double PRCCNum, PRCCDenom,PRCC,TTest ; 
	  //int i;  
     
	  int DoGlobal = 1;
	  int FirstCall = 1;

	  //  if (DoGlobal==1) printf("\nIf PRCC=1.0 or-1.0,t-test is effectively infinite\n"); 
	  //cout << "Strong Effect " <<StrongEffect <<  endl;
	  //	 cout << endl << "CurrOutcome = " << CurrOutcome  << " Call = " << Call << " \t";
	  if ((Invert==1) || (StrongEffect==1)) 
	    {
	      for(i=0; i< NumVarPar; i++)                  /* 1 */ 
		{ 
		  //	 	 cout << i << "\t";
		  TTest = 0.0;
		  PRCC=0.0; 
		  if ((StrongEffect !=1) && (Invert==1))            /*2*/ 
		    { 
		      PRCCNum = -1*(CijInverse[i][NumVarPar]); 
		      if((CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar])< 0.000) 
			{    
			  if((CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar])<-0.01) 
			    error("PRCC denominator <0 in GetPRCC-no square root"); 
			  else 
			    { 
			      PRCCDenom=-1*(CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
			      // printf("%g \t %g \n", PRCCNum,CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
			      // printf("PRCC Denom close to zero and negative, opposite sign given, check warning file"); 
			    } 
			}		
		      PRCCDenom= sqrt(CijInverse[i][i]*CijInverse[NumVarPar][NumVarPar]); 
		      if (PRCCDenom!=0.0000) 
			PRCC = PRCCNum/(PRCCDenom);
		      else { 
			PRCC = 999999;			 error("PRCCDenom is zero in GetPRCC"); }
		      if ((-1.0<=PRCC) && (PRCC<=1.0)) 
			TTest = fabs(PRCC)*(sqrt((NumSim-2)/(1-(fabs(PRCC)*fabs(PRCC)))));
		      else
			{ TTest = 999999;
			PRCC =    999999; }
		      //else error("Error in Calculating T-test Value"); 
		      //fileout << " !! PRCC, TTest = " << PRCC << "\t" << TTest << endl;
		      // fileout << "-\t-\t";
		    }					 
		  if ((StrongEffect==1) && (i==StrongPar)) PRCC=(double)Parity*1.0;
		  prop = 0.0;
		  if (PRCC <= 99999) {
		    if ((df+(TTest*TTest))!=0.000000)
		      prop = betai(0.5*df, 0.5, df/(df+(TTest*TTest)));
		    else 
		      prop = 0.0;
		    if (finite(PRCC)) 
		      fileout <<  PRCC << "\t";
		    else 
		      fileout << "-\t";
		    // fileout << TTest << "\t";
		    if (finite(TTest)) 
		      fileout << prop << "\t";
		    else
		      fileout << "-\t";
		  }					 
		  // if ((StrongEffect==1) && (StrongPar!=i))
		  //	 cout << pArray[NewI[i]].name << "\t 0.0 \t" << 0.0 << "\t";
		  //		cout << " 0.0 \t" << 0.0 << "\t";
		  // if (((StrongPar==i) || (StrongEffect!=1)) && (PRCC>=0.00))
		  //	 cout << pArray[NewI[i]].name << "\t" << PRCC << "\t" << TTest << "\t";
		  //		cout <<  PRCC << "\t" << TTest << "\t";
		  //if (((StrongPar==i) || (StrongEffect!=1)) && (PRCC<0.00)) 
		  //	 cout <<  pArray[NewI[i]].name << "\t" << PRCC << "\t" << TTest << "\t";
		  //	cout <<  PRCC << "\t" << TTest << "\t";
		}
	    }
	  else
	    {
	      for(i=0; i< NumVarPar; i++)      
		fileout << "-\t-\t";
				  
	      //				  cout << "n";
	    }
				
	  // plot 'prcc.0012' u 1:13, 'prcc.0012' u 1:14 t 'P' w d, 'prcc.0012' u 1:11, 'prcc.0012' u 1:12 t 'P' w d, 'prcc.0012' u 1:9 t 'A'
	}
      fileout.close();
      plotprcc <<"set title " << (char)34 << "PRCC : " << name_outcomes[CurrOutcome] << (char)34 << endl; 
      plotprcc << "plot ";
      for(i=0; i< NumVarPar; i=i+1)   
		  {
	  plotprcc << "'" << filename << "' u 1:" << (i*2)+9 <<  " t '" << pArray[NewI[i]].name << "', ";
		  }
      plotprcc << "0 t '' w l " << endl;


      plotavg << "plot '"<<filename << "' u 1:2:7:8 " << "t 'Avg' w errorbars, '" << filename <<"' u 1:6 t 'Min' w l, '" << filename <<"' u 1:5 t 'Max' w l" << endl;
      

      plotsignif <<"set title " << (char)34 << "PRCC Significance: " << name_outcomes[CurrOutcome] << (char)34 << endl; 
      plotsignif << "plot ";
      for(i=0; i< NumVarPar; i=i+1)   
		  {
			 plotsignif << "'" << filename << "' u 1:" << (i*2)+10 <<  " t '" << pArray[NewI[i]].name << "', ";
		  }
      plotsignif << "0 t '' w l, 0.01 t '' w l " << endl;

		filetemp.close();
		density.close();
		plot3d << "splot '" << filename2 << "' t '' w l" << endl; 

      
    }
  //  double p1 = betai(0.5*120, 0.5, 2);
  // double p2 = betai(0.5*120, 0.5, 0.99);
  // double p3 = betai(0.5*120, 0.5, 0.999);
  // cout << p1 << " \t " << p2 << "\t " << p3 << endl;
  plotavg << "exit" << endl;
  plotavg.close();
  plotprcc << "exit" << endl;
  plotprcc.close();
  plotsignif << "exit" << endl;
  plotsignif.close();
  plot3d << "exit" << endl;
  plot3d.close();
  myreturn = system("gnuplot plotprcc.gp"); if (myreturn != 0) { cout << "Error en la llamada a gnuplot plotprcc.gp"; }
  myreturn = system("gnuplot plotavg.gp");  if (myreturn != 0) { cout << "Error en la llamada a gnuplot plotavg.gp"; }
  myreturn = system("gnuplot plotsignif.gp"); if (myreturn != 0) { cout << "Error en la llamada a gnuplot plotsignif.gp"; }
  myreturn = system("gnuplot plot3d.gp"); if (myreturn != 0) { cout << "Error en la llamada a gnuplot plot3d.gp"; }
  cout << "OK!, prccvals.ps, stats.ps, signif.ps, and contour.ps created. " << endl;
  return 0;  
}






