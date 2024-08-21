#include <stdio.h>
#include <stdlib.h>
#include "datatype.h"
/* Modified on 12/30/02 */

/* This subroutine writes the gnuplot script that is fed to the external UNIX graphing program gnuplot 
      */  
/* The script that is created is called plot.scr, and if KeepGNU is set to 1, it is kept after the    */  
/* the data is plotted .                                                                                         */  
/* Variables that are used in creating this script are:                                                          */  

#define SMALL_NUMBER 0.001
#define LARGE_NUMBER 1e30


typedef struct {char name[50];  double max, min; } eqfig; 

double gradf (double g, double c)
{
  return (g*g)/((g*g)+c);
}

double ratio1 (double a, double b)
{
  if ((a + b) > SMALL_NUMBER)
    return (a /(a + b));
  else
    return 0.0;
}

double ratio2 (double a, double b, double p)
{
  if (b > SMALL_NUMBER)
    return ((a/b)/((a/b) + p));
  else
    return 0.0;
}


void WritePlotScript(int GraphMode, int NumEq, int nrows, eqfig Eqns[], int myrun, int printlhs)  
{    
  int     i;     
  FILE*   plotscript;   /* The file with script that plots the outcome graphs */  
  char graphfile[30] = "output1.pdf";
  char filename[15];
  float PointSize = 0.5;
  char TimeScale[30] = "Units";
  char plotline[30]="points";
    plotscript=fopen("plot.scr","w");  
  
  if ((GraphMode==0) || (GraphMode==1))   
     fprintf(plotscript,"set terminal pdf \n");
  // era ps fprintf(plotscript,"set term postscript color  10\n");  
 
  if ((GraphMode==2))   
      fprintf(plotscript,"set term postscript eps %s 10\n","portrait");  
  
  fprintf(plotscript," set output %c%s%c\n",'\"',graphfile,'\"');  
/*  fprintf(plotscript,"set key outside \n"); */
/*   fprintf(plotscript,"set key below  Left title 'Legend' box 3 \n");			*/
  if (printlhs > 0)
  fprintf(plotscript, "load \'legend.gp\' \n");
  fprintf(plotscript,"set pointsize %lf\n",PointSize);  
  fprintf(plotscript,"set xlabel %cTime(%s)%c\n",'\'',TimeScale,'\'');  
  
  for(i=1;i<=NumEq;i++)  
  {  
    fprintf(plotscript,"set title %c%s%c\n",'\"',Eqns[i].name,'\"');  
    /*  fprintf(plotscript,"show title \n");  */
    if ((Eqns[i].max - Eqns[i].min) <= 0.000001)
      Eqns[i].min = Eqns[i].min-0.00001;    
    if (printlhs <= 0)
      {
	fprintf(plotscript,"set yrange [%f:%f]\n",Eqns[i].min, Eqns[i].max); 
	fprintf(plotscript,"plot \'0001\' using 1:%d title %c%s %c with %s \n",i+1,'\'',Eqns[i].name,'\'',plotline); 
      }
    else
      { 
	sprintf(filename, "%04u.gp",i);
	fprintf(plotscript, "load \'%s\' \n",filename);} 
  }
  fprintf(plotscript,"print \'File output1.pdf created!' \n");
  fprintf(plotscript,"exit\n");  
  fclose(plotscript);  
  
}  
  
void WritePlotScript1(int GraphMode, int NumEq, int nrows,  eqfig Eqns[], int myrun, int printlhs)  
{    
  int     i;     
  FILE*   plotscript;   /* The file with script that plots the outcome graphs */  
  char graphfile[30] = "output1.pdf";
  float PointSize = 0.5;
   char filename[15];
  char TimeScale[30] = "Units";
  char plotline[30]="points";
    plotscript=fopen("plot.scr","w");  
  
  if ((GraphMode==0) || (GraphMode==1))   
   fprintf(plotscript,"set terminal pdf \n");
   // fprintf(plotscript,"set term postscript color  10\n");  
  
  if ((GraphMode==2))   
      fprintf(plotscript,"set term postscript eps %s 10\n","portrait");  
  
  fprintf(plotscript," set output %c%s%c\n",'\"',graphfile,'\"');  
  /*fprintf(plotscript,"set key outside \n"); */
  /*    fprintf(plotscript,"set key below  Left title 'Legend' box 3 \n"); */
  if (printlhs > 0)
     fprintf(plotscript, "load \'legend.gp\' \n");
  /* fprintf(plotscript,"set key outside \n"); */ 

  fprintf(plotscript,"set pointsize %lf\n",PointSize);  
  fprintf(plotscript,"set xlabel %cTime(%s)%c\n",'\'',TimeScale,'\'');  
  
  for(i=1;i<=NumEq;i++)  
  {  

    fprintf(plotscript,"set title %c%s%c\n",'\"',Eqns[i].name,'\"');  
    /* fprintf(plotscript,"show title \n");  */
    if ((Eqns[i].max - Eqns[i].min) <= 0.000001)
      Eqns[i].min = Eqns[i].min-0.00001;
    if (printlhs <=0)
      {
      fprintf(plotscript,"set yrange [%f:%f]\n",Eqns[i].min, Eqns[i].max); 
      fprintf(plotscript,"plot \'0001\' using 1:%d title %c%s %c with %s \n",i+1,'\'',Eqns[i].name,'\'',plotline);  
      }
    else
      {
	sprintf(filename, "%04u.gp",i);
	fprintf(plotscript, "load \'%s\' \n",filename); 
      }
  }
  fprintf(plotscript,"print \'File output1.pdf created!' \n");
  fprintf(plotscript,"exit\n");  
  fclose(plotscript);    
}  
