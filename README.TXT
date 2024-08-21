crunOrdinary Differential Equations Solver

ODES

This is the program solves ordinary differential equations with 
any of the following methods:

1.  Runge-Kutta Method ($Rk4)
2.  Adaptive Stepsize Control for Runge-Kutta  ($AdaptiveRk4)
3.  Polynomial  Extrapolation Bulirsch-Stoer Step Method ($BSStepPoly)
4.  Rational Extrapolation Bulirsch-Stoer Step Method ($BSStepRational)
5.  Stiff integrate stiff ODEs by fourth-order Rosenbrock ($StiffRosenbrock)
6.  Stiff: Semi Implicit Extrapolation Method ($StiffBSStep )

The following software is required: 

install: csh, g++, flex++, gnuplot, gv

sudo apt-get install csh
sudo apt-get install flex
sudo apt-get install bison
sudo apt-get install gnuplot
sudo apt-get install gv
sudo apt-get install xfonts-100dpi
sudo apt-get install xfonts-75dpi

then:

>./build

and the program "teval1" will be generated.
and to install the program in /usr/local/ :

>su 
>./install

or 

> sudo ./install


Then you feed different files to the solver with:

>cd test
>crun <model_filename>

where <model_filename> is a text file with the model description:

In the directory there are several test files (test09.txt  ..  testTB.txt
that can be used as examples).

1. Create a model description file (see example below)
2. Compile your model with crun (for GNU CC), or icrun (for Intel C++
compiler, in Linux):

> crun mymodel1.txt      // for gnu c compiler, or

> icrun mymodel1.txt    // for Intel C++ compilear


PRCC
----

If you want to obtain the PRCC values, you should run the model with
the $Rk4 method, and include the $lhs command. The program will create several
"00*" output and "lhs*" files. 

Call the PRCC program from the prompt:

> prcc N

Where N is an integer value, sampling each Nth row of the output files. 
If your computer does not have enough memory, the program "prcc"
will fail ("Segmentation fault"), then you should increase the value of N
(usually:  2 <= N <= 20).
  


// The input files should have the following instructions:

---------------------
// Examples of comments :
// Comment
% Comment
/* Comment */

$Initial_conditions              // You should have this line
                                 // Declare the variable and initial value 
T0 = 0.05 0.0 0.05;                 
T1 = 0.0 0.0 0.0;
T2 = 0.0 0.0 0.0;
T8= 0.0  0.0 0.0;
MR= 300.0  0.0 300.0;
MI= 0.0  0.0 0.0;
f1=0;
f2=0;
f3=0;
f4=0;
f5=0;
f6=0;
TINY = 1e-200;
f7= 10.0  0.0 10.0;
f8= 5.0   0.0 5.0;

$Equations			// You should have this line
                        
f1 = r4g*(((1.0/p1)-1.0)/(p1+1.0)) ; 
f2 = (1.0/r4g)*(((1.0/p2)-1.0)/(p2+1.0)) ;
f3 = r4g*(((1.0/p3)-1.0)/(p3+1.0));
Pw = 2.0;
int1 = IG/(IG+f1*I4+sc1);
TB = BI+BE;
T = T0+T1+T2+T8;
M = MI+MR+MA;
                              // the onese in derivative form:
T2' = x2+x7-mut2*T2 ;
T8' = t8recruit-t8death;
MR'  = sr+alpha4*(MA+w*MI)+x18+x5-x13-x3-mur*MR;
MI' =  x13-x12-(x12a)-x8-x52-mur*MI;
MA' = x3-x5-mur*MA;
IG' = ((TimeAct<(500.0))?1.0:1.0)*(x17+alpha5a*T1*(x24a))-muig*IG;
I12' = ((TimeAct<(500.0))?1.0:1.0)*(alpha8*MA+alpha22*MI*(IG/(IG+f7*I10+c22))+alpha23*MR*(x14/(x14+c23)))-mui12*I12 ;
I4' =  ((TimeAct<(500.0))?1.0:1.0)*(alpha11*T0+alpha12*T2+alpha12a*T8)-mui4*I4;
I10' = ((TimeAct<(500.0))?1.0:1.0)*(x19+alpha16*T1+alpha17*T2+alpha18*(I12/(I12+c2))*T0)-mui10*I10;
FA'  = ((TimeAct<(1000.0))?1.0:1.0)*(alpha30*MI+alpha31*MA*(x50/(x50+x51+s10))+(alpha32*T1+alpha33*T8)*(x14/(x14+c10a)))-mufa*FA;
BI' = x16+x11-x12*N-x8*N-mui*BI-N*(x12)-N*(x12a)-x52*N;
BE' = x8*N+alpha20*BE-x11-k15*(MA*BE)-x10+N2*x52+mui*BI+N1*(x12);

// Simulation parameters:
// You should have the following lines:

$StartTime 0.0                    
$StopTime  500 
$timestep  0.001
$timefigs  1                  
$plot TB T M MIp MITIratio;     // plot the variables, do not forget the SEMICOLON.

$Method $Rk4                    // here you do not need semicolon!

// Other Methods are:

// $AdaptiveRk4
// $BSStepPoly
// $BSStepRational
// $StiffRosenbrock
// $StiffBSStep

$lhs 20 alpha1 alpha2;        // for LHS
$seed $clock   		      // To set the clock with an integer number or $clock

------------------------------


