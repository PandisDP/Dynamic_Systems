#include "stdio.h"
#include "datatype.h"

/*   */
#define MAX_FIGS 250

typedef FILE* stream;
  
typedef struct {char name[50];  double max, min; } eqfig; 
float bessj1(float x);
float bessj0(float x);
float bessj(int n, float x);

void ludcmp(double **a, int n, int *indx, double *d);
void lubksb(double **a, int n, int *indx, double b[]);

void derivs(double TimeAct, double y[], double dydx[], double C[], double AUX[]) ;

int odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
	double hmin, int *nok, int *nbad,
	void (*derivs)(double, double [], double [], double C[], double AUX[]),
	int (*rkqs)(double [], double [], int, double *, double, double, double [],
	double *, double *, void (*)(double, double [], double [], double [], double[]),
	double [], double []), double C[], double AUX[], double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, stream myfile);


int rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
			 double yscal[], double *hdid, double *hnext,	void (*derivs)(double, double [], double [], double[], double[]), double C[], double AUX[]);


int  rkck(double y[], double dydx[], int n, double x, double h, double yout[],
			 double yerr[], void (*derivs)(double, double [], double [], double[], double []), double C[], double AUX[]);


void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
			 double yout[], void (*derivs)(double, double[], double[], double C[], double AUX[]), double C[], double AUX[]);

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);

int bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
				void (*derivs)(double, double [], double [], double C[], double AUX[]), double C[], double AUX[]);

void rzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);

int bsstep_rz(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
			       void (*derivs)(double, double [], double [], double C[], double AUX[]), double C[], double AUX[]);


int load_data_figures(double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, int j, double y[], double C[], double AUX[], double x, stream myfile); 

void jacobn(double x, double y[], double dfdx[], double **dfdy, int n,  double C[], double AUX[]) ;

int stiff(double y[], double dydx[], int n, double *x, double htry, double eps,
	  double yscal[], double *hdid, double *hnext,
	  void (*derivs)(double, double [], double [], double C[], double AUX[]), double C[], double AUX[]);

void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[],
	   void (*derivs)(double, double [], double [],double C[], double AUX[]),double C[], double AUX[] ) ;


int stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext,
				void (*derivs)(double, double [], double [] , double C[], double AUX[] ), double C[], double AUX[] );


/* void update_variables_lhs(int * c, double AUX[], double y[], double C[]); */


double gradf (double g, double c);

double ratio1 (double a, double b);

double ratio2 (double a, double b, double p);

void rk4( double y[],  double dydx[], int n,  double x,  double h,  double yout[],
	  void (*derivs)( double,  double [],  double [],  double[],  double[]),  double C[],  double AUX[]);

void WritePlotScript(int GraphMode, int NumEq, int nrows, eqfig Eqns[], int myrun, int printlhs)  ;

void WritePlotScript1(int GraphMode, int NumEq, int nrows,  eqfig Eqns[], int myrun, int printlhs)  ;


