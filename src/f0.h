#include <stdio.h>

#define NUMBER              3100
#define NUMBER_OF_EQUATIONS 3110
#define EQUATIONS           3112  
#define PLOT                3130
#define INITIAL_CONDITIONS  3140 
#define PARAMETERS          3150
#define PARAMETERS_FILE     3160
#define METHOD              3170
#define RK4                 3180
#define IDENTIFIER          3190
#define APOSTROPHE          3191
#define LEVEL               1001
#define AUX                 1002
#define CONST               1003
#define ASSIGNMENT          1005
#define OPERATOR            1006 
#define KEYWORD             1007
#define FUNCTION            1008
#define PAR_LEFT            1009
#define PAR_RIGHT           1010
#define SEMI_COLON          1011
//#define IF                  1012
//#define ELSE                1013
#define END                 1014
#define UNIFORM             1015
#define TRIANGULAR          1016
#define CINTERVALTYPE       1017
#define LOGTYPE             1018
#define ALLVARS             1020
#define LHSRUN              1021
#define SEED                1022
#define CLOCK               1023


#define FROM_T              2002
#define TO_T                2003
#define TIMESTEP            2004
#define TIMEFIGS            2005

#define MAXTS               50000  /* number of variables */
#define MAX_LHS_MATRIX      1000    /* number of simulations for lhs */
#define MAX_KEY_TS             1
#define MAXLINE             1000
#define MAX_FIGS            50000  /* max number of variables to plot */

/* methods */

#define RungeKutta4th          5001     /* $Rk4 */
#define CashKarpRungeKutta     5002     /* $AdaptiveRk4 */
#define BSStepPoly             5004     /* BulirshStoerStepPoly == $BSStepPoly */
#define BSStepRational         5005     /* $BSStepRational */ 
#define StiffRosenbrock           5007    /* $RosenStiffJacob */
#define StiffSemiImplicit         5008
#define StiffBSStep               5009

/*

# [16.1] rk4 integrate one step of ODEs, fourth-order Runge-Kutta (example)
# [16.1] rkdumb integrate ODEs by fourth-order Runge-Kutta (example)
# [16.2] rkqs integrate one step of ODEs with accuracy monitoring (example)
# [16.2] rkck Cash-Karp-Runge-Kutta step used by rkqs
# [16.2] odeint integrate ODEs with accuracy monitoring (example)
# [16.3] mmid integrate ODEs by modified midpoint method (example)
# [16.4] bsstep integrate ODEs, Bulirsch-Stoer step (example)
# [16.4] pzextr polynomial extrapolation, used by BSSTEP (example)
# [16.4] rzextr rational function extrapolation, used by BSSTEP (example)
# [16.5] stoerm integrate conservative second-order ODEs (example)
# [16.6] stiff integrate stiff ODEs by fourth-order Rosenbrock (example)
# [16.6] jacobn sample Jacobian routine for stiff
# [16.6] derivs sample derivatives routine for stiff
# [16.6] simpr integrate stiff ODEs by semi-implicit midpoint rule (example)
# [16.6] stifbs integrate stiff ODEs, Bulirsch-Stoer step (example)

*/
