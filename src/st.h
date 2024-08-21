#include <string>
using namespace std;
struct st_struct {
  string   name; // original name   A 
  string   eqname;  // equation name  EQ[1]
  string   constname;  // name for constants C[1] 
  int      object;
  int      type;         // attribute
  int      var_type;     // attribute for variables
  int      initialized;  // 0, not initilized
  string   dfequation;   // here to store df equations if any 
  string   constequation;  // here to store constant equations if any      
  int     p_diffeq;    // pointer to list of diff eqn
  int     declaration;   // 1=true, 0, false
  int     declared_here;  // pointer where it was declared
  double  begin_range;  // for contants : end and begin of range of values
  double  mid_range;    // for bias 
  double  end_range;
  string  type_range;   // for types: U == Uniform, T = Triangular
  int     lhsrun;  // 0 = norun, 1=runlsh.
};

struct st_mycode {
  int variable; // pointer to the variable
  string linecode;
};


struct dataxy {
  int lastid;
  int lasttok;
};


void FillHypercube(int element, double LHSMatrix[][MAX_LHS_MATRIX], st_struct pArray[], int NumSim);
void FillTriangularHypercube(double min, double peak, double max,int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim) ;
void FillUniformHypercube( double distmin, double distmax, int element, double LHSMatrix[][MAX_LHS_MATRIX], int NumSim);
void FillIntervalHypercube( double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim ) ;
void FillLogHypercube(double distmin, double distmax,int element,  double LHSMatrix[][MAX_LHS_MATRIX], int NumSim );

