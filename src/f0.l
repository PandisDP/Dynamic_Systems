/* Program for Flex++ to manage lexical analisis
   Programmer: Jose L. Segovia-Juarez
   2002-2022  */
%option noyywrap
%{
//#include <FlexLexer.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "f0.h"
#include "st.h"
using namespace std;
dataxy t;

//int yywrap();
int insert_st(st_struct st[], char * yytext);
int look_st(st_struct st[], char * yytext);
int insert_st(st_struct st[], string yytext);
int look_st(st_struct st[], string yytext);

void error(const char * error, int tok);
void error(std::string, int linen);
void error(string ab, int linen);
void error(string ab);

void warning(string ab);
st_struct st[MAXTS];   // store for variables/constants
st_mycode mc[MAXTS];   // array to store the lines of code
int levelvar[MAXTS];   // array of indexes to replace variables in diff eqn
int mylineno=1;
int lastid=1;
int nconstants=0;
int codelineno=0;
long double lastnumber;
int nomore = 0;
char myeq0[MAXLINE];
string myvar1;
int myreturn;

// FOR FIGURES
int nfigs =0;
int tfigs[MAXTS];
int LHSDecl[MAX_LHS_MATRIX];
// int nlhs=0;
%}

alpha    [A-Za-z_]
dig      [0-9]
id       ({alpha})({alpha}|{dig})*
num1     -?{dig}+\.?([eE][-+]?{dig}+)?
num2     -?{dig}*\.{dig}+([eE][-+]?{dig}+)?
number   {num1}|{num2}
idarray  {id}[(]{dig}+[)]
%%

{number}            {
                    lastnumber = atof(yytext);
                    return NUMBER;  }

$[Ee]quations           { return EQUATIONS; }
$[Pp]lot                { return PLOT; }
$[Ii]nitial_[Cc]onditions  { return INITIAL_CONDITIONS; }
$[Pp]arameters_file     { return PARAMETERS_FILE; }
$[Mm]ethod              { return METHOD;}
$[Ss]tartTime           { return FROM_T; }
$[Ss]topTime            { return TO_T; }
$[Tt]imestep            { return TIMESTEP; }
$[Tt]imefigs            { return TIMEFIGS; }
$[Rr]k4                 { return RungeKutta4th; }
$[Aa]daptiveRk4         { return CashKarpRungeKutta;}
$BSStepPoly             { return BSStepPoly; }
$BSStepRational         { return BSStepRational; }
$StiffRosenbrock        { return StiffRosenbrock; }
$StiffSemiImplicit      { return StiffSemiImplicit; }
$StiffBSStep            { return StiffBSStep; }
$[Uu]                   { return UNIFORM; }
$[Tt]                   { return TRIANGULAR; }
$[Cc]                   { return CINTERVALTYPE; }
$[Ll]                   { return LOGTYPE; }
$[Ll][Hh][Ss]           { return LHSRUN; }
$[Aa]ll[Vv]ars          { return ALLVARS; }
$[Ss]eed                { return SEED;  }
$[Cc]lock               { return CLOCK; }
";"                     { return SEMI_COLON; }

"sin"|"cos"|"tan"|"sincos"|"csin"|"ccos"|"ctan"|"asin"|"acos"|"atan"|"atan2"|"casin"|"cacos"|"catan"|"exp"|"exp2"|"exp10"|"pow10"|"log"|"log10"|"log2"|"logb"|"ilog"|"pow"|"sqrt"|"cbrt"|"hypot"|"expm1"|"log1p"|"cexp"|"clog"|"clog10"|"csqrt"|"cpow"|"sinh"|"cosh"|"tanh"|"csinh"|"ccosh"|"ctanh"|"asinh"|"acosh"|"atanh"|"casinh"|"cacosh"|"catanh"|"erf"|"erfc"|"lgamma"|"j1"|"y0"|"y1"|"yn"|"fabs"|"frexp"|"ldexp"|"ceil"|"floor"|"trunc"|"rint"|"nearbyint"|"round"|"lrint"|"llround"|"modf"|"fmod"|"drem"|"fmin"|"fmax"|"fdim"|"fmax"|"ratio1"|"ratio2"|"gradf"																																																																																							{ return FUNCTION; }

{idarray}             {
  int pos1;
  myvar1 = yytext;
  pos1 = myvar1.find('(',0);
  if (pos1 == 0) error("should be parleft at line",mylineno);
  myvar1.replace(pos1,1,"_");
  pos1 = myvar1.find(')',0);
  if (pos1 == 0) error("should be parclose at line",mylineno);
  myvar1.replace(pos1,1,"_");
  if (nomore==1) lastid = look_st(st,myvar1);
  else
    {
      lastid = insert_st(st, myvar1);
      st[lastid].type = IDENTIFIER;
    }
  return IDENTIFIER;
  }

{id}                    {
  if (nomore==1) lastid = look_st(st,yytext);
  else
    {
      lastid = insert_st(st, yytext);
      st[lastid].type = IDENTIFIER;
    }
  return IDENTIFIER;
   }

"'"                      { return APOSTROPHE; }

=                {  return ASSIGNMENT ; }

"{"[^}\n]*"}"     /* eat up one-line comments */

[ \t\r\f\b]+          /* eat up whitespace */

\n              mylineno++;

"//"|"%"    {  int c;
                   while ( (c = yyinput() ) != 0)
                   {  if (c == '\n')
                     { ++mylineno;
                       break;
                      }
                }
          /*  printf("Line %d comment", mylineno);  */
        }

"/*"    {
                   int c;
                   while((c = yyinput()) != 0)
                       {
                       if(c == '\n')
                           ++mylineno;
                       else if(c == '*')
                           {
                           if((c = yyinput()) == '/')
                               break;
                           else
                               unput(c);
                           }
                       }
                   }

","|"("|")"|"+"|"-"|"/"|"*"|">"|">="|"=="|"<"|"<="|"!="|"?"|":"|"&&"|"||"  { return OPERATOR; }

.          { printf( "Line %d:  Unrecognized character: %s\n", mylineno, yytext );
 return(yytext[0]);  error("Unrecognized character in line ", mylineno);
exit(1); }

%%
// using namespace std;

//int yywrap()  { return 1; }

void init_st(st_struct st[])
{
  int i=1;
  for (i=0; i<MAXTS; i++)
	 {
      st[i].name = "";
      st[i].dfequation = "";
      st[i].constequation = "";
      st[i].eqname = "";
      st[i].constname = "";
      st[i].object = MAX_KEY_TS;      // number of entries on the table
      st[i].type = 0;
      st[i].var_type = 0;
      st[i].p_diffeq = 0;
      st[i].initialized = 0;
      st[i].declaration = 0;
      st[i].declared_here=0;
      st[i].begin_range = 0;
      st[i].mid_range  = 0;
      st[i].end_range = 0;
      st[i].type_range = "";
      st[i].lhsrun = 0;
      levelvar[i] = 0;
    }
    st[0].name = "TimeAct";
    st[0].constname = "TimeAct";
    st[0].eqname = "TimeAct";
}

int look_st(st_struct st[], string yytext)
{ int i, foundhere;
   string a;
   foundhere=0;
   a = yytext;
   for (i=0; i<=t.lastid; i++)
	 {
	     if (a == st[i].name)
                  {
		     return i;
		  }
	 }
   return 0;
}

int look_st(st_struct st[], char * yytext)
{ int i, foundhere;
   string a;
   foundhere=0;
   a = yytext;
   for (i=0; i<=t.lastid; i++)
	 {
	     if (a == st[i].name)
                  {
		     return i;
		  }
	 }
   return 0;
}

int insert_st(st_struct st[], string  yytext)
{
   int i, foundhere;
   string a;
   foundhere=0;
   a = yytext;
   foundhere = look_st(st, yytext);
   if (foundhere == 0)
	{
	  if (t.lastid >= MAXTS) error("Symbol table too small, please increase the constant MAXTS, error found at ",mylineno);
	  t.lastid++;
	  st[t.lastid].name = a;
	  st[t.lastid].declared_here = foundhere;
	  //  cout << "A Inserting " << a << endl ;
	  return t.lastid;
   }
   else { return foundhere;  }
}

int insert_st(st_struct st[], char * yytext)
{
   int i, foundhere;
   string a;
   foundhere=0;
   a = yytext;
   foundhere = look_st(st, yytext);
   if (foundhere == 0)
     {
       if (t.lastid >= MAXTS) error("Symbol table too small, please increase the constant MAXTS, error found at ",mylineno);
       t.lastid++;
       st[t.lastid].name = a;
       st[t.lastid].declared_here = foundhere;
       // cout << "B Inserting " << a << endl ;
       return t.lastid;
     }
   else { return foundhere;  }
}

void error(const char * ab, int tok)
{
  printf("Error: %s %d \n", ab, tok);
  exit(1); }

void error(string ab, int linen)
{ cout << "Error: " << ab  << " " << linen << endl;
 exit(1); }

void error(string ab)
{ cout << "Error: " << ab  << endl;
 exit(1); }

void warning(string ab)
{ cout << "Warning: " << ab  << endl; }


#include "lhs.c"

int main( int argc, const char * * argv )
{
    ++argv, --argc;  /* skip over program name */
    FlexLexer* lexer = new yyFlexLexer;

	 //	if ( argc > 0 )
	 //		  yyin = fopen( argv[0], "r" );
	 //else
	 //	  yyin = stdin;
    cout << "Creating temp1.c " << endl;
    int token;
	 int mymethod=0;
    float neq;
    t.lastid=0;
    char myeq[MAXLINE];
    string buf;
    string myequation1;
    ostringstream out(buf);
    init_st(st);
    int i, pointerid;
    // checking if all equations' have initial conditions.
    // build variables list
    string list_var;
    int ninitial_conditions = 0;
	 // while here !
    token = lexer->yylex();
    nomore = 0;
    if (token == INITIAL_CONDITIONS)
      {
	token = lexer->yylex();
	if (token!= IDENTIFIER)
	  error("Expecting an variable in line ", mylineno);
	while (token == IDENTIFIER)
	  {
	    myequation1="";
	    pointerid = lastid;
	    token = lexer->yylex();
	    if (token != ASSIGNMENT) error("Expecting  =  in line ", mylineno);
	    token = lexer->yylex();
	    while ((token == IDENTIFIER) || (token == OPERATOR) || (token == NUMBER)|| (token == FUNCTION) || (token == UNIFORM) || (token == TRIANGULAR)   || (token == CINTERVALTYPE)  || (token == LOGTYPE) )
	      {  if (token == IDENTIFIER)
		{
		  myequation1 = myequation1+st[lastid].name+" ";
		}
	      else
		myequation1 = myequation1+lexer->YYText()+" ";
	      token = lexer->yylex();
	      }
	    if (token != SEMI_COLON) error("Expecting ; in line ", mylineno);
	    st[pointerid].constequation = myequation1;
	    st[pointerid].initialized = 1;
	    ostringstream out(buf);
	    out << pointerid;
	    st[pointerid].constname = "C["+ out.str()+"]";
	    //	cout << pointerid << " : " << st[pointerid].name << ", " << st[pointerid].constname << " " << st[pointerid].constequation << endl;
	    if (st[pointerid].var_type == 0)
	      st[pointerid].var_type = CONST;
	    ninitial_conditions++;
	    mc[codelineno].variable = pointerid;
	    mc[codelineno].linecode = myequation1;
	    codelineno++;
	    token = lexer->yylex();
	  }
      }
    int nequations = 0;
    if (token == EQUATIONS)
      {
	token = lexer->yylex();
	i = 0;
	while (token == IDENTIFIER)
	  {
	    nequations++;
	    pointerid = lastid;
	    myequation1="";
	    token =lexer->yylex();
	    if (token == APOSTROPHE) // Level equation
	      {
		st[pointerid].var_type = LEVEL;
		token = lexer->yylex();
		if (token != ASSIGNMENT)
		  error("Expecting = in line", mylineno);
		token = lexer->yylex();
		//   cout << token << endl;
		while ((token == IDENTIFIER) || (token == OPERATOR) || (token == NUMBER) || (token == FUNCTION))
		  {
		    if (token == IDENTIFIER)
		      myequation1 = myequation1+st[lastid].name+"";
		    else
		      myequation1 = myequation1+lexer->YYText()+"";
		    token = lexer->yylex();
		  }
		if (token != SEMI_COLON) error("Expecting ; in line ", mylineno);
		st[pointerid].dfequation = myequation1;
		if (st[pointerid].eqname == "")
		  {
		    ostringstream out(buf);
		    out <<  nequations;
		    st[pointerid].eqname = "y["+out.str()+"]";
		  }
		else error("You are modifing a derivative equation at line ", mylineno);
	      }
	    else
	      {
		if (token != ASSIGNMENT)
		  error(" Expecting = in line ", mylineno);
		token = lexer->yylex();
		while ((token == IDENTIFIER) || (token == OPERATOR) || (token == NUMBER) || (token == FUNCTION))
		  {
		    if (token == IDENTIFIER)
		      myequation1 = myequation1+st[lastid].name+"";
		    else
		      myequation1 = myequation1+lexer->YYText()+"";
		    token = lexer->yylex();
		  }
		if (token != SEMI_COLON) error(" expecting ; in line ", mylineno);
		st[pointerid].dfequation = myequation1;
		ostringstream out(buf);
		out <<  nequations;
		if (st[pointerid].eqname == "")
		  {
		    st[pointerid].var_type = AUX;
		    st[pointerid].eqname = "AUX["+out.str()+"]";
		  }
		else
		  error("You are modifing an auxiliar equation at line  ", mylineno);
	      }
	    st[pointerid].declaration = 1;
	    mc[codelineno].variable = pointerid;
	    mc[codelineno].linecode = st[pointerid].dfequation;
	    codelineno++;
	    levelvar[nequations] = pointerid;
	    st[pointerid].p_diffeq = nequations;
	    token = lexer->yylex();
	    i++;
	  }
      }
    // check simulation parameters
    double xfromt=1.0, xtot=1.0, xtimestep=1.0, xtimefigs = 1.0;
    if (token == FROM_T)
      {
	token = lexer->yylex();
	if (token != NUMBER) error("Expecting a number at line", mylineno);
	xfromt = lastnumber;
	token = lexer->yylex();
      }
    if (token == TO_T)
      {
	token = lexer-> yylex();
	if (token != NUMBER) error("Expecting a number in to-t ");
	xtot = lastnumber;
	token = lexer->yylex();
      }
    if (xtot < xfromt) error(" error in time interval, fromt should be more than tot ", mylineno);
    if (token == TIMESTEP)
      {
	token = lexer->yylex();
	if (token != NUMBER) error("Expecting a timestep ");
	xtimestep = lastnumber;
	if (xtimestep <= 0.000000001) error ("Timestep too small :", mylineno);
	token = lexer->yylex();
      }
    if (token == TIMEFIGS)
      {
	token = lexer->yylex();
	if (token != NUMBER) error("Expecting a time for figures ");
	xtimefigs = lastnumber;
	if (xtimefigs <= 0.0000001) error("number too small to plot intervals: ", mylineno);
	token = lexer->yylex();
      }
    nomore = 1;
    if (token == PLOT)  // syntax: $display variable ... ;
      {  // time is the first implicit variable
	//  cout << " Plot" << endl;
	tfigs[0] = 0;
	token = lexer->yylex();
	// cout << "Token =" << token << endl;
	while(token == IDENTIFIER)
	  {
	    nfigs++;
	    if (lastid == 0) error("at $Plot: undeclared variable in line ", mylineno);
	    tfigs[nfigs] = lastid;
	    token = lexer->yylex();
	    //cout  << st[nfigs].name <<  " ";
	  }
	if (token != SEMI_COLON) error("at $Plot: Expecting an identifier or ; in line",mylineno);
      }
    token = lexer->yylex();
    // check method
    if (token == METHOD)
      {
	token = lexer->yylex();
	if (token == RungeKutta4th)
	  mymethod = RungeKutta4th;
	else
	  if (token ==  CashKarpRungeKutta)
	    mymethod = CashKarpRungeKutta;
	  else
	    if (token ==  BSStepPoly)
	      mymethod =  BSStepPoly;
	    else
	      if (token == BSStepRational)
			  mymethod = BSStepRational;
	      else
			  if (token == StiffRosenbrock)
				 mymethod =   StiffRosenbrock;
			  else
				 if (token == StiffSemiImplicit)
					mymethod = StiffSemiImplicit;
				 else
					if (token == StiffBSStep)
					  mymethod = StiffBSStep;
					else
					  error("Unknown method at line", mylineno);
      }
    if (ninitial_conditions == 0) error("Please add $Initial_conditions.  ", mylineno);
    if (nequations == 0) error("Please add $equations at line ", mylineno);
    if (mymethod == 0) error("Please include a method ($Rk4, $AdaptiveRk4, $BSStepPoly, or $BSStepRational) ", mylineno);
      // print LHS data into a File
    token = lexer->yylex();
    // random number
    int loc=0; int ii;
    long myneg = 0;
    int lhscounter = 0; int n_sample_lhs =0;
    if (token == LHSRUN)
      {
	//tfigs[0] = 0;
	token = lexer->yylex();
	if (token != NUMBER) error("at $LHS, please input a number (number of samples) at line ", mylineno);
	n_sample_lhs = (int)floor(lastnumber);
	if (n_sample_lhs <= 0) error ("at $LHS Number should be a postive number at line", mylineno);
	if (n_sample_lhs > MAX_LHS_MATRIX) error("at $LHS Number too larger for number of samples/simulations at line",mylineno);
	token = lexer->yylex();
	if (token != IDENTIFIER) error("$LHS, Expecting a parameter defined with a range in line ", mylineno);
	while(token == IDENTIFIER)
	  {
	    if (lastid == 0) error("at $LHS, undeclared variable in line ", mylineno);
	    lhscounter++;
	    //if (lhscounter >=  MAX_LHS_MATRIX) error ("too many $lhs variables ");
	    LHSDecl[lhscounter] = lastid;
	    st[lastid].lhsrun=lhscounter;
            //cout << st[lastid].name << "\t" << st[lastid]. endl;
 //  << st[lastid].begin_range << "\t" << st[lastid].mid_range << " \t" << st[lastid].end_range << "\t" << st[lastid].type_range << endl;
	    token = lexer->yylex();
	  }
	if (token != SEMI_COLON) error("at $LHS, expecting an identifier or ; in line",mylineno);
	token = lexer->yylex();
      }
    if (token == SEED)
      {
	token = lexer->yylex();
	if (token == NUMBER)
	  {
	    myneg = -(int)floor(lastnumber);
	    long *midum = (&myneg);
	    //cout << "Seed = " << myneg << endl;
	  }
	else
	  if (token == CLOCK)
	    {
	      myneg = -time(0);
	      long *midum = &myneg;
	      //cout << "Seed Clock= " << myneg << endl;

	    }
	  else error("Seed: should be a number or $clock ");
	token = lexer->yylex();
      }
    else
      if (token > 0)
	error("Unknown command at line ", mylineno);
    if ((myneg >=0) && (lhscounter > 0)) error("Please include a positive seed number or $seed command ");




    double LHSMatrix[lhscounter+1][MAX_LHS_MATRIX];

    // If OK then read the file and continue using flex
    // here to replace the variables
    // reset
    // build arrays for va

    nomore =1; // no more insertions in symbol table
    int j;
    // create the file temp1.c where the output of the code is:
    //
    ofstream temp1("temp1.c");
    temp1 << "/*  Program generated by teval1, a translator and LHS/PRCC evaluator written by Jose L. Segovia,  2002-03, for the Kirschner Lab.  */ " << endl;
    temp1 << "#include <stdio.h> " << endl;
    temp1 << "#include <stdlib.h> " << endl;
    temp1 << "#include <math.h>  " << endl;
    temp1 << "#include <time.h>  " << endl;
    temp1 << "#include <string.h> " << endl;
    temp1 << "#include " << (char)34 << "datatype.h" << (char)34 << endl;
    temp1 << "#define NRANSI " << endl;
      temp1 << "#include "<<(char)34 << "my_nr.h" << (char)34 << endl;
      temp1 << "#include "<<(char)34 << "nrutil.h" << (char)34 << endl;
      temp1 << "double dxsav,*xp,**py;  /* defining declarations */ " << endl;
      temp1 << "int kmax,kount; " << endl;
      temp1 << "int nrhs;   /* counts function evaluations */ "<< endl;

    temp1 << "#define NINC " << ninitial_conditions+2 << endl;
    temp1 << "#define N " << nequations+2 << endl;
    temp1 << "#define NFigs " << nfigs+2 << endl;
    temp1 << "\t int myreturnsyscall = 0;" << endl;
    temp1 << endl;

 temp1 << "void derivs(double TimeAct, double y[], double dydx[], double C[], double AUX[]) " << endl;
 temp1 << "{" << endl;
 temp1 << "	nrhs++; " << endl;
 for (i=ninitial_conditions; i< codelineno; i++)
    {
      if (st[mc[i].variable].var_type == LEVEL)
	temp1 << "\t" <<  "dydx[" << st[mc[i].variable].p_diffeq << "]=  ";
      else
	temp1 << "\t" << st[mc[i].variable].eqname << "=  ";
      //  cout << i << ": " << mc[i].linecode << endl;
      istringstream StringFlow(mc[i].linecode);
      token = lexer->yylex(&StringFlow, &cout);
      while(token != 0)
	{
	  if ((token == IDENTIFIER))
	    {
	      if (st[lastid].eqname != "")
		temp1 << "" << st[lastid].eqname;
	      else
	       if (st[lastid].constname != "")
	          temp1 << "" <<  st[lastid].constname;
	       else
	          {
		   temp1 << "   ---ERROR IN SOURCE FILE!--- " << endl;
		  temp1.close();
		  error("The variable =>" + st[lastid].name + "<= was not declared! ");
//	          temp1 << " UNDECLARED_VARIABLE_HERE_PLEASE_CHECK_YOUR_SOURCE_FILE ";
		  }
	    }
	  else
	    {
	      temp1 << ""  << lexer->YYText();
	    }
	  token = lexer->yylex();
	}
      temp1 << ";  ";
      temp1 << "\t\t /* " << st[mc[i].variable].name << "=" << mc[i].linecode << "; */" << endl;
	 }
	temp1 << "}" << endl;
	temp1 << "void load_constants(double C[]) {  " << endl;
	temp1 << "\t /* Filling initial conditions */ " << endl;
	int printthis = 0;
	int flagids = 0;
	int my_lastid = 0; // For the constant identifier
	for (i=0; i< codelineno; i++)
	  {
	    if (i < ninitial_conditions)
	      {
	 // temp1 << "/* " << i << " */ " ;
         printthis = 0;
	 temp1 << "\t " << st[mc[i].variable].constname << " =  ";
	// cout << mc[i].linecode <<  " || " << endl;
	 istringstream StringFlow(mc[i].linecode);
	 token = lexer->yylex(&StringFlow, &cout);
	 flagids = 0; // to check just for one token to be assigned to the variable.
	 if  ((token == IDENTIFIER) && (st[lastid].var_type == CONST))
	   temp1 << st[my_lastid].constname << " ";
	 else
	   temp1 << lexer->YYText() << " ";
	 token = lexer->yylex();
	// cout << token << endl;
	 if (token != 0)
	   {
	     if (token == NUMBER)
	       st[mc[i].variable].begin_range = lastnumber;
	     else error("Expecting a number as initial range value of variable: "+st[mc[i].variable].name);
	     token = lexer->yylex();
	     if (token == NUMBER)
	       st[mc[i].variable].end_range = lastnumber;
	     else error("Expecting a number as a range value of variable: "+st[mc[i].variable].name);
	     token = lexer->yylex();
	     if (token == NUMBER)
	       {
		 st[mc[i].variable].mid_range = st[mc[i].variable].end_range;
		 st[mc[i].variable].end_range = lastnumber;
	       }
	     else error("Expecting a number as end range value of variable: "+st[mc[i].variable].name);
	     token = lexer->yylex();
	     // cout << token << endl;
	     if (token == TRIANGULAR)
	       st[mc[i].variable].type_range = "T";
	     else if (token == UNIFORM)
	       st[mc[i].variable].type_range = "U";
	     else if (token == CINTERVALTYPE )
	       st[mc[i].variable].type_range = "C";
	     else if (token == LOGTYPE )
	       st[mc[i].variable].type_range = "L";
	     else
	       error("Unknown range type, should be $U, $T, $C, or $L, at variable :"+st[mc[i].variable].name);
	     if (st[mc[i].variable].begin_range > st[mc[i].variable].end_range)
	       error("Initial range value larger than end range of "+st[mc[i].variable].name);
	     if (st[mc[i].variable].begin_range > st[mc[i].variable].mid_range)
	       error("Initial range value larger than midvalue  of "+st[mc[i].variable].name);
	     if (st[mc[i].variable].mid_range > st[mc[i].variable].end_range)
	       error("Middle range  value larger than end range of "+st[mc[i].variable].name);
	   }
	  // 	cout << i << " : " ;

	//  cout << st[mc[i].variable].name << "  B=" << st[mc[i].variable].begin_range << " M=" << st[mc[i].variable].mid_range << " E=" << st[mc[i].variable].end_range << ", " << st[mc[i].variable].type_range << endl;
	 temp1 << "; ";
	 temp1 << "\t\t /* " << st[mc[i].variable].name << "=" << mc[i].linecode << "; */" <<   endl;
       }
   }
  temp1 << "\t } " << endl;
  temp1 << "/* End of loading constants */ " << endl;

  // Filling auxiliary variables
  // Filling level variables with initial conditions

  temp1 << "void load_auxlevel(double AUX[], double y[], double C[]) {" << endl;
  temp1 << "\t /* Filling leval and aux variables with initial conditions */ " << endl;
  for (i=0; i< codelineno; i++)
    {
      if (i >= ninitial_conditions)
       {
	 if  (st[mc[i].variable].constname > "")
	   temp1 << "\t " << st[mc[i].variable].eqname << " =  " << st[mc[i].variable].constname << "; \t /* " << st[mc[i].variable].name << " = */ " << endl;
	 else
	   if (st[mc[i].variable].var_type == LEVEL)
	     {
	       temp1 << "  ---ERROR IN SOURCE FILE!--- " << endl;
	       temp1.close();
	       error("The variable =>"+ st[mc[i].variable].name +"<= should be initialized.");
	     }
	   else
	     temp1 << "\t " << st[mc[i].variable].eqname << " =  0.0; \t /* " << st[mc[i].variable].name << " = */ " << endl;
       }
    }
  temp1 << "} " << endl;
  temp1 << "\n"  << endl;


  // Here code for lhs evaluation
  // checking if variables have been initalized
  int NumSim = n_sample_lhs;
  for (i=1; i<codelineno; i++)
       if ((st[i].lhsrun >= 1)&&(st[i].type_range == ""))
	 {
	   error("The range of variable: "+st[i].name + " was not initialized.");
	 }
       else
	 if (st[i].lhsrun >= 1)
	   {
	     int mypointer = st[i].lhsrun;
	     double temp;
	     FillHypercube(i, LHSMatrix, st, NumSim);
	     if (lhscounter >= 2)
	     for(ii=0; ii<NumSim; ii++)
	       {
		 loc = rand_int(ii, NumSim-1, &myneg);
		 //cout << LHSMatrix[mypointer][ii] << " == " <<  LHSMatrix[mypointer][loc] << " -> " ;
		 temp = LHSMatrix[mypointer][ii];
		 LHSMatrix[mypointer][ii] =  LHSMatrix[mypointer][loc];
		 LHSMatrix[mypointer][loc] = temp;
		 //cout << LHSMatrix[mypointer][ii] << " == " <<  LHSMatrix[mypointer][loc] << " -> " ;
		 //cout << " Swapping i, ii " << mypointer << ", " << ii <<  ", " << loc << "\t";
	       }
	     //	  cout << endl;
	     /* if (element==LabelPar) PrintLabelSample(element); */
	   }

  // print LHS matrix into a File
  ofstream mylhs("lhsmatrix");
  mylhs <<lhscounter << "\t" << NumSim << "\t" << nfigs << "\t" << (int)(xtot/xtimefigs)+1 << endl;  // nfigs
  for (i=1; i<=lhscounter; i++)
    {
      for (ii=0; ii<NumSim ; ii++)
	mylhs << LHSMatrix[i][ii] << " ";
      mylhs << endl;

    }
   mylhs.close();
  // filling outcome/figs variables
  ofstream lhsoutcome("lhsoutcome");
  lhsoutcome << nfigs << endl;
  for (i=1; i<=nfigs; i++)
	 lhsoutcome << st[tfigs[i]].name << endl;
  lhsoutcome.close();

  temp1 << "void update_variables_lhs(int c, double C[]) { " << endl;
  temp1 << "\t /* Updates the following variables at point given by c */ " << endl;
  temp1 << "\t " << endl;
  int uptoNUMSIM;
  for (uptoNUMSIM = 0; uptoNUMSIM<NumSim; uptoNUMSIM++)
    {
      temp1 << "\t if (c == " << uptoNUMSIM << ")   { " << endl;
      for (i=1; i<codelineno; i++)
	if (st[i].lhsrun >= 1)
	  {
	    int mypointer = st[i].lhsrun;
	    temp1 << " \t\t " << st[i].constname <<  " = " << LHSMatrix[mypointer][uptoNUMSIM] << "; " << endl;
	  }
      temp1 << "\t } " << endl;
      //      temp1 << "\t " << i << ", " << st[i].lhsrun << ", ";
      //      for (int jj=0; jj<NumSim; jj++)
      //	temp1 << LHSMatrix[mypointer][jj] << " ";
      //      temp1 << endl;
    }
  temp1 << " } /*  End of update_variables_lhs  */ " << endl;
  temp1 << endl << endl;

  myreturn = system("rm -f 0*");    if (myreturn != 0) { cout << "Error calling rm -f 0*";   }
  myreturn = system("rm -f 0*gp");  if (myreturn != 0) { cout << "Error calling rm -f 0*gp"; }

  // function to load figures data:
  char filename[15];
  temp1 << "static int nfigsloaded; " << endl;
  temp1 << "int load_data_figures(double FIGS[][MAX_FIGS], eqfig figdata[], int nfigs, int jj, double y[], double C[], double AUX[], double x, FILE * myfile) { " << endl;
  temp1 << "\t int j = nfigsloaded; " << endl;
  temp1 << "\t int i; int errorover=0; " << endl;
  // temp1 << "\t int myreturnsyscall = 0;" << endl;
  temp1 << "\t\t  FIGS[0][j] = x; " << endl;
  for (i = 1; i<=nfigs; i++)
	 {
	   if (st[tfigs[i]].eqname <= "") { temp1.close(); myreturn=system("rm -f temp1.c"); error("You cannot plot a constant/parameter: "+ st[tfigs[i]].name); }
      temp1 << "\t\t  if (! finite(" <<  st[tfigs[i]].eqname << ")) { printf(" <<(char)34 << "Overflow error on variable "<<st[tfigs[i]].name<<" \\n" << (char)34 << "); errorover++; } " << endl;
      temp1 << "\t\t  FIGS["<<i<<"][j] = "<< st[tfigs[i]].eqname << ";" << endl;
	   temp1 << "\t\t  if (figdata["<<i<<"].max < FIGS["<<i<<"][j]) { figdata["<<i<<"].max =  FIGS["<<i<<"][j]; } "<< endl;
	   temp1 << "\t\t  if (figdata["<<i<<"].min > FIGS["<<i<<"][j]) { figdata["<<i<<"].min =  FIGS["<<i<<"][j]; } "<< endl;
	   if (NumSim > 0)
	     {
	       sprintf(filename, "%04u.gp",i);
	       ofstream plotlhs(filename);
	       plotlhs.width(3);
	       plotlhs << setiosflags(ios::scientific);
	       //plotlhs << setiosflags(ios::fixed);
	       //plotlhs << setiosflags(ios::showpoint);
	       plotlhs << "plot ";
		    // legend // lhsdata
	       ofstream  lhsdata("lhsdata");
	       lhsdata << lhscounter << endl;
	       ofstream  filelegend ("legend.gp");
//	       filelegend << "set key below  Left title \'\t"; // box 4 \n";
//	       filelegend << "set key right outside Right title \'LHS Matrix " ;
          if (lhscounter <= 1)
			 filelegend << "set key right outside Right title \'";
			 else
	       filelegend << "set nokey" << endl;
	       for (int uptoNUMSIM = 0; uptoNUMSIM<NumSim; uptoNUMSIM++)
						 {
						 sprintf(filename, "%04u",uptoNUMSIM+1);
						 plotlhs << "'" <<filename << "' u 1:" << i+1 << " t '" << uptoNUMSIM+1 << "\t";
						  plotlhs << setiosflags(ios::scientific);
						  //plotlhs << setiosflags(ios::fixed);
						  //plotlhs << setiosflags(ios::showpoint);
						 if (lhscounter>=1)
						 for (int xlhs = 1; xlhs <= lhscounter; xlhs++)
							{
							  if (lhscounter <=1 ) {
							     plotlhs <<  LHSMatrix[xlhs][uptoNUMSIM] << "\t";
                                 // cout <<  xlhs << " --- "<< LHSMatrix[xlhs][uptoNUMSIM] << "\t";
							      }
							  if (uptoNUMSIM == 0)
								 {
								 if (lhscounter <=1)
							          filelegend << st[LHSDecl[xlhs]].name << "\'" << endl;
					lhsdata << st[LHSDecl[xlhs]].name  << "\t" <<  st[LHSDecl[xlhs]].begin_range << "\t" ;
                    lhsdata << st[LHSDecl[xlhs]].mid_range << " \t" << st[LHSDecl[xlhs]].end_range << "\t" << st[LHSDecl[xlhs]].type_range << endl;
                    			 }
							}
						 plotlhs << resetiosflags(ios::scientific);
						 //plotlhs << resetiosflags(ios::fixed);
						 //plotlhs << resetiosflags(ios::showpoint);
						 plotlhs << "' w p ";
						 if (uptoNUMSIM < (NumSim-1)) plotlhs << ", ";
						 }
	       plotlhs << endl;
	       plotlhs.close();
//	       filelegend << "\' box 4 \n " << endl;
	       filelegend.close();
	       lhsdata.close();

	     }
	 }
  temp1 << "\t nfigsloaded ++;" << endl;
  temp1 << "\t if (nfigsloaded >= MAX_FIGS ) { int xx, yy; " << endl;
  temp1 << "\t\t for (xx=0;xx<MAX_FIGS;xx++) { " << endl;
  temp1 << "\t\t fprintf(myfile," << (char)34 << "%f\\t"<< (char)34 << ",FIGS[0][xx]); " << endl;
  temp1 << "\t\t  for (yy=1;yy<=" <<nfigs <<";yy++) " << endl;
  temp1 << "\t\t\t  fprintf(myfile, " << (char)34 << "%.6e\\t" << (char)34 << ",FIGS[yy][xx]); " << endl;
  temp1 << "\t\t fprintf(myfile,"<<(char)34<<"\\"<< "n" << (char)34 << "); " << endl;
  temp1 << "} nfigsloaded =0; } " << endl;
  temp1 << " return errorover; " << endl;
  temp1 << "} " << endl << endl << endl;
  if (NumSim <= 0) NumSim =1;
  // End of code for lhs evaluation
  temp1 << "int main (void) " << endl;
  temp1 << "{ " << endl;
  temp1 << "\t time_t init_time=time(NULL), end_time; double mydifftime;" << endl;
  // temp1 << "\t int borra = gettimeofday(init_time); " << endl;
  // from here functions
  if (mymethod == RungeKutta4th)
	 {
  temp1 << "\t int i,j, resode=0;  int myc=0; " << endl;
  temp1 << "\t double h, x=1.0, *y, *dydx, *yout, *C, *AUX; " << endl;
  temp1 << "\t FILE *stream; " << endl;
  temp1 << "\t char filename[5]; " <<endl;
  temp1 << "\t int nfigs = " << nfigs << ";" <<endl;
  temp1 << "\t eqfig figdata["<<nfigs+2<<"];" << endl;
  int nrowsfigs= (int)((xtot-xfromt)/xtimefigs);
  nrowsfigs = (int)(nrowsfigs*1.2);
  temp1 << "\t double FIGS[" <<  nfigs+2 << "][MAX_FIGS]; " << endl;
  temp1 << "\t int maxloop =0, nloop =0; int figres=0; " << endl;
  temp1 << "\t y    =  dvector(1,N);  " << endl;
  temp1 << "\t dydx =  dvector(1,N); " << endl;
  temp1 << "\t yout =  dvector(1,N);  " << endl;
  temp1 << "\t C   =   dvector(1,NINC); " << endl;
  temp1 << "\t AUX =   dvector(1,N); " <<endl;
  temp1 << "\t nfigsloaded = 0; " << endl;
  temp1 << "\t for(i=1;i<=nfigs;i++) { figdata[i].max = -HUGE_VAL; figdata[i].min=HUGE_VAL; }; " << endl;
  for (i=1; i<=nfigs; i++)
	 temp1 << "\t  strcpy(figdata[" << i << "].name," << (char)34 << st[tfigs[i]].name << (char)34 << ");" << endl;

  ////
  temp1 << "\t while (myc < " <<  NumSim << ") { " << endl;
  temp1 << "\t load_constants(C); " << endl;
  temp1 << "\t update_variables_lhs(myc, C); " << endl;
  temp1 << "\t load_auxlevel(AUX,y,C); " << endl;
   temp1 << "\t nfigsloaded = 0; " << endl;
  temp1 << "\t /* ---rk4--- */ " << endl;
  temp1 << "\t h = (double)" << xtimestep<< ";"  << endl;
  temp1 << "\t x = (double)" << xfromt << "; " << endl;
  temp1 << "\t maxloop = (int) (" << xtot << "/h)+1;" << endl;
  temp1 << "\t nloop = " << xfromt << "; " << endl;
  temp1 << "\t j=0; "<< endl;
  temp1 << "\t printf(" << (char)34 << "Solving run %d " << "\\" << "n" << (char)34 << ", myc+1);" << endl;
  temp1 << "\t sprintf(filename, " << (char)34 << "%04u" << (char)34 << ", myc+1); " << endl;
  temp1 << "\t stream = fopen(filename," <<   char(34) << "w" <<  char(34) << "  ); " << endl;

//  temp1 << "\t figres=load_data_figures(FIGS, figdata, nfigs, j, y, C, AUX, x, stream);" << endl;
//  temp1 << "\t j++; " << endl;
  // here maxloop = (int)(xtot/xtimesep)+1
  temp1 << "\t while(nloop <= maxloop) {" << endl;
//
  temp1 << "\t\t if (nloop % (int)ceil("<<xtimefigs << "/h) == 0) {" << endl;
  temp1 << "\t\t resode=load_data_figures(FIGS, figdata, nfigs, j, y, C, AUX, x,  stream);" << endl;
  temp1 << "\t\t if (resode > 0) break; " << endl;
  temp1 << "\t\t\t j++; }  " << endl;
//
  temp1 << "\t\t derivs(x,y,dydx,C, AUX); " << endl;
  temp1 << "\t\t rk4(y,dydx,N,x,h,yout,derivs, C, AUX); " <<endl;
  temp1 << "\t\t for(i=0; i<N; i++) " << endl;
  temp1 << "\t\t y[i] = yout[i]; " << endl;

  //

 //
  temp1 << "\t\t nloop++; " << endl;
  temp1 << "\t\t x = x + h;" << endl;
   temp1 << "\t }  " << endl;
  temp1 << "\tfor(i=0; i<nfigsloaded; i++) { " << endl;
  temp1 << "\t\t fprintf(stream,"<< (char)34 <<"%f\\t"<< (char)34 << ",FIGS[0][i]);" << endl;
  temp1 << "\t\t for(j=1; j<=nfigs; j++) fprintf(stream," << (char)34 << "%.5e\\t" << (char)34 << ",FIGS[j][i]); " << endl;
  temp1 << "\t\t fprintf(stream,"<<(char)34<<"\\"<< "n" << (char)34 << "); }; " << endl;
  temp1 << "\t\t fclose(stream); " << endl;
 // temp1 << " } " << endl;
  temp1 << "\t if (resode > 0) { printf( " << (char)34 << "Error solving run %d "<< "\\" << "n" << (char)34 << ", myc+1); break;	} " << endl;
  temp1 << "\t myc++; } " << endl;
  temp1 << "\t WritePlotScript(1, nfigs, " <<nrowsfigs<<", figdata, 0, myc-1);" << endl;
  temp1 << "\t myreturnsyscall = system(" << char(34) << "rm output1.pdf; gnuplot plot.scr &" << char(34) << "); " << endl;
  temp1 << "\t free_dvector(yout,1,N); " << endl;
  temp1 << "\t free_dvector(dydx,1,N); " << endl;
  temp1 << "\t free_dvector(y,1,N); " << endl;
  temp1 << "\t free_dvector(C,1,NINC); " << endl;
  temp1 << "\t free_dvector(AUX,1,N); " << endl;
//  temp1 << "\t free_dmatrix(FIGS, 0," <<  nrowsfigs+1 << ", 0, MAX_FIGS);" << endl;
	 }
  else
	 {
		if ((mymethod == CashKarpRungeKutta) || (mymethod == BSStepPoly) || (mymethod == BSStepRational) ||
			 (mymethod == StiffRosenbrock ) || (mymethod == StiffBSStep))
	 {
		temp1 << "\t int i,j,nbad,nok, resode; int myc=0; " << endl;
		temp1 << "\t double eps=1.0e-8,h1="<<xtimestep<<",hmin=0.0,x1=" << xfromt << "; " << endl;
		temp1 << "\t double x2=" << xtot <<" ,*y, *C, *AUX; " << endl;
		temp1 << "\t FILE *stream; " << endl;
		temp1 << "\t char filename[5]; " << endl;
		temp1 << "\t int nfigs = " << nfigs << ";" <<endl;
		int nrowsfigs= (int)((xtot-xfromt)/xtimestep);
		nrowsfigs = (int)(nrowsfigs*1.2);
		temp1 << "\t double FIGS["<< nfigs+2 << "][MAX_FIGS]; " << endl;
		temp1 << "\t eqfig figdata["<<nfigs+2<<"];" << endl;
		temp1 << "\t int NPOINTS = " << nrowsfigs << "; " << endl;
		temp1 << "\t y=  dvector(1,N); " << endl;
		temp1 << "\t xp= dvector(1,NPOINTS); " << endl;
		temp1 << "\t C = dvector(1,NINC);  " << endl;
		temp1 << "\t AUX=dvector(1,N);  " << endl;
		temp1 << "\t py=dmatrix(1,N,1,NPOINTS);   " << endl;

		temp1 << "\t for(i=1;i<=nfigs;i++) { figdata[i].max = -HUGE_VAL; figdata[i].min=HUGE_VAL; }; " << endl;
		for (i=1;i<=nfigs; i++)
		  temp1 << "\t  strcpy(figdata[" << i << "].name," << (char)34 << st[tfigs[i]].name << (char)34 << ");" << endl;
		temp1 << "\t while (myc < " <<  NumSim << ") { " << endl;
		temp1 << "\t eps=1.0e-8; h1="<<xtimestep<<"; hmin=0.0; x1=" << xfromt << "; " << endl;
		temp1 << "\t x2=" << xtot << ";" << endl;
		temp1 << "\t load_constants(C); " << endl;
		temp1 << "\t update_variables_lhs(myc, C); " << endl;
		temp1 << "\t load_auxlevel(AUX,y,C); " << endl;
		temp1 << "\t nfigsloaded = 0; " << endl;
                temp1 << "\t nrhs=0.0;  " << endl;
		temp1 << "\t kmax=NPOINTS-2;  " << endl;
		temp1 << "\t dxsav="<<xtimestep<<"; " << endl;
		temp1 << "\t\t printf(" << (char)34 << "Solving run %d " << "\\" << "n" << (char)34 << ", myc+1);" << endl;
		temp1 << "\t\t  sprintf(filename, " << (char)34 << "%04u" << (char)34 << ", myc+1); " << endl;
		temp1 << "\t\t stream = fopen(filename," <<   char(34) << "w" <<  char(34) << "  ); " << endl;
	 }
  if (mymethod == CashKarpRungeKutta)
	  temp1 << "\t resode=odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs, rkqs, C, AUX, FIGS, figdata, nfigs, stream); " << endl;
  if (mymethod ==  BSStepPoly)
	   temp1 << "\t resode=odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep, C, AUX, FIGS, figdata, nfigs, stream); " << endl;
  if (mymethod ==  BSStepRational)  // seems not a good method
	   temp1 << "\t resode=odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,bsstep_rz, C, AUX, FIGS, figdata, nfigs, stream); " << endl;
  if (mymethod == StiffRosenbrock)
      temp1 << "\t resode=odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,stiff, C, AUX, FIGS, figdata, nfigs, stream); " << endl;
  if (mymethod == StiffBSStep)
	 temp1 << "\t resode=  odeint(y,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,stifbs, C, AUX, FIGS, figdata, nfigs, stream); " << endl;
  if ((mymethod == CashKarpRungeKutta) || (mymethod == BSStepPoly)  || (mymethod == BSStepRational) || (mymethod==StiffRosenbrock) || (mymethod == StiffBSStep))
	 {

	   temp1 << "\t\t for(i=0; i<nfigsloaded; i++) { " << endl;
	   temp1 << "\t\t fprintf(stream,"<< (char)34 <<"%f\\t"<< (char)34 << ",FIGS[0][i]);" << endl;
	   temp1 << "\t\t for(j=1; j<=nfigs; j++) fprintf(stream," << (char)34 << "%.5e\\t" << (char)34 << ",FIGS[j][i]); " << endl;
	   temp1 << "\t\t fprintf(stream,"<<(char)34<<"\\"<< "n" << (char)34 << "); }; " << endl;

    temp1 << "\t\t fclose(stream); " << endl;

//    temp1 << "} " << endl;
    temp1 << "\t if (resode > 0) { printf( " << (char)34 << "Error solving run %d "<< "\\" << "n" << (char)34 << ", myc+1);	 break; } " << endl;
    temp1 << "\t myc++; } " << endl;
    temp1 << "\t  WritePlotScript1(1, "<<nfigs<<", kount, figdata, 0, myc-1); " << endl;
    temp1 << "\t myreturnsyscall = system(" << char(34) << "rm output1.pdf; gnuplot plot.scr &" << char(34) << "); " << endl;

//    temp1 << "\t if (myreturnsyscall != 0) { printf ( " << (char)34 << " Error in calling rm output1.ps gnuplot plot.scr" << (char)34 <<   ); << endl;

    temp1 << "\t free_dmatrix(py,1,N,1,NPOINTS);  " << endl;
    temp1 << "\t free_dvector(xp,1,NPOINTS); " << endl;
    temp1 << "\t free_dvector(y,1,N); " << endl;
    temp1 << "\t free_dvector(C,1,NINC); " << endl;
    temp1 << "\t free_dvector(AUX,1,N); " << endl;
  //  temp1 << "\t free_dmatrix(FIGS,0,MPOINTS,0,MAX_FIGS);" << endl;
  };
	 }
  temp1 << "\t end_time = time(NULL); " << endl;
  temp1 << "\t mydifftime = difftime(end_time, init_time); " << endl;
//  temp1 << "\t printf(" << (char)34 << "\\" << "n" << (char)34 << ");" << endl;
  temp1 << "\t printf("<<(char)34 <<"##  %f seconds \\" << "n" << (char)34 << ",mydifftime); " << endl;
  temp1 << "\t return 0; " << endl;
  temp1 << "} " << endl;
  temp1.close();
//  cout << "Compiling .. " << endl;
  return 0;
}


