/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "modelParser.h"

/*
int lineNum = 1;

ContinuousReachability continuousProblem;
HybridReachability hybridProblem;
ReachabilitySetting mode_local_setting;

extern int yyparse();

void parseError(const char *str, int lnum)
{
	cerr << "Error @line " << lineNum << ":" << string(str) << endl;
	exit(1);
}*/

int main(int argc, const char *argv[])
{
  
  //compare two logs
  //first 2 arguments are input file names
  //3rd argument is the name of the TM that will be compared
  //if no third argument then it will print the names of the first 2
  //(optional) 4th is the index of the TM in the first file (0 if not present)
  //(optional) 5th is the index in the second file (same as 4th if no present)
  if(argc >= 3) {
    mforce1("print");
    
    int index1 = 0;
    int index2 = 0;
    if (argc >= 5) {
      index1 = atoi(argv[4]);
      index2 = index1;
    }
    if (argc >= 6) {
      index2 = atoi(argv[5]);
    }
    if (argc == 3) {//just print names
      printTMVFiles(argv[1], argv[2], "", -1, -1);
    } else {
      printTMVFiles(argv[1], argv[2], argv[3], index1, index2);
    }
    exit(0);
  }
  
  //given a single argument will convert the TMs in that file into Mathematica's
  //format
  if(argc == 2) {
    //transform to mathematica
    toMathematica(argv[1]);
    exit(0); 
  }
  mdisable();
  yyparse();
	//simpleImplMain();
  //compMain();
	
	return 0;
}



