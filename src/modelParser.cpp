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
  if(argc >= 3) {
    logger.force("print");
    int index1 = -1;
    int index2 = -1;
    if(argc == 4) {
      index1 = index2 = atoi(argv[3]);
    } else if (argc == 5) {
      index1 = atoi(argv[3]);
      index2 = atoi(argv[4]);
    }
    logger.log(index1);
    logger.log(index2);
    printTMVFiles(argv[1], argv[2], index1, index2);
    exit(0);
  }
  logger.disable();
  yyparse();
	//simpleImplMain();
  //compMain();
	
	return 0;
}



