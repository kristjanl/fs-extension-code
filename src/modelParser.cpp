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
  logger.disable();
  yyparse();
	//simpleImplMain();
  //compMain();
	
	return 0;
}



