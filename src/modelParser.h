/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef MODELPARSER_H_
#define MODELPARSER_H_

#include "Hybrid.h"
#include "SimpleImplApp.h"
#include "CompApp.h"
#include "MyLogger.h"
#include "ParsingUtils.h"
#include "Utils.h"
#include "OutputWriter.h"

class MySettings2;

//extern int lineNum;

extern mpfr_prec_t intervalNumPrecision;

extern ContinuousReachability continuousProblem;
extern HybridReachability hybridProblem;
extern MySettings2 *settings2;
extern bool usePlainFlowstar;

//extern ParseSetting parseSetting;
//extern ParseResult parseResult;

//extern int yyparse();

//void parseError(const char *str, int lnum);

extern ReachabilitySetting mode_local_setting;

#endif /* MODELPARSER_H_ */
