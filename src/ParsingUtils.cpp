#include "ParsingUtils.h"
#include "Utils2.h"

int lineNum = 1;

ContinuousReachability continuousProblem;
HybridReachability hybridProblem;
ReachabilitySetting mode_local_setting;

MySettings *settings2 = new MySettings();

extern int yyparse();

void parseError(const char *str, int lnum) {
	cerr << "Error @line " << lineNum << ":" << string(str) << endl;
	exit(1);
}
