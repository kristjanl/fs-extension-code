#ifndef UNITTESTING_H_
#define UNITTESTING_H_

#include "include.h"
#include "Interval.h"
#include "Monomial.h"
#include "Polynomial.h"
#include "MyComponent.h"
#include "OutputWriter.h"
#include "Continuous.h"
#include "Hybrid.h"


extern int lineNum;

extern ParseSetting parseSetting;
extern ParseResult parseResult;

extern int yyparse();

void parseError(const char *str, int lnum);


#endif /* UNITTESTING_H_ */
