#ifndef PARSINGUTILS_H_
#define PARSINGUTILS_H_

#include "include.h"
#ifdef DInt
  #include "DoubleInterval.h"
#else
  #include "Interval.h"
#endif
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


#endif /* PARSINGUTILS_H_ */
