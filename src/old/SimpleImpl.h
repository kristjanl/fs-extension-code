#ifndef SIMPLEIMPL_H_
#define SIMPLEIMPL_H_

#include "Continuous.h"
#include "MyLogger.h"
#include "OutputWriter.h"


class SimpleImplSystem: public ContinuousSystem  {
  public:
    SimpleImplSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input);
    SimpleImplSystem(const ContinuousSystem & system);
    void my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const;
};

class SimpleImplReachability: public ContinuousReachability {
  public:
    SimpleImplReachability();
    void myRun();
};

#endif /* SIMPLEIMPL_H_ */
