#ifndef SIMPLECOMP_H_
#define SIMPLECOMP_H_

#include <typeinfo>

#include "Continuous.h"
#include "MyLogger.h"
#include "OutputWriter.h"



class SimpleCompSystem: public ContinuousSystem  {
  public:
    SimpleCompSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input);
    SimpleCompSystem(const ContinuousSystem & system);
    void my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const;
  private:
    void integrateComponent(const vector<int> comp, vector<TaylorModelVec> & pipes, 
    const OutputWriter writer, const vector<HornerForm> & ode, 
    const TaylorModelVec & init, const vector<Interval> & domain, 
    int order, double step, double time, vector<Interval> step_end_exp_table) const;
    void addEmptyTM(vector<TaylorModelVec> & pipes) const;
};

class SimpleCompReachability: public ContinuousReachability {
  public:
    SimpleCompReachability();
    void myRun();
};

#endif /* SIMPLECOMP_H_ */
