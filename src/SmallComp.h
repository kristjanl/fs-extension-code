#ifndef SMALLCOMP_H_
#define SMALLCOMP_H_

#include <typeinfo>
#include <ctime>

#include "Continuous.h"
#include "MyLogger.h"
#include "OutputWriter.h"
#include "MyComponent.h"
#include "Matrix.h"
#include "Exceptions.h"
#include "Utils.h"
#include "Transformer.h"

namespace smallComp {
  void bar12();
  double shrinkWrap(MyComponent & component, vector<Interval> domain, 
      vector<Interval> step_end_exp_table);
  void shrinkWrapSet(MyComponent & all, MyComponent * component, double factor, 
        vector<Interval> domain);
        
  double applyShrinkWrapping(MyComponent & all, vector<Interval> domain, 
      vector<Interval> step_end_exp_table, vector<MyComponent *> comps,
      OutputWriter & writer);
}

class SmallCompSystem: public ContinuousSystem  {
  public:
    ShrinkWrappingCondition *swChecker;
    SmallCompSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input);
    SmallCompSystem(const ContinuousSystem & system, vector< vector<int> > components);
    void my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const;
    void my_reach_picard_old(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const;
    vector< vector<int> > components;
  private:
    void integrateComponent(const vector<int> comp, vector<TaylorModelVec> & pipes, 
    const OutputWriter writer, const vector<HornerForm> & ode, 
    const TaylorModelVec & init, const vector<Interval> & domain, 
    int order, double step, double time, vector<Interval> step_end_exp_table) const;
    //void addEmptyTM(vector<TaylorModelVec> & pipes) const;
};

class SmallCompReachability: public ContinuousReachability {
  public:
    SmallCompReachability();
    void myRun();
};


#endif /* SMALLCOMP_H_ */
