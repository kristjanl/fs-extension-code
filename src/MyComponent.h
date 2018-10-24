#ifndef MYCOMPONENT_H_
#define MYCOMPONENT_H_

#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <algorithm>

#include "TaylorModel.h"
#include "MyLogger.h"

class CompDependency;
class PrecondModel;

class MySettings;
class MySettings2;

class MyComponent {
  public:
    MyComponent();
    MyComponent(vector<int> vs, vector<int> tps);
    
    //indexes to be solved in this component
    vector<int> varIndexes;
    //linking variables, probably not needed
    vector<int> linkVars;
    //all variable indexes in component (solve + links)
    vector<int> compVars;
    //indexes of variables that need to be solved (wrt component variable position)
    vector<int> solveIndexes;
    
    //variables introduced in component and dependencies (including implicit)
    vector<int> allVars;
    bool isSolved;
    bool isPrepared;
    bool isPreconditioned;
    bool firstPrecondition;
    bool usingPreconditioning;
    static int nextFreeParam;
    
    //vector of variables that will introduce parameter
    vector<int> varsToBeIntroduced;
    
    //parameters that are introduced in this component
    vector<int> tpIndexes; //TODO rename
    //parameters that are used in flowpipes (component's own + dependencies)
    //ordered by parameter order in the whole system
    vector<int> allTMParams;
    vector<CompDependency *> dependencies;
    
    vector<MyComponent> previous;
    
    //first are component's variables, then come dependencies
    TaylorModelVec timeStepPipe;
    TaylorModelVec initSet;
    TaylorModelVec swInput;
    
    //the part of right taylor model that corresponds to this components params
    TaylorModelVec unpairedRight;
    
    
    
    vector<HornerForm> odes;
    vector<Interval> dom;
    vector<TaylorModelVec> pipes;
    vector<PrecondModel *> pipePairs;
    vector<TaylorModelVec> output;
    
    
    void addDependency(int linkVar, MyComponent *pComp);
    
    vector< vector<int> > previousMappers();
    vector< vector<int> > previousMappers2();
    void log();
    
    void addVar(int var);
    void prepareComponent(TaylorModelVec init, const vector<HornerForm> & ode, 
        vector<Interval> domain, bool discardEmptyParams);
    void prepareVariables(TaylorModelVec tmv, const vector<HornerForm> & ode, 
        bool discardEmptyParams);
    void prepareMappers();
    void remapTimeStepPipe();
    TaylorModelVec orderedTSPRemap(bool first);
    
    bool isSolveVar(int var);
    bool belongsToComp(int param);
    
    int getIntergrationParamCount();
    
    PrecondModel *lastPre();
    
    TaylorModelVec lastPipe();
    
    TaylorModel getRightModelForVar(vector<int> mapper, int var);
    void getIthPipePair(vector<int> lMapper, vector<int> rMapper, 
        TaylorModel & left, TaylorModel & right, int var, int i);
    
    void serializeFlows();
    void deserializeFlows();
    
    void computeMappingPositions(int variable, int *depPos, int *dLinkPos, 
        int *linkPos);
    
  private:
    void remapIVP(TaylorModelVec tmv, const vector<HornerForm> & ode, 
        vector<Interval> domain);
};

class CompDependency {
  public:
    CompDependency(int link, MyComponent *pComp);
    int linkVar;
    MyComponent *pComp;
    vector<int> mapper;
    vector<int> rightMapper;
    vector<int> leftMapper;
};


vector<int> concateMapper(vector<int> & smaller, vector<int> & bigger);

//TODO REMOVE after completing refactoring
vector<MyComponent *> createComponents(MySettings *settings, 
    const vector<HornerForm> & ode);
void makeCompIndexes(MySettings *settings, const vector<HornerForm> & ode);

vector<MyComponent *> createComponents(MySettings2 *settings, 
    const vector<HornerForm> & ode);
void makeCompIndexes(MySettings2 *settings, const vector<HornerForm> & ode);


void prepareComponents(vector<MyComponent *> & comps, TaylorModelVec init, 
    const vector<HornerForm> & ode, vector<Interval> domain, 
    bool discardEmptyParams);

//TODO remove after refactoring
MyComponent getSystemComponent(vector<MyComponent *> comps, 
    TaylorModelVec init, const vector<HornerForm> & ode, 
    vector<Interval> domain, bool discardEmptyParams);

MyComponent* pGetSystemComponent(vector<MyComponent *> comps, 
    TaylorModelVec init, const vector<HornerForm> & ode, 
    vector<Interval> domain, bool discardEmptyParams);

#endif /* MYCOMPONENT_H_ */
