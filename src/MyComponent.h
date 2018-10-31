#ifndef MYCOMPONENT_H_
#define MYCOMPONENT_H_

#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <algorithm>
#include <vector>

#include "TaylorModel.h"
#include "MyLogger.h"

class CompDependency;
class PrecondModel;

class MySettings;

class MyComponent {
  public:
    MyComponent();
    MyComponent(std::vector<int> vs, std::vector<int> tps);
    
    //indexes to be solved in this component
    std::vector<int> varIndexes;
    //linking variables, probably not needed
    std::vector<int> linkVars;
    //all variable indexes in component (solve + links)
    std::vector<int> compVars;
    //indexes of variables that need to be solved (wrt component variable position)
    std::vector<int> solveIndexes;
    
    //variables introduced in component and dependencies (including implicit)
    std::vector<int> allVars;
    bool isSolved;
    bool isPrepared;
    bool isPreconditioned;
    bool firstPrecondition;
    bool usingPreconditioning;
    
    //vector of variables that will introduce parameter
    std::vector<int> varsToBeIntroduced;
    
    //parameters that are introduced in this component
    std::vector<int> tpIndexes; //TODO rename
    //parameters that are used in flowpipes (component's own + dependencies)
    //ordered by parameter order in the whole system
    std::vector<int> allTMParams;
    std::vector<CompDependency *> dependencies;
    
    std::vector<MyComponent> previous;
    
    //first are component's variables, then come dependencies
    TaylorModelVec timeStepPipe;
    TaylorModelVec initSet;
    TaylorModelVec swInput;
    
    //the part of right taylor model that corresponds to this components params
    TaylorModelVec unpairedRight;
    
    
    
    std::vector<HornerForm> odes;
    std::vector<Interval> dom;
    std::vector<TaylorModelVec> pipes;
    std::vector<PrecondModel *> pipePairs;
    std::vector<TaylorModelVec> output;
    
    
    void addDependency(int linkVar, MyComponent *pComp);
    
    std::vector< std::vector<int> > previousMappers();
    std::vector< std::vector<int> > previousMappers2();
    void log();
    
    void addVar(int var);
    void prepareComponent(TaylorModelVec init, const std::vector<HornerForm> & ode, 
        MySettings *settings);
    void prepareVariables(TaylorModelVec tmv, const std::vector<HornerForm> & ode, 
        bool discardEmptyParams);
    void prepareMappers();
    void remapTimeStepPipe();
    TaylorModelVec orderedTSPRemap(bool first);
    
    bool isSolveVar(int var);
    bool belongsToComp(int param);
    
    int getIntergrationParamCount();
    
    PrecondModel *lastPre();
    
    TaylorModelVec lastPipe();
    
    TaylorModel getRightModelForVar(std::vector<int> mapper, int var);
    void getIthPipePair(std::vector<int> lMapper, std::vector<int> rMapper, 
        TaylorModel & left, TaylorModel & right, int var, int i);
    
    void serializeFlows();
    void deserializeFlows();
    
    void computeMappingPositions(int variable, int *depPos, int *dLinkPos, 
        int *linkPos);
    string getVarName(MySettings *settings);
    
  private:
    void remapIVP(TaylorModelVec tmv, const std::vector<HornerForm> & ode, 
        std::vector<Interval> domain);
};

class CompDependency {
  public:
    CompDependency(int link, MyComponent *pComp);
    int linkVar;
    MyComponent *pComp;
    std::vector<int> mapper;
    std::vector<int> rightMapper;
    std::vector<int> leftMapper;
};


std::vector<int> concateMapper(std::vector<int> & smaller, std::vector<int> & bigger);

//TODO REMOVE after completing refactoring
std::vector<MyComponent *> createComponents(MySettings *settings, 
    const std::vector<HornerForm> & ode);
void makeCompIndexes(MySettings *settings, const std::vector<HornerForm> & ode);

MyComponent* pGetSystemComponent(std::vector<MyComponent *> comps, 
    TaylorModelVec init, const std::vector<HornerForm> & ode, 
    MySettings *settings);

#endif /* MYCOMPONENT_H_ */
