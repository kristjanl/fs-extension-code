#include <iostream>

#include "include.h"
#include "Continuous.h"

#include "MyLogger.h"
#include "ExtractedPicard.h"
//#include "SimpleImpl.h"
#include "SimpleComp.h"
#include "CompApp.h"
#include <typeinfo>

using namespace std;

namespace compApp {
  //x1' = -0.1*x1
  Polynomial* dX1(ContinuousReachability continuousProblem) {
    int numVars1 = continuousProblem.stateVarNames.size()+1;
    Interval coefI(-0.1);
    Polynomial *coef = new Polynomial(coefI, numVars1);
    
    vector<int> degrees;
    for(int i=0; i<numVars1; ++i) {
      degrees.push_back(0);
    }
    
    Interval ONE(1);
    Interval ZERO;
    
    degrees[continuousProblem.getIDForStateVar("x1")+1] = 1;
    Monomial monomial(ONE, degrees);

    Polynomial *p = new Polynomial(monomial);
    *p *= (*coef);
    return p;
  }

  //x2' = 0.1*x1 - 0.1*x2
  Polynomial* dX2(ContinuousReachability continuousProblem) {
    int numVars1 = continuousProblem.stateVarNames.size()+1;
    Interval ONE(1);
    Interval ZERO;
    
    Polynomial *coef1 = new Polynomial(Interval(0.1), numVars1);
    vector<int> degrees1;
    for(int i=0; i<numVars1; ++i) {
      degrees1.push_back(0);
    }
    
    degrees1[continuousProblem.getIDForStateVar("x1")+1] = 1;
    Monomial monomial1(ONE, degrees1);

    Polynomial *p1 = new Polynomial(monomial1);
    *p1 *= (*coef1);
    
    Polynomial *coef2 = new Polynomial(Interval(-0.1), numVars1);
    vector<int> degrees2;
    for(int i=0; i<numVars1; ++i) {
      degrees2.push_back(0);
    }
    
    degrees2[continuousProblem.getIDForStateVar("x2")+1] = 1;
    Monomial monomial2(ONE, degrees2);

    Polynomial *p2 = new Polynomial(monomial2);
    *p2 *= (*coef2);
    
    Polynomial *ret = new Polynomial(); //*p1 + *p2;
    *ret = *p1 + *p2;
    return ret;
  }

  //x3' = 0.1*x2
  Polynomial* dX3(ContinuousReachability continuousProblem) {
    int numVars1 = continuousProblem.stateVarNames.size()+1;
    Interval coefI(0.1);
    Polynomial *coef = new Polynomial(coefI, numVars1);
    
    vector<int> degrees;
    for(int i=0; i<numVars1; ++i) {
      degrees.push_back(0);
    }
    
    Interval ONE(1);
    Interval ZERO;
    
    degrees[continuousProblem.getIDForStateVar("x2")+1] = 1;
    Monomial monomial(ONE, degrees);

    Polynomial *p = new Polynomial(monomial);
    *p *= (*coef);
    return p;
  }


  //construct a derivative TaylorModelVector
  TaylorModelVec* derivativeTMV(ContinuousReachability continuousProblem) {
    int numVars = continuousProblem.stateVarNames.size();
    logger.log("derTMV");
    TaylorModelVec *myTMV = new TaylorModelVec;
    TaylorModel tmTemp;

    for(int i=0; i<numVars; ++i)
    {
      myTMV->tms.push_back(tmTemp);
    }
    
    Interval intZeroX1;
    int idX1 = continuousProblem.getIDForStateVar("x1");
    Polynomial *p1 = compApp::dX1(continuousProblem);
    TaylorModel tmX1(*p1, intZeroX1);
    myTMV->tms[idX1] = tmX1;
    
    Interval intZeroX2;
    int idX2 = continuousProblem.getIDForStateVar("x2");
    Polynomial *p2 = compApp::dX2(continuousProblem);
    TaylorModel tmX2(*p2, intZeroX2);
    myTMV->tms[idX2] = tmX2;
    
    Interval intZeroX3;
    int idX3 = continuousProblem.getIDForStateVar("x3");
    Polynomial *p3 = compApp::dX3(continuousProblem);
    TaylorModel tmX3(*p3, intZeroX3);
    myTMV->tms[idX3] = tmX3;
    //is it safe to delete tmXi?
    
    return myTMV;
  }
}

int compMain() {
	logger.log("compapp <");
	logger.inc();
	
	SimpleCompReachability continuousProblem;
	
	//declare all variables
	continuousProblem.declareStateVar("x1");
	continuousProblem.declareStateVar("x2");
	continuousProblem.declareStateVar("x3");
	
	//setttings
	Interval estimationI(-0.001, 0.001);
	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i) {
    logger.log(i);
		continuousProblem.estimation.push_back(estimationI);
	}
	continuousProblem.precondition = QR_PRE;
	continuousProblem.bAdaptiveSteps = false;
	continuousProblem.step = 0.2;
	continuousProblem.time = 0.4;
	continuousProblem.bAdaptiveOrders = false;
	continuousProblem.orderType = UNIFORM;
	continuousProblem.orders.push_back(4);
	continuousProblem.globalMaxOrder = 4;
	
	Interval cutoff_threshold(-1e-15,1e-15);
	continuousProblem.cutoff_threshold = cutoff_threshold;

	intervalNumPrecision = 53;
	strcpy(continuousProblem.outputFileName, "and_gate2");

  //plotting
	continuousProblem.outputAxes.push_back(continuousProblem.getIDForStateVar("x1"));
	continuousProblem.outputAxes.push_back(continuousProblem.getIDForStateVar("x2"));
	continuousProblem.plotSetting = PLOT_OCTAGON;
	continuousProblem.plotFormat = PLOT_GNUPLOT;
	//end of settings



	
	
	//odes
	TaylorModelVec *myTMV = compApp::derivativeTMV(continuousProblem);
	
	
	continuousProblem.declareTMVar("local_t");
	continuousProblem.declareTMVar("local_var_1");
	continuousProblem.declareTMVar("local_var_2");
	continuousProblem.declareTMVar("local_var_3");
	
	vector<Interval> *initVec = new vector<Interval>(continuousProblem.stateVarNames.size());
	(*initVec)[continuousProblem.getIDForStateVar("x1")] = Interval(0.9,1.1);
	(*initVec)[continuousProblem.getIDForStateVar("x2")] = Interval(0.0,0.0);
	(*initVec)[continuousProblem.getIDForStateVar("x3")] = Interval(0.0,0.0);
	
  
  logger.logTMV("tmv", *myTMV);
  logger.logVI("init", *initVec);
  
	Interval intZero;
	Flowpipe *initCond = new ExtractedPicard(*initVec, intZero);
  
	
	SimpleCompSystem system(*myTMV, *initCond);
  
  //logger.log(sbuilder() << "type1: " << typeid(system).name());
  continuousProblem.system = system;
  continuousProblem.pSystem = &system;
  //logger.log(sbuilder() << "type2: " << typeid(continuousProblem.system).name());
  //logger.log(sbuilder() << "type3: " << typeid(*(continuousProblem.pSystem)).name());
	continuousProblem.integrationScheme = ONLY_PICARD;
  
	//continuousProblem.run();
  continuousProblem.myRun();
	
	logger.dec();
	logger.log("compapp >");
	return 0;
}
