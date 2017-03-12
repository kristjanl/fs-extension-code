#include <iostream>

#include "include.h"
#include "Continuous.h"

#include "MyLogger.h"
#include "ExtractedPicard.h"
#include "SimpleImpl.h"
//#include "SimpleComp.h"

using namespace std;


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

//x2' = 0.2*x1*x2
Polynomial* dX2(ContinuousReachability continuousProblem) {
	int numVars1 = continuousProblem.stateVarNames.size()+1;
	Interval coefI(0.2);
	Polynomial *coef = new Polynomial(coefI, numVars1);
	
	vector<int> degrees;
	for(int i=0; i<numVars1; ++i) {
		degrees.push_back(0);
	}
	
	Interval ONE(1);
	Interval ZERO;
	
	degrees[continuousProblem.getIDForStateVar("x1")+1] = 1;
	degrees[continuousProblem.getIDForStateVar("x2")+1] = 1;
	Monomial monomial(ONE, degrees);

	Polynomial *p = new Polynomial(monomial);
	*p *= (*coef);
	return p;
}

//x3' = -0.3*x3
Polynomial* dX3(ContinuousReachability continuousProblem) {
	int numVars1 = continuousProblem.stateVarNames.size()+1;
	Interval coefI(-0.3);
	Polynomial *coef = new Polynomial(coefI, numVars1);
	
	vector<int> degrees;
	for(int i=0; i<numVars1; ++i) {
		degrees.push_back(0);
	}
	
	Interval ONE(1);
	Interval ZERO;
	
	degrees[continuousProblem.getIDForStateVar("x3")+1] = 1;
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
	Polynomial *p1 = dX1(continuousProblem);
	TaylorModel tmX1(*p1, intZeroX1);
	myTMV->tms[idX1] = tmX1;
	
	Interval intZeroX2;
	int idX2 = continuousProblem.getIDForStateVar("x2");
	Polynomial *p2 = dX2(continuousProblem);
	TaylorModel tmX2(*p2, intZeroX2);
	myTMV->tms[idX2] = tmX2;
	
	Interval intZeroX3;
	int idX3 = continuousProblem.getIDForStateVar("x3");
	Polynomial *p3 = dX3(continuousProblem);
	TaylorModel tmX3(*p3, intZeroX3);
	myTMV->tms[idX3] = tmX3;
	//is it safe to delete tmXi?
	
	logger.log(sbuilder() << "&tmX1: " << &(myTMV->tms[idX1]));
	logger.log(sbuilder() << "&tmX2: " << &(myTMV->tms[idX2]));
	logger.log(sbuilder() << "&tmX3: " << &(myTMV->tms[idX3]));
	return myTMV;
}

int simpleImplMain() {

	logger.log("simpleapp <");
	logger.inc();
	
	SimpleImplReachability continuousProblem;
	
	//declare all variables
	continuousProblem.declareStateVar("x1");
	continuousProblem.declareStateVar("x2");
	continuousProblem.declareStateVar("x3");
	
	//setttings
	Interval estimationI(-0.001, 0.001);
	for(int i=0; i<continuousProblem.stateVarNames.size(); ++i) {
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
	TaylorModelVec *myTMV = derivativeTMV(continuousProblem);
	logger.log(sbuilder() << "&tmX1: " << &(myTMV->tms[0]));
	logger.log(sbuilder() << "&tmX2: " << &(myTMV->tms[1]));
	logger.log(sbuilder() << "&tmX3: " << &(myTMV->tms[2]));
	
	
	continuousProblem.declareTMVar("local_t");
	continuousProblem.declareTMVar("local_var_1");
	continuousProblem.declareTMVar("local_var_2");
	continuousProblem.declareTMVar("local_var_3");
	
	vector<Interval> *initVec = new vector<Interval>(continuousProblem.stateVarNames.size());
	(*initVec)[continuousProblem.getIDForStateVar("x1")] = Interval(0.9,1.1);
	(*initVec)[continuousProblem.getIDForStateVar("x2")] = Interval(0.25,0.25);
	(*initVec)[continuousProblem.getIDForStateVar("x3")] = Interval(2.0,2.2);
	
	Interval intZero;
	Flowpipe *initCond = new ExtractedPicard(*initVec, intZero);
  
	
	SimpleImplSystem system(*myTMV, *initCond);
  
	continuousProblem.pSystem = &system;
	continuousProblem.integrationScheme = ONLY_PICARD;	
	
	//continuousProblem.run();
  continuousProblem.myRun();
	
	logger.dec();
	logger.log("simpleapp >");
	return 0;
}
