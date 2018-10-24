/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef CONTINUOUS_H_
#define CONTINUOUS_H_

#include "TaylorModel.h"
#include "Geometry.h"
#include "MyLogger.h"
#include "OutputWriter.h"
#include "Utils.h"


// extern vector<LinearConstraint> constraint_template;

class Flowpipe					// A flowpipe is represented by a composition of two Taylor models. The left Taylor model is the preconditioning part.
{
public:
	TaylorModelVec tmvPre;		// preconditioning Taylor model
	TaylorModelVec tmv;
public:
	vector<Interval> domain;	// domain of TMV_right, the first variable is t
public:
	Flowpipe();
	Flowpipe(const TaylorModelVec & tmvPre_input, const TaylorModelVec & tmv_input, const vector<Interval> & domain_input);
	Flowpipe(const vector<Interval> & box, const Interval & I);								// represent a box
	Flowpipe(const TaylorModelVec & tmv_input, const vector<Interval> & domain_input);		// construct a flowpipe from a Taylor model
	Flowpipe(const Flowpipe & flowpipe);
	~Flowpipe();

	void clear();
	void dump(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, const Interval & cutoff_threshold) const;
	void dump_normal(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames, vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;
	void composition(TaylorModelVec & result, const Interval & cutoff_threshold) const;	// apply the preconditioning part to the Taylor model
	void composition_normal(TaylorModelVec & result, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;

	void intEval(vector<Interval> & result, const Interval & cutoff_threshold) const;
	void intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;

	void normalize();


	// Taylor model integration by only using Picard operation
	// fixed step sizes and orders
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_picard(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;


	// fast integration scheme for low-degree ODEs
	// fixed step sizes and orders
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_low_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const vector<HornerForm> & taylorExpansion, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// integration scheme for high-degree ODEs
	// fixed step sizes and orders
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_high_degree(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;


	// integration scheme for non-polynomial ODEs (using Taylor approximations)
	// fixed step sizes and orders
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const double miniStep, const vector<int> & orders, const int globalMaxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, int & order, const int maxOrder, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;
	int advance_non_polynomial_taylor(Flowpipe & result, const vector<string> & strOde, const vector<string> & strOde_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, vector<int> & orders, const int localMaxOrder, const vector<int> & maxOrders, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;


  virtual int advance_picard2(Flowpipe & result, const vector<HornerForm> & ode, const vector<HornerForm> & ode_centered, const int precondition, vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const int order, const vector<Interval> & estimation, const vector<PolynomialConstraint> & invariant, const Interval & cutoff_threshold) const;

	Flowpipe & operator = (const Flowpipe & flowpipe);

	friend class ContinuousSystem;
	friend class ContinuousReachability;
	friend class HybridSystem;
	friend class HybridReachability;
  friend class ExtractedPicard;
  friend class SimpleCompSystem;
	friend class SimpleCompReachability;
  friend class SmallCompSystem;
	friend class SmallCompReachability;
  friend class SimpleImplSystem;
	friend class SimpleImplReachability;
};

class MySettings;

class ContinuousSystem
{
public:
	TaylorModelVec tmvOde;
	TaylorModelVec tmvOde_centered;
	vector<HornerForm> hfOde;		// a Horner form of the ode
	vector<HornerForm> hfOde_centered;
	Flowpipe initialSet;			// the initial set
	vector<string> strOde;
	vector<string> strOde_centered;

public:
	ContinuousSystem();
	ContinuousSystem(const TaylorModelVec & ode_input, const Flowpipe & initialSet_input);
	ContinuousSystem(const vector<string> & strOde_input, const Flowpipe & initialSet_input);
	ContinuousSystem(const ContinuousSystem & system);
	~ContinuousSystem();
  
  
	MySettings *settings;

	// efficient integration method for linear ODEs
	// since we can always compute a safe remainder, adaptive techniques are not needed
	// since the size of a Taylor series is linear w.r.t. the number of state variables, we may simply use uniform orders
	void reach_linear(list<Flowpipe> & results, const double step, const double time, const int order, const bool bPrint, const Interval & cutoff_threshold) const;

	// only use Picard operation
	// fixed step sizes and orders
	void reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_picard(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	void reach_picard(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_picard(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	void reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_picard(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;


	// for low-degree ODEs
	// fixed step sizes and orders
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	void reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_low_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// for high-degree ODEs
	// fixed step sizes and orders
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	void reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_high_degree(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// for non-polynomial ODEs (using Taylor approximations)
	// fixed step sizes and orders
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;
	void reach_non_polynomial_taylor(list<Flowpipe> & results, const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	ContinuousSystem & operator = (const ContinuousSystem & system);

  virtual void my_reach_picard(list<Flowpipe> & results, const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames, const Interval & cutoff_threshold, OutputWriter & writer) const;
  
	friend class SimpleImplSystem;
	friend class SimpleCompSystem;
  friend class SmallCompSystem;
	friend class ContinuousReachability;
	friend class SimpleCompReachability;
	friend class SmallCompReachability;
	friend class SimpleImplReachability;
};

class SimpleCompReachability;
class SimpleCompSystem;
class SmallCompReachability;
class SmallCompSystem;
class SimpleImplReachability;
class SimpleImplSystem;
class ShrinkWrappingCondition;
class Transformer;

class ContinuousReachability		// The reachability analysis of continuous systems
{
public:
	ContinuousSystem system;		// the continuous system
	ContinuousSystem *pSystem;		// the continuous system
	double step;					// the step size used in the reachability analysis
	double time;					// the time horizon for the reachability analysis
	int precondition;				// the preconditioning technique
	MySettings *settings;
	vector<int> outputAxes;			// the output axes
	int plotSetting;
	int plotFormat;
	int numSections;				// the number of sections in each dimension

	int orderType;
	bool bAdaptiveSteps;
	bool bAdaptiveOrders;

	vector<Interval> estimation;	// the remainder estimation for varying time step
	double miniStep;				// the minimum step size
	vector<int> orders;				// the order(s)
	vector<int> maxOrders;			// the maximum orders
	int globalMaxOrder;

	bool bPrint;
	bool bSafetyChecking;

	int integrationScheme;
  int algorithm;

	Interval cutoff_threshold;

	list<Flowpipe> flowpipes;
	list<TaylorModelVec> flowpipesCompo;
	list<vector<Interval> > domains;

	vector<PolynomialConstraint> unsafeSet;

	map<string,int> stateVarTab;
	vector<string> stateVarNames;

	map<string,int> tmVarTab;
	vector<string> tmVarNames;

	map<string,int> parTab;
	vector<string> parNames;
	vector<Interval> parRanges;

	char outputFileName[NAME_SIZE];
public:
	ContinuousReachability();
	~ContinuousReachability();
  

	void dump(FILE *fp) const;

	void run();
	void contRun();
  virtual void myRun();
  
	void composition();
	int safetyChecking() const;
	unsigned long numOfFlowpipes() const;

	void dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, const list<Interval> & globalTimes) const;

	void plot_2D() const;

	void plot_2D_GNUPLOT(FILE *fp) const;
	void plot_2D_interval_GNUPLOT(FILE *fp) const;
	void plot_2D_octagon_GNUPLOT(FILE *fp) const;
	void plot_2D_grid_GNUPLOT(FILE *fp) const;

	void plot_2D_MATLAB(FILE *fp) const;
	void plot_2D_interval_MATLAB(FILE *fp) const;
	void plot_2D_octagon_MATLAB(FILE *fp) const;
	void plot_2D_grid_MATLAB(FILE *fp) const;

	bool declareStateVar(const string & vName);
	int getIDForStateVar(const string & vName) const;
	bool getStateVarName(string & vName, const int id) const;

	bool declareTMVar(const string & vName);
	int getIDForTMVar(const string & vName) const;
	bool getTMVarName(string & vName, const int id) const;

	bool declarePar(const string & pName, const Interval & range);
	int getIDForPar(const string & pName) const;
	bool getParName(string & pName, const int id) const;
	bool getRangeForPar(Interval & range, const string & pName) const;
  SimpleCompReachability createSimpleComp();
  SmallCompReachability createSmallComp();
  SimpleImplReachability createSimpleImpl();
  friend class SimpleCompReachability;
  friend class SmallCompReachability;
  friend class SimpleImplReachability;
};
#include "SimpleComp.h"
#include "SmallComp.h"
#include "SimpleImpl.h"

void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const int order, const Interval & cutoff_threshold);
void computeTaylorExpansion(TaylorModelVec & result, const TaylorModelVec & first_order_deriv, const TaylorModelVec & ode, const vector<int> & orders, const Interval & cutoff_threshold);

void construct_step_exp_table(vector<Interval> & step_exp_table, vector<Interval> & step_end_exp_table, const double step, const int order);
void construct_step_exp_table(vector<Interval> & step_exp_table, const Interval & step, const int order);

void preconditionQR(Matrix & result, const TaylorModelVec & tmv, const int rangeDim, const int domainDim);

Interval rho(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & domain);
Interval rhoNormal(const TaylorModelVec & tmv, const vector<Interval> & l, const vector<Interval> & step_end_exp_table);

Interval rho(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & domain);
Interval rhoNormal(const TaylorModelVec & tmv, const RowVector & l, const vector<Interval> & step_end_exp_table);

void templatePolyhedron(Polyhedron & result, const TaylorModelVec & tmv, const vector<Interval> & domain);
void templatePolyhedronNormal(Polyhedron & result, const TaylorModelVec & tmv, vector<Interval> & step_end_exp_table);

/*
int intersection_check_interval_arithmetic(const list<PolynomialConstraint> & pcs, const list<HornerForm> & objFuncs, const list<Interval> & remainders, const vector<Interval> & domain, list<bool> & bNeeded);
bool boundary_intersected_collection(const vector<PolynomialConstraint> & pcs, const vector<HornerForm> & objFuncs, const vector<Interval> & remainders, const vector<Interval> & domain, vector<bool> & boundary_intersected);
*/

// domain contraction by using interval arithmetic
int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const vector<PolynomialConstraint> & pcs, vector<bool> & boundary_intersected, const Interval & cutoff_threshold);

int contract_interval_arithmetic(TaylorModelVec & flowpipe, vector<Interval> & domain, const Polyhedron & inv, vector<bool> & boundary_intersected);

int contract_remainder(const vector<Interval> & polyRange, vector<Interval> & remainders, const vector<HornerForm> & hfs, const vector<Interval> & b);

void gridBox(list<vector<Interval> > & grids, const vector<Interval> & box, const int num);

void exp_int_mat(Interval_matrix & result_ts, Interval_matrix & result_rem, const Interval_matrix & A, const int order);
void int_exp_int_mat(Interval_matrix & result_ts, Interval_matrix & result_rem, const Interval_matrix & A, const double step, const int order);


class MySettings {
  public:
    OutputWriter *writer;
    vector< vector<int> > intComponents;
    bool autoComponents;
    int order;
    double step; 
    double time;
    vector<Interval> estimation;
    vector<Interval> step_exp_table;
    vector<Interval> step_end_exp_table;
    vector<Interval> domain;
    Interval cutoff;
    bool useFlow;
    bool discardEmptyParams;
    vector<string> varNames;
    Transformer *transformer; // determines how are initials sets transformed for each timestep
    MySettings();
    MySettings(OutputWriter *writer, int order, double step, 
        double time, vector<Interval> estimation, 
        vector<Interval> step_exp_table, 
        vector<Interval> step_end_exp_table, 
        vector<Interval> domain, const Interval cutoff);
    void log();
  
};

//void parseODE();


// void build_constraint_template(const int d);

#endif /* CONTINUOUS_H_ */
