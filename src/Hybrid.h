/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef HYBRID_H_
#define HYBRID_H_

#include "Continuous.h"

extern ParseSetting parseSetting;
extern ParseResult parseResult;

class ResetMap
{
public:
	TaylorModelVec tmvReset;
public:
	ResetMap();
	ResetMap(const TaylorModelVec & tmv);
	ResetMap(const ResetMap & reset);
	~ResetMap();

	void reset(TaylorModelVec & result, const TaylorModelVec & tmv, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const;

	ResetMap & operator = (const ResetMap & reset);
};

class DiscTrans
{
public:
	int jumpID;
	int startID;
	int targetID;
	vector<PolynomialConstraint> guard;
	ResetMap resetMap;
public:
	DiscTrans();
	DiscTrans(const int id, const int start, const int target, const vector<PolynomialConstraint> & lcs, const ResetMap & reset);
	DiscTrans(const DiscTrans & trans);
	~DiscTrans();

	DiscTrans & operator = (const DiscTrans & trans);
};

class TreeNode
{
public:
	int jumpID;
	int modeID;
	Interval localTime;
	TreeNode *parent;
	list<TreeNode *> children;

	TreeNode(const int jump, const int mode, const Interval & t);
	TreeNode(const TreeNode & node);
	~TreeNode();

	void dump(FILE *fp, const string & prefix, const vector<string> & modeNames) const;

	TreeNode & operator = (const TreeNode & node);
};


class ReachabilitySetting
{
public:
	double step;
	int precondition;
	int orderType;
	bool bAdaptiveSteps;
	bool bAdaptiveOrders;
	vector<Interval> estimation;
	double miniStep;
	vector<int> orders;
	vector<int> maxOrders;
	int globalMaxOrder;
	Interval cutoff_threshold;

public:
	ReachabilitySetting();
	~ReachabilitySetting();
	ReachabilitySetting(const ReachabilitySetting & setting);
	void clear();

	ReachabilitySetting & operator = (const ReachabilitySetting & setting);
};



class HybridSystem
{
private:
	vector<int> modes;
	vector<TaylorModelVec> odes;
	vector<vector<HornerForm> > hfOdes;
	vector<vector<string> > strOdes;
	vector<TaylorModelVec> odes_centered;
	vector<vector<HornerForm> > hfOdes_centered;
	vector<vector<string> > strOdes_centered;
	vector<vector<PolynomialConstraint> > invariants;
	vector<Polyhedron> linear_invariants;	// only used by linear systems
	vector<vector<DiscTrans> > transitions;
	int initialMode;
	Flowpipe initialSet;
public:
	HybridSystem();
	HybridSystem(const HybridSystem & hybsys);
	~HybridSystem();


	// efficient flowpipe construction method for linear continuous dynamics
	void reach_continuous_linear(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const bool bPrint, vector<bool> & invariant_boundary_intersected,
			const vector<string> & modeNames, const vector<string> & stateVarNames, const Interval & cutoff_threshold) const;

	// only use Picard operation
	// fixed step sizes and orders
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint,
			const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const int order, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
			const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected,
			const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_picard(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;


	// for low-degree ODEs
	// fixed step sizes and orders
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint,
			const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const int order, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
			const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected,
			const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_low_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// for high-degree ODEs
	// fixed step sizes and orders
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint,
			const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const int order, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
			const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected,
			const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_high_degree(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// for non-polynomial ODEs (using Taylor approximations)
	// fixed step sizes and orders
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int precondition, const vector<Interval> & estimation, const bool bPrint,
			const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const int globalMaxOrder, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive step sizes and fixed orders
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const int order, const int precondition,
			const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double miniStep, const double time, const vector<int> & orders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// adaptive orders and fixed step sizes
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const int order, const int maxOrder, const int precondition, const vector<Interval> & estimation,
			const bool bPrint, const vector<string> & stateVarNames, vector<bool> & invariant_boundary_intersected,
			const vector<string> & modeNames, const Interval & cutoff_threshold) const;
	bool reach_continuous_non_polynomial_taylor(list<TaylorModelVec> & resultsCompo, list<vector<Interval> > & domains, const int mode, const Flowpipe & initFp,
			const double step, const double time, const vector<int> & orders, const vector<int> & maxOrders, const int globalMaxOrder,
			const int precondition, const vector<Interval> & estimation, const bool bPrint, const vector<string> & stateVarNames,
			vector<bool> & invariant_boundary_intersected, const vector<string> & modeNames, const Interval & cutoff_threshold) const;

	// hybrid reachability
	void reach_hybrid(list<list<TaylorModelVec> > & resultsCompo, list<list<vector<Interval> > > & domains, list<int> & modeIDs,
			list<TreeNode *> & traceNodes, TreeNode * & traceTree, const vector<int> & integrationSchemes, const double time, const int maxJmps,
			const ReachabilitySetting & global_setting, const vector<ReachabilitySetting> & local_settings, const vector<vector<int> > & aggregType,
			const vector<vector<vector<RowVector> > > aggregationTemplate_candidates, const vector<RowVector> default_aggregation_template,
			const vector<vector<Matrix> > & weightTab, const vector<vector<vector<bool> > > & linear_auto, const vector<vector<vector<RowVector> > > & template_auto,
			const bool bPrint, const vector<string> & stateVarNames, const vector<string> & modeNames, const vector<string> & tmVarNames) const;


	HybridSystem & operator = (const HybridSystem & hybsys);

	friend class HybridReachability;
};




class HybridReachability
{
public:
	HybridSystem system;		// the hybrid system
//	double step;				// the step size used in the reachability analysis
	double time;				// the time horizon for the reachability analysis
//	int precondition;			// the preconditioning technique
	vector<int> outputAxes;		// the output axes
	int plotSetting;
	int plotFormat;
	int numSections;			// the number of sections in each dimension

//	int orderType;
//	bool bAdaptiveSteps;
//	bool bAdaptiveOrders;

//	vector<Interval> estimation;// the remainder estimation for varying time step
//	double miniStep;			// the minimum step size
//	vector<int> orders;			// the order(s)
//	vector<int> maxOrders;		// the maximum orders
//	int globalMaxOrder;

	int maxJumps;
	int numOfJumps;

	bool bPrint;
	bool bSafetyChecking;

	vector<int> integrationSchemes;


//	vector<double> mode_steps;
//	vector<int> model_preconditions;
//	vector<int> mode_orderTypes;
//	vector<bool> mode_bAdaptiveSteps;
//	vector<bool> mode_bAdaptiveOrders;
//	vector<vector<Interval> > mode_estimations;
//	vector<double> mode_miniSteps;
//	vector<vector<int> > mode_orders;
//	vector<vector<int> > mode_maxOrders;
//	vector<int> mode_globalMaxOrder;
//	vector<Interval> mode_cutoff_threshold;


	ReachabilitySetting global_setting;
	vector<ReachabilitySetting> local_settings;

	TreeNode *traceTree;

	vector<bool> bVecUnderCheck;

	vector<vector<int> > aggregationType;
	vector<RowVector> default_aggregation_template;
	vector<vector<vector<RowVector> > > aggregationTemplate_candidates;

	// information for template selection
	vector<vector<vector<bool> > > linear_auto;
	vector<vector<vector<RowVector> > > template_auto;

	vector<vector<Matrix> > weightTab;

//	vector<vector<TaylorModelVec> > aggregationTemplate_TaylorModel;	// will be considered later

	list<list<TaylorModelVec> > flowpipesCompo;
	list<list<vector<Interval> > > domains;
	list<int> modeIDs;
	list<TreeNode *> traceNodes;

	map<string,int> stateVarTab;
	vector<string> stateVarNames;

	map<string,int> tmVarTab;
	vector<string> tmVarNames;

	map<string,int> modeTab;
	vector<string> modeNames;

	map<string,int> parTab;
	vector<string> parNames;
	vector<Interval> parRanges;

	vector<vector<PolynomialConstraint> > unsafeSet;

	char outputFileName[NAME_SIZE];
public:
	HybridReachability();
	~HybridReachability();

	void dump(FILE *fp) const;

	void run();

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

	bool declareMode(const string & mName, const TaylorModelVec & ode, const vector<PolynomialConstraint> & inv, const int integrationScheme, const ReachabilitySetting & local_setting);
	bool declareMode(const string & mName, const vector<string> & strOde, const vector<PolynomialConstraint> & inv, const int integrationScheme, const ReachabilitySetting & local_setting);

	int getIDForMode(const string & mName) const;
	bool getModeName(string & mName, const int id) const;

	void declareTrans(const int start, const int end, const vector<PolynomialConstraint> & guard, const ResetMap & reset, const int aggregType, const vector<vector<double> > & candidates);
	void declareTrans();

	void initialConfig(const int modeID, const Flowpipe & initialSet);
	void set_default_template();
	void constructWeightTab();

	int safetyChecking();
	unsigned long numOfFlowpipes() const;
	void dump_potential_counterexample(FILE *fp, const list<TaylorModelVec> & flowpipes, const list<vector<Interval> > & domains, TreeNode * const node, const list<Interval> & globalTimes) const;

	// parallelotopic aggregation
	friend void aggregate_flowpipes_by_Parallelotope(TaylorModelVec & tmvAggregation, vector<Interval> & doAggregation, const vector<TaylorModelVec> & flowpipes,
			const vector<vector<Interval> > & domains, const vector<PolynomialConstraint> & invariant, const DiscTrans & jump, vector<bool> & boundary_intersected,
			const vector<RowVector> & template_candidates, const vector<RowVector> & template_default, const vector<vector<Matrix> > & weightTab,
			const vector<vector<vector<bool> > > & linear_auto, const vector<vector<vector<RowVector> > > & template_auto, const int globalMaxOrder, const int rangeDim);
};

class FactorTab
{
public:
	int index;
	Interval factor;
	Interval intercept;
public:
	FactorTab();
	FactorTab(const int index_input, const Interval & factor_input, const Interval & intercept_input);
	~FactorTab();

	friend bool compareFactor(const FactorTab & a, const FactorTab & b);
	friend bool compareIntercept(const FactorTab & a, const FactorTab & b);
};

void generateNodeSeq(list<TreeNode *> & result, TreeNode *root);

// interval aggregation
void aggregate_flowpipes_by_interval(TaylorModelVec & tmvAggregation, vector<Interval> & doAggregation, const vector<TaylorModelVec> & flowpipes, const vector<vector<Interval> > & domains);

bool vector_selection(FactorTab & lst_selected, list<FactorTab> & lst_unselected, Matrix & matTemplate, const vector<RowVector> & rowVecs, int & rank);
bool check_validity(Matrix & matTemplate, const RowVector & rowVec, const int rank);

void aggregate_flowpipes_by_Parallelotope(TaylorModelVec & tmvAggregation, vector<Interval> & doAggregation, const vector<TaylorModelVec> & flowpipes,
		const vector<vector<Interval> > & domains, const vector<PolynomialConstraint> & invariant, const DiscTrans & jump, vector<bool> & boundary_intersected,
		const vector<RowVector> & template_candidates, const vector<RowVector> & template_default, const Matrix & weightTab,
		const vector<bool> & linear_auto, const vector<RowVector> & template_auto, const int globalMaxOrder, const int rangeDim);


#endif /* HYBRID_H_ */
