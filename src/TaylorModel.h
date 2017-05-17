/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef TAYLORMODEL_H_
#define TAYLORMODEL_H_

#include "Polynomial.h"

class TaylorModel			// Taylor models: R^n -> R. We use t to denote the time variable and x to denote the state variable.
{
public:
	Polynomial expansion;	// Taylor expansion
	Interval remainder;		// remainder interval
public:
	TaylorModel();			// empty Taylor model.
	TaylorModel(const Interval & I, const int numVars);				// constant
	TaylorModel(const Polynomial & polyExp, const Interval & I);	// Taylor model (P,I)
	TaylorModel(const RowVector & coefficients);
	TaylorModel(const RowVector & coefficients, const Interval & I);
	TaylorModel(const vector<Interval> & coefficients);
	TaylorModel(const vector<Interval> & coefficients, const Interval & I);
	TaylorModel(const TaylorModel & tm);
	~TaylorModel();

	void clear();
	void dump_interval(FILE *fp, const vector<string> & varNames) const;
	void dump_constant(FILE *fp, const vector<string> & varNames) const;
	void constant(Interval & result) const;									// Return the constant part of the expansion.

	void intEval(Interval & result, const vector<Interval> & domain) const;
	void intEvalNormal(Interval & result, const vector<Interval> & step_exp_table) const;

	void ctrunc(const vector<Interval> & domain, const int order);		// conservative truncation
	void nctrunc(const int order);										// non-conservative truncation
	void ctrunc_normal(const vector<Interval> & step_exp_table, const int order);

	void inv(TaylorModel & result) const;									// additive inverse
	void inv_assign();

	void add(TaylorModel & result, const TaylorModel & tm) const;			// addition
	void sub(TaylorModel & result, const TaylorModel & tm) const;			// subtraction
	void add_assign(const TaylorModel & tm);
	void sub_assign(const TaylorModel & tm);

	void mul_ctrunc(TaylorModel & result, const TaylorModel & tm, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const;
	void mul_ctrunc_normal(TaylorModel & result, const TaylorModel & tm, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const;
	void mul_no_remainder(TaylorModel & result, const TaylorModel & tm, const int order, const Interval & cutoff_threshold) const;
	void mul_no_remainder_no_cutoff(TaylorModel & result, const TaylorModel & tm, const int order) const;
	void mul(TaylorModel & result, const Interval & I) const;

	void mul_ctrunc_assign(const TaylorModel & tm, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold);
	void mul_ctrunc_normal_assign(const TaylorModel & tm, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);
	void mul_no_remainder_assign(const TaylorModel & tm, const int order, const Interval & cutoff_threshold);
	void mul_no_remainder_no_cutoff_assign(const TaylorModel & tm, const int order);
	void mul_assign(const Interval & I);

	void mul_insert(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void mul_insert_normal(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;
	void mul_insert_ctrunc(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const;
	void mul_insert_ctrunc_normal(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const;
	void mul_insert_ctrunc_normal(TaylorModel & result, Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const;

	void mul_insert_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold);
	void mul_insert_normal_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const Interval & cutoff_threshold);
	void mul_insert_ctrunc_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold);
	void mul_insert_ctrunc_normal_assign(const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);
	void mul_insert_ctrunc_normal_assign(Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);

	void div(TaylorModel & result, const Interval & I) const;
	void div_assign(const Interval & I);

	void derivative(TaylorModel & result, const int varIndex) const;		// derivative with respect to a variable

	// Lie derivative, the vector field is given by f
	void LieDerivative_no_remainder(TaylorModel & result, const TaylorModelVec & f, const int order, const Interval & cutoff_threshold) const;

	void integral(TaylorModel & result, const Interval & I) const;				// Integral with respect to t
	void integral_no_remainder(TaylorModel & result) const;

	void linearCoefficients(vector<Interval> & result) const;
	void toHornerForm(HornerForm & result, Interval & I) const;

	void insert(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void insert_ctrunc(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const;
	void insert_no_remainder(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const;
	void insert_ctrunc_normal(TaylorModel & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;

	void evaluate_t(TaylorModel & result, const vector<Interval> & step_exp_table) const;			// evaluate the Taylor model at time t

	void mul(TaylorModel & result, const int varIndex, const int degree) const;		// multiplied by a term x^d
	void mul_assign(const int varIndex, const int degree);

	void rmConstant();
	void cutoff_normal(const vector<Interval> & step_exp_table, const Interval & cutoff_threshold);
	void cutoff(const vector<Interval> & domain, const Interval & cutoff_threshold);
	void cutoff(const Interval & cutoff_threshold);

	int degree() const;
	bool isZero() const;

	void center_nc();

	void rmZeroTerms(const vector<int> & indices);

	void normalize(vector<Interval> & domain);

	void polyRange(Interval & result, const vector<Interval> & domain) const;
	void polyRangeNormal(Interval & result, const vector<Interval> & step_exp_table) const;

	void exp_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void rec_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void sin_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void cos_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void log_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void sqrt_taylor(TaylorModel & result, list<Interval> & ranges, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;

	Interval getRemainder() const;
	void getExpansion(Polynomial & P) const;

  string toString(const vector<string> & varNames) const;
  TaylorModel transform(map<int, int> lookup, int size);
  TaylorModel transform(vector<int> indexes);
  TaylorModel* ptransform(vector<int> indexes);
  TaylorModel prepareSecondary(int prefix) const;
  int getParamCount() const;
  TaylorModel addNParams(int n) const;
  vector<int> getParams() const;
  void getLinearPart(Interval ret[], vector<int> variables) const;
  void getLinearPart2(vector<int> variables) const;
  
  //partitions the TaylorModel into constant, linear, nonlinear and remainder parts
  void getParts(TaylorModel & constant, TaylorModel & linear, 
      TaylorModel & nonLinear, TaylorModel & remainder) const;
  
  //
  void removeHighTerms(int order);
  
  bool isClose(const TaylorModel & tm, double d) const;

	TaylorModel & operator = (const TaylorModel & tm);

	friend class HornerForm;
	friend class Polynomial;
	friend class TaylorModelVec;
	friend class Flowpipe;
	friend class ContinuousSystem;
	friend class ContinuousReachability;
	friend class HybridSystem;
	friend class HybridReachability;
};

class TaylorModelVec			// Taylor models: R^n -> R^m
{
public:
	vector<TaylorModel> tms;
public:
	TaylorModelVec();
	TaylorModelVec(const vector<TaylorModel> & tms_input);
	TaylorModelVec(const vector<Interval> & constants, const int numVars);
	TaylorModelVec(const Matrix & coefficients);
	TaylorModelVec(const Matrix & coefficients, const vector<Interval> & remainders);
	TaylorModelVec(const vector<vector<Interval> > & coefficients);
	TaylorModelVec(const vector<vector<Interval> > & coefficients, const vector<Interval> & remainders);
	TaylorModelVec(const vector<Interval> & intVec, vector<Interval> & domain);
	TaylorModelVec(const TaylorModelVec & tmv);
	~TaylorModelVec();

	void clear();
	void dump_interval(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const;
	void dump_constant(FILE *fp, const vector<string> & stateVarNames, const vector<string> & tmVarNames) const;
	void constant(vector<Interval> & result) const;

	void intEval(vector<Interval> & result, const vector<Interval> & domain) const;
	void intEvalNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const;

	void ctrunc(const vector<Interval> & domain, const int order);		// conservative truncation
	void nctrunc(const int order);										// non-conservative truncation
	void ctrunc_normal(const vector<Interval> & step_exp_table, const int order);

	void ctrunc(const vector<Interval> & domain, const vector<int> & orders);		// conservative truncation
	void nctrunc(const vector<int> & orders);										// non-conservative truncation
	void ctrunc_normal(const vector<Interval> & step_exp_table, const vector<int> & orders);

	void inv(TaylorModelVec & result) const;									// additive inverse
	void add(TaylorModelVec & result, const TaylorModelVec & tmv) const;		// addition
	void sub(TaylorModelVec & result, const TaylorModelVec & tmv) const;		// subtraction

	void add_assign(const TaylorModelVec & tmv);
	void sub_assign(const TaylorModelVec & tmv);

	void mul(TaylorModelVec & result, const Interval & I) const;				// Multiplication of a Taylor model and a constant.
	void mul_assign(const Interval & I);

	void div(TaylorModelVec & result, const Interval & I) const;				// Divide the Taylor model by a constant.
	void div_assign(const Interval & I);

	void derivative(TaylorModelVec & result, const int varIndex) const;

	void LieDerivative_no_remainder(TaylorModelVec & result, const TaylorModelVec & f, const int order, const Interval & cutoff_threshold) const;
	void LieDerivative_no_remainder(TaylorModelVec & result, const TaylorModelVec & f, const vector<int> & orders, const Interval & cutoff_threshold) const;

	void integral(TaylorModelVec & result, const Interval & I) const;
	void integral_no_remainder(TaylorModelVec & result) const;

	void linearCoefficients(vector<vector<Interval> > & result) const;

	void rmZeroTerms(const vector<int> & indices);

	void insert(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const Interval & cutoff_threshold) const;

	void insert_ctrunc(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const int order, const Interval & cutoff_threshold) const;
	void insert_no_remainder(TaylorModelVec & result, const TaylorModelVec & vars, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void insert_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;

	void insert_ctrunc(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & domain, const vector<int> & orders, const Interval & cutoff_threshold) const;
	void insert_no_remainder(TaylorModelVec & result, const TaylorModelVec & vars, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const;
	void insert_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & vars, const vector<Interval> & varsPolyRange, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const;

	void evaluate_t(TaylorModelVec & result, const vector<Interval> & step_exp_table) const;

	void mul(TaylorModelVec & result, const int varIndex, const int degree) const;
	void mul_assign(const int varIndex, const int degree);

  //multiplies TMV from the left using A (A*tmv)
	void linearTrans(TaylorModelVec & result, const Matrix & A) const;		// linear transformation
	void linearTrans_assign(const Matrix & A);

	void rmConstant();
	void cutoff_normal(const vector<Interval> & step_exp_table, const Interval & cutoff_threshold);
	void cutoff(const vector<Interval> & domain, const Interval & cutoff_threshold);
	void cutoff(const Interval & cutoff_threshold);

	void center_nc();

	void Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void Picard_no_remainder_assign(const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const int order, const Interval & cutoff_threshold);
	void Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold);

	void Picard_ctrunc_normal(TaylorModelVec & result, vector<RangeTree *> & trees, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void Picard_ctrunc_normal(TaylorModelVec & result, vector<RangeTree *> & trees, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const;
	void Picard_only_remainder(vector<Interval> & result, vector<RangeTree *> & trees, const TaylorModelVec & x0, const vector<HornerForm> & ode, const Interval & timeStep) const;

	void Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold) const;
	void Picard_no_remainder_assign(const TaylorModelVec & x0, const vector<HornerForm> & ode, const int numVars, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold);
	void Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold) const;
	void Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<Interval> & polyRange, const vector<HornerForm> & ode, const vector<Interval> & step_exp_table, const int numVars, const vector<int> & orders, const Interval & cutoff_threshold);

	// using Taylor approximation
	void Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const int order, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const vector<string> & strOde, const int order, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<int> & orders, const vector<bool> & bIncreased, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const vector<int> & orders, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const vector<string> & strOde, const vector<Interval> & step_exp_table, const vector<int> & orders, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_only_remainder(vector<Interval> & result, const TaylorModelVec & x0, const vector<string> & strOde, const Interval & timeStep, const int order) const;
	void Picard_non_polynomial_taylor_only_remainder(vector<Interval> & result, const TaylorModelVec & x0, const vector<string> & strOde, const Interval & timeStep, const vector<int> & orders) const;

	void normalize(vector<Interval> & domain);		// we assume that the original domain is full-dimensional

	void polyRange(vector<Interval> & result, const vector<Interval> & domain) const;
	void polyRangeNormal(vector<Interval> & result, const vector<Interval> & step_exp_table) const;

  string toString(const vector<string> & varNames) const;
  TaylorModelVec prepareSecondary(int prefix) const;
  TaylorModelVec transform(vector<int> indexes);
  bool compare(const TaylorModelVec & tmv, const vector<Interval> & domain) const;
  
  bool isClose(const TaylorModelVec & tmv, double d) const;
  TaylorModelVec addNParams(int n) const;
  
  //partitions the TaylorModelVec into constant, linear, nonlinear and remainder parts
  void getParts(TaylorModelVec & constant, TaylorModelVec & linear, 
      TaylorModelVec & nonLinear, TaylorModelVec & remainder) const;
  
  //remove higher terms
  void removeHighTerms(int order);
  
  vector<Interval> getRemainders();

	TaylorModelVec & operator = (const TaylorModelVec & tmv);
};

class ParseSetting
{
public:
	string strODE;

	list<Interval> ranges;
	list<Interval>::iterator iterRange;
	vector<string> variables;

	vector<Interval> step_exp_table;
	TaylorModelVec flowpipe;
	int order;

	Interval cutoff_threshold;

	ParseSetting();
	ParseSetting(const ParseSetting & setting);
	~ParseSetting();

	void clear();
	
	void addVar(string v);

	ParseSetting & operator = (const ParseSetting & setting);
};

class ParseResult					// the data structure of a parsed non-polynomial ODE with its remainder
{
public:
	Polynomial	expansion;
	Interval	remainder;
	TaylorModel model;
	TaylorModelVec tmv;
	vector<Polynomial> *polys;
	Monomial mono;
	vector<int> integerVec;
	string		strExpansion;

	ParseResult();
	ParseResult(const ParseResult & result);
	~ParseResult();

	ParseResult & operator = (const ParseResult & result);
};


class Interval_matrix
{
private:
	vector<vector<Interval> > elements;
public:
	Interval_matrix();
	Interval_matrix(const vector<vector<Interval> > & elements_input);
	Interval_matrix(const vector<Interval> & intVector);
	Interval_matrix(const int m, const int n);
	Interval_matrix(const Interval_matrix & A);
	~Interval_matrix();

	int rows() const;
	int cols() const;

	void get(Interval & I, const int i, const int j) const;
	Interval get(const int i, const int j) const;

	void set(const Interval & I, const int i, const int j);

	void mul(Interval_matrix & result, const double v) const;
	void mul_assign(const double v);

	void mul(Interval_matrix & result, const Interval & I) const;
	void mul_assign(const Interval & I);

	void div(Interval_matrix & result, const double v) const;
	void div_assign(const double v);

	void decompose(Interval_matrix & det, Interval_matrix & nondet) const;

	void pow(Interval_matrix & result, const int degree) const;
	void pow_assign(const int degree);

	double max_norm() const;

	void linear_trans(vector<Polynomial> & result, const vector<Polynomial> & polynomial) const;
	void transpose(Interval_matrix & result) const;

	void output(FILE *fp) const;

	Interval_matrix & operator += (const Interval_matrix & A);
	Interval_matrix & operator -= (const Interval_matrix & A);
	Interval_matrix & operator *= (const Interval_matrix & A);

	Interval_matrix operator + (const Interval_matrix & A) const;
	Interval_matrix operator - (const Interval_matrix & A) const;
	Interval_matrix operator * (const Interval_matrix & A) const;

	Interval_matrix & operator = (const Interval_matrix & A);
};

// need a class for Taylor model matrices


void exp_taylor_remainder(Interval & result, const Interval & tmRange, const int order);
void rec_taylor_remainder(Interval & result, const Interval & tmRange, const int order);
void sin_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const int order);
void cos_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const int order);
void log_taylor_remainder(Interval & result, const Interval & tmRange, const int order);
void sqrt_taylor_remainder(Interval & result, const Interval & tmRange, const int order);

void exp_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);
void rec_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);
void sin_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);
void cos_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);
void log_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);
void sqrt_taylor_only_remainder(Interval & result, const Interval & remainder, list<Interval>::iterator & iterRange, const int order);

extern ParseSetting parseSetting;
extern ParseResult parseResult;

#endif /* TAYLORMODEL_H_ */
