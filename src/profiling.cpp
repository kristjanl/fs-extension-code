#include "profiling.h"

using namespace std;

int decreasing() {
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("x1");
	parseSetting.addVar("x2");
	parseSetting.addVar("x3");
  //need to have variables for each of the one present in component
  vector<HornerForm> odes = parseHFFromPoly("my hfs {x1 + x2 + x3, 0, 0}");
  
	parseSetting.addVar("x4");
	parseSetting.addVar("x5");
  TaylorModelVec init = parseTMV("my models {x1, x2, x3}");
  TaylorModelVec p = parseTMV("my models {x1, x2 + x3, x4}");
  
  //number of parameters on left 
  int paramCount = parseSetting.variables.size(); //substract time
  int varCount = init.tms.size();// + parsedP.tms.size();
  
  mlog1(sbuilder() << "paramCount: " << paramCount);
  mlog1(sbuilder() << "varCount: " << varCount);
  
  int order = 3;
  Interval cutoff = Interval(1e-15);
  
  MyComponent component;
  component.solveIndexes.push_back(0);
  component.initSet = init;
  component.odes = odes;
  
  
  
  mlog("initSet", component.initSet);
  mlog("odes", component.odes);
  
  for(int i = 1; i <= order; i++) {
    p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    mlog("pic_no_rem", p);
  }
  
  mlog("f_no_rem", p);
  
	vector<Interval> pPolyRange;
	vector<RangeTree *> trees;
  
  vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 
      Interval(0,0.1).sup(), 2*order);
      
  
  p.polyRangeNormal(pPolyRange, step_exp_table);
	  
    
  mlog("pRange", pPolyRange);
    
  
  
  vector<Interval> guess;
  for(int i = 0; i < varCount; i++) {
    guess.push_back(Interval(-1,1));
  }
  for(int i = 0; i < varCount; i++) {
    if(component.isSolveVar(i) == false)
      continue;
    p.tms[i].remainder = guess[i];
  }
  
  mlog("g_rem", p);
  
  TaylorModelVec compTemp;
  p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
    step_exp_table, paramCount, order, cutoff);
  mlog("final", p);
}

TaylorModelVec idModel(int dim) {
  parseSetting.clear();
	parseSetting.addVar("t");
  for(int i = 0; i < dim; i++) {
	  stringstream ss;
    ss << "x" << (i+1);
    mlog1(ss.str());
	  parseSetting.addVar(ss.str());	
	}
	
  stringstream ss;
  ss << "my models {";
  ss << "x1";
  for(int i = 1; i < dim; i++) {
    ss << ", x" << (i+1);
  }
  ss << "}";
	
	mlog1(ss.str());
  return parseTMV(ss.str());
}


TaylorModel oneModel(int dim, int params, int paramShift, int order) {
  parseSetting.clear();
	parseSetting.addVar("t");
  for(int i = 0; i < dim; i++) {
	  stringstream ss;
    ss << "x" << (i+1);
    //mlog1(ss.str());
	  parseSetting.addVar(ss.str());	
	}
	
  stringstream modelS;
  modelS << "my models {";
  modelS << "1";  
  
  //skip constant (0), include highest power
  for(int i = 1; i < pow(order + 1, params); i++) {
    stringstream ss;
    int temp = i;
    int totalOrder = 0;
    for(int param = 0; param < params; param++) {
      //power of parameter)
      int paramOrder = temp%(order + 1);
      //cout << paramOrder;
      totalOrder += paramOrder;
      if(totalOrder > order) {
        break;
      }
      
      if(param != 0) {
        ss << "*";
      }
      ss << "x" << (param + 1 + paramShift) << "^" << paramOrder;
      temp /= (order + 1);
    }
    //cout << endl;
    if(totalOrder > order)
      continue;
    modelS << " + " << ss.str(); 
    //cout << ss.str() << " - " << totalOrder << endl;
  }
  modelS << "}";
  
  //cout << modelS.str() << endl;
	
	//mlog1(ss.str());
  return parseTMV(modelS.str()).tms[0];
}


TaylorModel oneFullModel(int dim, int order) {
  return oneModel(dim, dim, 0, order);
}

int precond() {

  mlog("id", idModel(4));
  mlog("one", oneModel(3, 2, 1, 3));
  mlog("full", oneFullModel(3, 3));
  
  exit(0);
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("x1");
	parseSetting.addVar("x2");
	parseSetting.addVar("x3");
  TaylorModelVec leftToRight = parseTMV("my models {x1, 2*x2 + 13*x1, 3*x3}");
  
  
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
  int rightParamCount = parseSetting.variables.size();
  TaylorModelVec prevRight = parseTMV("my models {5*a1, 7*a2, 11*a2}");
  
  
  int order = 3;
  Interval cutoff = Interval(1e-15);
  
  
  vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 
      Interval(0,0.1).sup(), 2*order);
  
  
  
  vector<Interval> prevRightPolyRange;
	prevRight.polyRangeNormal(prevRightPolyRange, step_end_exp_table);
	
  tstart1(pre_insert);
	TaylorModelVec currentRight;
	leftToRight.insert_ctrunc_normal(currentRight, prevRight, prevRightPolyRange,
      step_end_exp_table, rightParamCount, order, cutoff);
      
  mlog("final", currentRight);
}

int main() {
  //decreasing();
  precond();
}
