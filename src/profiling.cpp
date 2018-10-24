
#include <algorithm> 

#include "profiling.h"



using namespace std;

TaylorModelVec idModel(int paramCount) {
  parseSetting.clear();
	parseSetting.addVar("t");
  for(int i = 0; i < paramCount; i++) {
	  stringstream ss;
    ss << "x" << (i+1);
    //mlog1(ss.str());
	  parseSetting.addVar(ss.str());	
	}
	
  stringstream ss;
  ss << "my models {";
  ss << "x1";
  for(int i = 1; i < paramCount; i++) {
    ss << ", x" << (i+1);
  }
  ss << "}";
	
	//mlog1(ss.str());
  return parseTMV(ss.str());
}


//makes 1-dimensional TM with paramCount parameters, paramsPresent parameters
//will have nonzeroe coefficents (all monomials possible will be present)
//paramShift skips the first parameters
//order denotes the order of monomials present
//(3,1,0,2) -> 1 + a1 + a1^2
//(3,1,1,2) -> 1 + a2 + a2^2
//(3,2,1,2) -> 1 + a1 + a1^2 + a2 + a1*a2 + a2^2
TaylorModel oneModel(int paramCount, int paramsPresent, int paramShift, 
    int order) {
  parseSetting.clear();
	parseSetting.addVar("t");
  for(int i = 0; i < paramCount; i++) {
	  stringstream ss;
    ss << "a" << (i+1);
    //mlog1(ss.str());
	  parseSetting.addVar(ss.str());	
	}
	
  stringstream modelS;
  modelS << "my models {";
  
  bool notFirst = false;
  //skip constant (0), include highest power
  for(int i = 1; i < pow(order + 1, paramsPresent); i++) {
    stringstream ss;
    int temp = i;
    int totalOrder = 0;
    for(int param = 0; param < paramsPresent; param++) {
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
      ss << "a" << (param + 1 + paramShift) << "^" << paramOrder;
      temp /= (order + 1);
    }
    //cout << endl;
    if(totalOrder > order)
      continue;
    if(notFirst)
      modelS << " + ";
    modelS << "5e-2*" << ss.str();
    notFirst = true;
    //cout << ss.str() << " - " << totalOrder << endl;
  }
  modelS << "}";
  
  //cout << modelS.str() << endl;
	
	//mlog1(ss.str());
  return parseTMV(modelS.str()).tms[0];
}


TaylorModel oneFullModel(int paramCount, int order) {
  return oneModel(paramCount, paramCount, 0, order);
}



int polyPart() {
  cout << "polyPart" << endl;
  for(int k = 1; k < 100; k++) {
    parseSetting.clear();
	  parseSetting.addVar("t");
	  //cout << i << endl;
	  for(int j = 1; j < k + 1; j++) {
      stringstream ss;
      ss << "x" << j;
      //cout << ss.str() << endl;
  	  parseSetting.addVar(ss.str());
	  }
	  
	  
    //need to have variables for each of the one present in component
    vector<HornerForm> odes = parseHFFromPoly("my hfs {x1 + x1^2}");
    
    TaylorModelVec init = parseTMV("my models {x1}");
    TaylorModelVec p = parseTMV("my models {x1}");
    
    //number of parameters on left 
    int paramCount = parseSetting.variables.size(); //substract time
    int varCount = init.tms.size();
    
    //mlog1(sbuilder() << "paramCount: " << paramCount);
    //mlog1(sbuilder() << "varCount: " << varCount);
    
    int order = 20;
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(0);
    component.initSet = init;
    component.odes = odes;
    
    
    
    //mlog("initSet", component.initSet);
    //mlog("odes", component.odes);
    
    
    clock_t start = clock();
    
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
  	clock_t end = clock();
    cout << k << ": " << double(end - start) / CLOCKS_PER_SEC << endl;
    //mlog("pic_no_rem", p);
  }
}

int decPart() {
  cout << "decPart" << endl;
  
  for(int k = 1; k < 100; k++) {
    parseSetting.clear();
	  parseSetting.addVar("t");  
    stringstream pStr;
    pStr << "my models{x1";
    
	  //cout << i << endl;
	  for(int j = 1; j < k + 1; j++) {
      stringstream ss;
      ss << "x" << j;
      //cout << ss.str() << endl;
  	  parseSetting.addVar(ss.str());
  	  if(j!=1)
  	    pStr << ",0";
	  }
	  pStr << "}";
	  //cout << pStr.str() << endl;
	  
    //need to have variables for each of the one present in component
    vector<HornerForm> odes = parseHFFromPoly("my hfs {x1}");
    
    TaylorModelVec init = parseTMV("my models {x1}");
    TaylorModelVec p = parseTMV(pStr.str());
    
    //number of parameters on left 
    int paramCount = parseSetting.variables.size(); //substract time
    int varCount = init.tms.size();
    
    //mlog1(sbuilder() << "paramCount: " << paramCount);
    //mlog1(sbuilder() << "varCount: " << varCount);
    
    int order = 200;
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(0);
    component.initSet = init;
    component.odes = odes;
    
    
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
    
    //mlog("f_no_rem", p);
  
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
    
    vector<Interval> step_exp_table, step_end_exp_table;
	  construct_step_exp_table(step_exp_table, step_end_exp_table, 
        Interval(0,0.1).sup(), 2*order);
        
    
    p.polyRangeNormal(pPolyRange, step_exp_table);
	    
      
    //mlog("pRange", pPolyRange);
      
    
    
    vector<Interval> guess;
    for(int i = 0; i < varCount; i++) {
      guess.push_back(Interval(-1,1));
    }
    for(int i = 0; i < varCount; i++) {
      if(component.isSolveVar(i) == false)
        continue;
      p.tms[i].remainder = guess[i];
    }
    
    //mlog("g_rem", p);
    
    TaylorModelVec compTemp;
    
    clock_t start = clock();
    p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
      step_exp_table, paramCount, order, cutoff);
  	clock_t end = clock();
    cout << k << ": " << double(end - start) / CLOCKS_PER_SEC << endl;
  }
}



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

void precondMultPossibleOnePresentParameter() {
  int order = 35;
  Interval cutoff = Interval(1e-15);
  
  vector<Interval> step_exp_table, step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 
      Interval(0,0.1).sup(), 2*order);
  
  
  for(int i = 99; i < 100; i++) {
    int possibleParams = i;
    TaylorModelVec leftToRight;
    leftToRight.tms.push_back(oneModel(possibleParams, 1, 0, order));
    
    TaylorModelVec prevRight;
    prevRight.tms.push_back(oneModel(possibleParams, 1, 0, order));
    
    int rightParamCount = prevRight.getParamCount();
        
    vector<Interval> prevRightPolyRange;
	  prevRight.polyRangeNormal(prevRightPolyRange, step_end_exp_table);
	
	
  	clock_t start = clock();
    TaylorModelVec currentRight;
	  leftToRight.insert_ctrunc_normal(currentRight, prevRight, prevRightPolyRange,
        step_end_exp_table, rightParamCount, order, cutoff);
  	clock_t end = clock();
    cout << i << ": " << double(end - start) / CLOCKS_PER_SEC << endl;
  }  
  
  
}



void precond() {
  /*
  mlog("id", idModel(4));
  mlog("one", oneModel(3, 2, 1, 3));
  mlog("full", oneFullModel(3, 3));
  //*/
  
  int tmOrder = 6;
  
  TaylorModelVec tmv;
  tmv.tms.push_back(oneModel(3, 1, 0, tmOrder));
  tmv.tms.push_back(oneModel(3, 2, 0, tmOrder));
  tmv.tms.push_back(oneModel(3, 3, 0, tmOrder));
  
  
  TaylorModelVec leftToRight;
  leftToRight.tms.push_back(oneModel(3,3,0,tmOrder));
  //leftToRight = idModel(3);
  
  TaylorModelVec prevRight;
  prevRight = tmv;
  int rightParamCount = prevRight.getParamCount();
  
  mlog("l", leftToRight);
  mlog("r", prevRight);
  
  
  /*
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
	parseSetting.addVar("a3");
  TaylorModelVec leftToRight = parseTMV("my models {a1, 2*a2 + 13*a1, 3*a3}");
  
  
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
  int rightParamCount = parseSetting.variables.size();
  TaylorModelVec prevRight = parseTMV("my models {5*a1, 7*a2, 11*a2}");
  //*/
  
  
  int order = tmOrder;
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
  //polyPart();
  //decPart();
  //decreasing();
  //precond();
  precondMultPossibleOnePresentParameter();
}
