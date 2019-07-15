
#include <algorithm> 
#include <utility>

#include <sys/time.h>
#include <sys/resource.h>

#include "Exceptions.h"
#include "profiling.h"

using namespace std;

#include "Utils.h"



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

//generates polynomial that has variables x{startVar} to x{endVar} (exclusive)
//with all the terms up to {order}
string pomS(int startVar, int endVar, int order) {
  bool notFirst = false;
  //skip constant (0), include highest power
  int varCount = endVar - startVar;
  //mlog1(sbuilder() << "varCount: " << varCount);
  
  stringstream ss;
  //skip constant by starting from 1
  for(int i = 1; i < pow(order + 1, varCount); i++) {
    //mlog1(sbuilder() << "i: " << i);
    int temp = i;
    stringstream mono;
    int totalOrder = 0;
    for(int varShift = 0; varShift < varCount; varShift++) {
      int varOrder = temp%(order + 1);
      if(varShift != 0)
        mono << "*";
      mono << "x" << (startVar + varShift) << "^" << varOrder;
      totalOrder += varOrder;
      temp /= order + 1;

      if(totalOrder > order) {
        break;
      }
    }
    if(totalOrder > order) {
      continue;
    }
    if(i!= 1)
      ss << " + ";
    ss << mono.str();
    //cout << mono.str() << endl;
  }
  //cout << ss.str() << endl;
  return ss.str();
}

int polyPart() {
  cout << "polyPart" << endl;
  int order = 7;

  for(int k = 0; k < 11; k++) {
    clock_t start2 = clock();
    int varCount = k;
    if(varCount == 0)
      varCount = 3;
    else
      varCount *= 10;

    //int solveVar = varCount;
    int solveVar = 1;
    int paramCount = varCount + 1;
    
    parseSetting.clear();
	  parseSetting.addVar("t");
	  //cout << i << endl;
	  for(int j = 1; j < varCount + 1; j++) {
      stringstream ss;
      ss << "x" << j;
  	  parseSetting.addVar(ss.str());
	  }
	  
    stringstream hss; hss << "my hfs {";
    stringstream iss; iss << "my models {";
    stringstream pss; pss << "my models {";
    for(int j = 1; j < varCount + 1; j++) {
      //cout << "var: " << j << endl;
      stringstream ss;

      if(j != 1) {hss << ","; iss << ","; pss << ",";}

      //hss << "x" << j << " + x" << j << "^2";
      hss << pomS(j, j+1, 1);
      iss << pomS(j, j+1, order);
      pss << pomS(j, j+1, order);

      if(j == solveVar && varCount > 2) {
        hss << " + " << pomS(j, j+3, order);
        //hss << " + x" << j-1 << " + x" << j-1 << "^2";
        //hss << " + x" << j-2 << " + x" << j-2 << "^2";
        //iss << " + " << pomS(j-1, j+1, 5);
        //pss << " + " << pomS(j-1, j+1, 5);
      }
	  }
    hss << "}"; iss << "}"; pss << "}";
	  //cout << hss.str() << endl;
	  //cout << iss.str() << endl;
	  //cout << pss.str() << endl;
	  
    //need to have variables for each of the one present in component
    //vector<HornerForm> odes = parseHFFromPoly("my hfs {x1 + x1^2}");
    vector<HornerForm> odes = parseHFFromPoly(hss.str());
    
    //TaylorModelVec init = parseTMV("my models {x1}");
    //TaylorModelVec p = parseTMV("my models {x1}");
    TaylorModelVec init = parseTMV(iss.str());
    TaylorModelVec p = parseTMV(pss.str());

        
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(solveVar - 1);
    component.initSet = init;
    component.odes = odes;
    
    clock_t start = clock();
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
    cout << varCount << "\t" << (long double)(clock() - start) / CLOCKS_PER_SEC << endl;
    //mlog("pic_no_rem", p);
    //cout << "    " << (long double)(clock() - start2) / CLOCKS_PER_SEC << endl;
  }
}



class Se {
  public:
    Se(int order, bool first, int hfVars, int hfOrder, int iOrder, int pOrder);
    Se(int order, bool first, int hfVars, int hfOrder, int iOrder, int pOrder, int refCount);
    int order;
    bool first;
    int hfVars;
    int hfOrder;
    int iOrder;
    int pOrder;
    int refCount;
    string toString();
};

Se::Se(int order, bool first, int hfVars, int hfOrder, int iOrder, int pOrder) 
    : order(order), first(first), hfVars(hfVars), hfOrder(hfOrder), iOrder(iOrder)
    , pOrder(pOrder), refCount(-1) {
}

Se::Se(int order, bool first, int hfVars, int hfOrder, int iOrder, int pOrder, int refCount) 
    : order(order), first(first), hfVars(hfVars), hfOrder(hfOrder), iOrder(iOrder)
    , pOrder(pOrder), refCount(refCount) {
}

string Se::toString() {
  stringstream ss;
  ss << "order: " << order << ", var: " << (first?"first":"last") << ", hfVars: " << hfVars 
      << ", hfOrder: " << hfOrder << ", iOrder: " << iOrder << ", pOrder: " << pOrder;
  if(refCount != -1)
    ss << ", refCount: " << refCount;
  return ss.str();
}


typedef pair<string, vector<long double> > se_data;

extern bool stop;

se_data polyPart(Se se) {
  cout << "polyPart" << endl;
  int order = se.order;

  se_data ret;
  ret.first = se.toString();

  for(int k = 0; k < 11; k++) {
    clock_t start2 = clock();
    int varCount = k;
    if(varCount == 0)
      varCount = 3;
    else
      varCount *= 10;

    int solveVar = -1;
    if(se.first) {
      solveVar = 1;
    } else {
      solveVar = varCount;
    }
    int paramCount = varCount + 1;
    
    parseSetting.clear();
	  parseSetting.addVar("t");
	  //cout << i << endl;
	  for(int j = 1; j < varCount + 1; j++) {
      stringstream ss;
      ss << "x" << j;
  	  parseSetting.addVar(ss.str());
	  }
	  
    stringstream hss; hss << "my hfs {";
    stringstream iss; iss << "my models {";
    stringstream pss; pss << "my models {";
    for(int j = 1; j < varCount + 1; j++) {
      //cout << "var: " << j << endl;
      stringstream ss;

      if(j != 1) {hss << ","; iss << ","; pss << ",";}

      //hss << "x" << j << " + x" << j << "^2";
      //hss << pomS(j, j+1, 1);
      iss << pomS(j, j+1, se.iOrder);
      pss << pomS(j, j+1, se.pOrder);

      if(j == solveVar) {
        if(se.first) {
          hss << pomS(j, j + se.hfVars, se.hfOrder);
        } else {
          hss << pomS(j + 1 - se.hfVars, j+1, se.hfOrder);
        }
        //hss << " + x" << j-1 << " + x" << j-1 << "^2";
        //hss << " + x" << j-2 << " + x" << j-2 << "^2";
        //iss << " + " << pomS(j-1, j+1, 5);
        //pss << " + " << pomS(j-1, j+1, 5);
      } else {
        hss << pomS(j, j+1, 1);
      }
	  }
    hss << "}"; iss << "}"; pss << "}";
	  //cout << hss.str() << endl;
	  //cout << iss.str() << endl;
	  //cout << pss.str() << endl;
	  
    //need to have variables for each of the one present in component
    vector<HornerForm> odes = parseHFFromPoly(hss.str());
    TaylorModelVec init = parseTMV(iss.str());
    TaylorModelVec p = parseTMV(pss.str());

    if(stop) {
      mlog("odes", odes);
      mlog("init", init);
      mlog("p", p);
      exit(0);
    }


        
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(solveVar - 1);
    component.initSet = init;
    component.odes = odes;
    
    clock_t start = clock();
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
    cout << varCount << "\t" << (long double)(clock() - start) / CLOCKS_PER_SEC << endl;
    ret.second.push_back((long double)(clock() - start) / CLOCKS_PER_SEC);
    //mlog("pic_no_rem", p);
    //cout << "    " << (long double)(clock() - start2) / CLOCKS_PER_SEC << endl;
  }
  return ret;
}





void applyAndPrintSet(vector<Se> set, se_data (*f)(Se)) {
  vector<se_data> res;
  for(int i = 0; i < set.size(); i++) {
    res.push_back(f(set[i]));
  }

  cout << "# var";
  for(int j = 0; j < res.size(); j++) {
    cout << "\t" << res[j].first;
  }
  cout << endl;
  for(int i = 0; i < res[0].second.size(); i++) {
    int varCount = i;
    if(varCount == 0)
      varCount = 3;
    else
      varCount *= 10;
    cout << varCount;
    for(int j = 0; j < res.size(); j++) {
      cout << "\t" << res[j].second[i];
    }
    cout << endl;
  }
}
bool stop = false;
void poly() {
  //order, first, hfVars, hfOrder, iOrder, pOrder
  //Se(7, true, 3, 7, 7, 7)

  vector<Se> set1;
  set1.push_back(Se(30, true, 1, 2, 1, 1));
  set1.push_back(Se(15, true, 1, 15, 1, 1));
  set1.push_back(Se(15, true, 1, 15, 15, 15));
  applyAndPrintSet(set1, polyPart);

  vector<Se> set2;
  set2.push_back(Se(30, false, 1, 2, 1, 1));
  set2.push_back(Se(15, false, 1, 15, 1, 1));
  set2.push_back(Se(15, false, 1, 15, 15, 15));
  //applyAndPrintSet(set2, polyPart);

  vector<Se> set3;
  set3.push_back(Se(15, true, 2, 2, 1, 1));
  set3.push_back(Se(10, true, 2, 10, 1, 1));
  set3.push_back(Se(10, true, 3, 2, 1, 1));
  set3.push_back(Se(7, true, 3, 7, 1, 1));
  //applyAndPrintSet(set3, polyPart);

  vector<Se> set4;
  set4.push_back(Se(15, false, 2, 2, 1, 1));
  set4.push_back(Se(10, false, 2, 10, 1, 1));
  set4.push_back(Se(10, false, 3, 2, 1, 1));
  set4.push_back(Se(7, false, 3, 7, 1, 1));
  applyAndPrintSet(set4, polyPart);
}



se_data treePart(Se se) {
  cout << "treePart" << endl;
  int order = se.order;

  se_data ret;
  ret.first = se.toString();

  for(int k = 0; k < 11; k++) {
    clock_t start2 = clock();
    int varCount = k;
    if(varCount == 0)
      varCount = 3;
    else
      varCount *= 10;

    int solveVar = -1;
    if(se.first) {
      solveVar = 1;
    } else {
      solveVar = varCount;
    }
    int paramCount = varCount + 1;
    
    parseSetting.clear();
	  parseSetting.addVar("t");
	  //cout << i << endl;
	  for(int j = 1; j < varCount + 1; j++) {
      stringstream ss;
      ss << "x" << j;
  	  parseSetting.addVar(ss.str());
	  }
	  
    stringstream hss; hss << "my hfs {";
    stringstream iss; iss << "my models {";
    stringstream pss; pss << "my models {";
    for(int j = 1; j < varCount + 1; j++) {
      //cout << "var: " << j << endl;
      stringstream ss;

      if(j != 1) {hss << ","; iss << ","; pss << ",";}

      //hss << "x" << j << " + x" << j << "^2";
      //hss << pomS(j, j+1, 1);
      iss << pomS(j, j+1, se.iOrder);
      pss << pomS(j, j+1, se.pOrder);

      if(j == solveVar) {
        if(se.first) {
          hss << pomS(j, j + se.hfVars, se.hfOrder);
        } else {
          hss << pomS(j + 1 - se.hfVars, j+1, se.hfOrder);
        }
        //hss << " + x" << j-1 << " + x" << j-1 << "^2";
        //hss << " + x" << j-2 << " + x" << j-2 << "^2";
        //iss << " + " << pomS(j-1, j+1, 5);
        //pss << " + " << pomS(j-1, j+1, 5);
      } else {
        hss << pomS(j, j+1, 1);
      }
	  }
    hss << "}"; iss << "}"; pss << "}";
	  //cout << hss.str() << endl;
	  //cout << iss.str() << endl;
	  //cout << pss.str() << endl;
	  
    //need to have variables for each of the one present in component
    vector<HornerForm> odes = parseHFFromPoly(hss.str());
    TaylorModelVec init = parseTMV(iss.str());
    TaylorModelVec p = parseTMV(pss.str());

    if(stop) {
      mlog("odes", odes);
      mlog("init", init);
      mlog("p", p);
      exit(0);
    }


        
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(solveVar - 1);
    component.initSet = init;
    component.odes = odes;
    
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
    
    vector<Interval> step_exp_table, step_end_exp_table;
	  construct_step_exp_table(step_exp_table, step_end_exp_table, 
        Interval(0,0.1).sup(), 2*order);
        
    
    p.polyRangeNormal(pPolyRange, step_exp_table);

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

    
		throw runtime_error("fix it");
    /*TODO fix it
    p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
      step_exp_table, paramCount, order, cutoff);
    clock_t end = clock();
    cout << varCount << "\t" << (long double)(end - start) / CLOCKS_PER_SEC << endl;
    ret.second.push_back((long double)(end - start) / CLOCKS_PER_SEC);
    */
  }
  return ret;
}

void tree() {
  //order, first, hfVars, hfOrder, iOrder, pOrder
  //Se(7, true, 3, 7, 7, 7)

  vector<Se> set1;
  set1.push_back(Se(30, true, 1, 2, 1, 1));
  set1.push_back(Se(15, true, 1, 15, 1, 1));
  set1.push_back(Se(15, true, 1, 15, 15, 15));
  //applyAndPrintSet(set1, treePart);

  vector<Se> set2;
  set2.push_back(Se(30, false, 1, 2, 1, 1));
  set2.push_back(Se(15, false, 1, 15, 1, 1));
  set2.push_back(Se(15, false, 1, 15, 15, 15));
  //applyAndPrintSet(set2, treePart);

  vector<Se> set3;
  set3.push_back(Se(15, true, 2, 2, 1, 1));
  set3.push_back(Se(10, true, 2, 10, 1, 1));
  set3.push_back(Se(10, true, 3, 2, 1, 1));
  set3.push_back(Se(7, true, 3, 7, 1, 1));
  //applyAndPrintSet(set3, treePart);

  vector<Se> set4;
  set4.push_back(Se(15, false, 2, 2, 1, 1));
  set4.push_back(Se(10, false, 2, 10, 1, 1));
  set4.push_back(Se(10, false, 3, 2, 1, 1));
  set4.push_back(Se(7, false, 3, 7, 1, 1));
  //applyAndPrintSet(set4, treePart);
}


se_data refPart(Se se) {
  cout << "refPart" << endl;
  int order = se.order;

  se_data ret;
  ret.first = se.toString();

  for(int k = 0; k < 11; k++) {
    clock_t start2 = clock();
    int varCount = k;
    if(varCount == 0)
      varCount = 3;
    else
      varCount *= 10;

    int solveVar = -1;
    if(se.first) {
      solveVar = 1;
    } else {
      solveVar = varCount;
    }
    int paramCount = varCount + 1;
    
    parseSetting.clear();
	  parseSetting.addVar("t");
	  //cout << i << endl;
	  for(int j = 1; j < varCount + 1; j++) {
      stringstream ss;
      ss << "x" << j;
  	  parseSetting.addVar(ss.str());
	  }
	  
    stringstream hss; hss << "my hfs {";
    stringstream iss; iss << "my models {";
    stringstream pss; pss << "my models {";
    for(int j = 1; j < varCount + 1; j++) {
      //cout << "var: " << j << endl;
      stringstream ss;

      if(j != 1) {hss << ","; iss << ","; pss << ",";}

      //hss << "x" << j << " + x" << j << "^2";
      //hss << pomS(j, j+1, 1);
      iss << pomS(j, j+1, se.iOrder);
      pss << pomS(j, j+1, se.pOrder);

      if(j == solveVar) {
        if(se.first) {
          hss << pomS(j, j + se.hfVars, se.hfOrder);
        } else {
          hss << pomS(j + 1 - se.hfVars, j+1, se.hfOrder);
        }
        //hss << " + x" << j-1 << " + x" << j-1 << "^2";
        //hss << " + x" << j-2 << " + x" << j-2 << "^2";
        //iss << " + " << pomS(j-1, j+1, 5);
        //pss << " + " << pomS(j-1, j+1, 5);
      } else {
        hss << pomS(j, j+1, 1);
      }
	  }
    hss << "}"; iss << "}"; pss << "}";
	  //cout << hss.str() << endl;
	  //cout << iss.str() << endl;
	  //cout << pss.str() << endl;
	  
    //need to have variables for each of the one present in component
    vector<HornerForm> odes = parseHFFromPoly(hss.str());
    TaylorModelVec init = parseTMV(iss.str());
    TaylorModelVec p = parseTMV(pss.str());

    if(stop) {
      mlog("odes", odes);
      mlog("init", init);
      mlog("p", p);
      exit(0);
    }


        
    Interval cutoff = Interval(1e-15);
    
    MyComponent component;
    component.solveIndexes.push_back(solveVar - 1);
    component.initSet = init;
    component.odes = odes;
    
    for(int i = 1; i <= order; i++) {
      p.Picard_no_remainder_assign(&component, paramCount, i, cutoff);
    }
    
	  vector<Interval> pPolyRange;
	  vector<RangeTree *> trees;
    
    vector<Interval> step_exp_table, step_end_exp_table;
	  construct_step_exp_table(step_exp_table, step_end_exp_table, 
        Interval(0,0.1).sup(), 2*order);
        
    
    p.polyRangeNormal(pPolyRange, step_exp_table);

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
    
		throw runtime_error("fix it");
    //p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
    //  step_exp_table, paramCount, order, cutoff);
    
	  vector<Interval> cutoffInt;
    for(int i=0; i < varCount; i++) {
		  Polynomial polyTemp;
		  polyTemp = compTemp.tms[i].expansion - p.tms[i].expansion;

		  Interval intTemp;
		  polyTemp.intEvalNormal(intTemp, step_exp_table);
		  
		  cutoffInt.push_back(intTemp);
		  
      compTemp.tms[i].remainder += intTemp;
	  }


    clock_t start = clock();
		vector<Interval> newRemainders;
    if(se.refCount == -1) {
      throw IntegrationException(sbuilder() << "refCount not set");
    }

	  for(int l = 0; l < se.refCount; l++) {
		  for(int i = 0; i < varCount; i++) {        
        if(component.isSolveVar(i) == false)
          continue;
		    p.tms[i].remainder = guess[i];
		  }
  		p.Picard_only_remainder(newRemainders, trees, &component, step_exp_table[1]);
		}
    clock_t end = clock();
    cout << varCount << "\t" << (long double)(end - start) / CLOCKS_PER_SEC << endl;
    ret.second.push_back((long double)(end - start) / CLOCKS_PER_SEC);
  }
  return ret;
}


void refining() {
  //order, first, hfVars, hfOrder, iOrder, pOrder, refCount
  //Se(7, true, 3, 7, 7, 7, 100)

  vector<Se> set1;
  set1.push_back(Se(30, true, 1, 2, 1, 1, 100));
  set1.push_back(Se(15, true, 1, 15, 1, 1, 100));
  set1.push_back(Se(15, true, 1, 15, 15, 15, 100));
  //applyAndPrintSet(set1, refPart);

  vector<Se> set2;
  set2.push_back(Se(30, false, 1, 2, 1, 1, 100));
  set2.push_back(Se(15, false, 1, 15, 1, 1, 100));
  set2.push_back(Se(15, false, 1, 15, 15, 15, 100));
  //applyAndPrintSet(set2, refPart);

  vector<Se> set3;
  set3.push_back(Se(15, true, 2, 2, 1, 1, 100));
  set3.push_back(Se(10, true, 2, 10, 1, 1, 10));
  set3.push_back(Se(10, true, 3, 2, 1, 1, 100));
  set3.push_back(Se(7, true, 3, 7, 1, 1, 10));
  //applyAndPrintSet(set3, refPart);

  vector<Se> set4;
  set4.push_back(Se(15, false, 2, 2, 1, 1, 100));
  set4.push_back(Se(10, false, 2, 10, 1, 1, 10));
  set4.push_back(Se(10, false, 3, 2, 1, 1, 100));
  set4.push_back(Se(7, false, 3, 7, 1, 1, 10));
  applyAndPrintSet(set4, refPart);
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
    
		throw runtime_error("fix it");
    //p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
    //  step_exp_table, paramCount, order, cutoff);
  	clock_t end = clock();
    cout << k << ": " << (long double)(end - start) / CLOCKS_PER_SEC << endl;
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
  throw runtime_error("fix it");
  //p.Picard_ctrunc_normal(compTemp, trees, &component, pPolyRange, 
  //  step_exp_table, paramCount, order, cutoff);
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
    cout << i << ": " << (long double)(end - start) / CLOCKS_PER_SEC << endl;
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

void andOr() {
  cout << "here\n";
}

void foo(Se se) {
  cout << "foo: " << se.toString() << endl; 
}





Monomial getMono(int dim, bool first=true) {
  vector<int> d;
  if(first)
    d.push_back(2);
  for(int j = 0; j < dim - 1; j++) {
    d.push_back(0);
  }
  if(first == false)
    d.push_back(2);
  return Monomial(Interval(2), d);
}


//#include <fenv.h>
void monoArith() {
  for(int step = 0; step < 30 + 1; step += 1) {
    int dim = step * 20;
    if(dim == 0)
      dim = 1;
    Monomial m1 = getMono(dim - 1, false);
    Monomial m2 = getMono(dim - 1, false);

    cout << dim;

    long t1, t2;

    /*
    t1 = getTime();
    for(int i = 0; i < 10000; i++) {
      Monomial m3 = m1 * m2;
      //m1 *= m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false);
    */
    t1 = getTime();
    for(int i = 0; i < 100000; i++) {
      Monomial m3 = m1 * m2;
      //m1 *= m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false) << endl;
  }
}

void monoOperations() {
  int count = 100000;
  for(int step = 0; step < 30 + 1; step += 1) {
    int dim = step * 20;
    if(dim == 0)
      dim = 1;
    Monomial m1 = getMono(dim - 1, false);
    Monomial m2 = getMono(dim - 1, false);

    cout << dim;

    long t1, t2;

    t1 = getTime();
    for(int i = 0; i < count; i++) {
      Monomial m3 = m1 + m2;
      //m1 *= m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false);
    t1 = getTime();
    for(int i = 0; i < count; i++) {
      Monomial m3 = m1 * m2;
      //m1 *= m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false);


    t1 = getTime();
    for(int i = 0; i < count; i++) {
      m1 += m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false);


    t1 = getTime();
    for(int i = 0; i < count; i++) {
      m1 *= m2;
    }
    t2 = getTime();
    cout << "\t " << getDiff(t1, t2, false);

    cout << endl;
  }
}

void vectorSize() {
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("x1");
  TaylorModelVec init = parseTMV("my models {x1 + x1^2 + 1^3 + x1^4}");

  mlog("init", init);
  cout << init.toString(parseSetting.variables) << endl;

  cout << "here" << endl;

  int n = 20;
  long t1,t2;
  t1 = getTime();
  for(int j = 0; j < 100; j++) {
    TaylorModelVec tmv(n);

    for(int i = 0; i < n; i++) {
      //tmv.tms.push_back(init.tms[0]);r
      tmv.tms[i] = init.tms[0];
    }
  }
  t2 = getTime();
  cout << "time: " << getDiff(t1, t2, false) << endl;

  t1 = getTime();
  for(int j = 0; j < 100; j++) {
    TaylorModelVec tmv;
    tmv.tms.reserve(n);
    for(int i = 0; i < n; i++) {
      tmv.tms.push_back(init.tms[0]);
    }
  }
  t2 = getTime();
  cout << "time: " << getDiff(t1, t2, false) << endl;
  
  t1 = getTime();
  for(int j = 0; j < 100; j++) {
    TaylorModelVec tmv;
    for(int i = 0; i < n; i++) {
      tmv.tms.push_back(init.tms[0]);
    }
  }
  t2 = getTime();
  cout << "time: " << getDiff(t1, t2, false) << endl;
  
}

int main() {
  //vectorSize();

  /*
  vector<int> v;
  v.reserve(50);
  v.clear();
  cout << v.capacity() << endl;
  */

  //monoArith();
  monoOperations();

  //polyPart();
  //poly();
  //tree();
  //refining();
  //decPart();
  //decreasing();
  //precond();
  //precondMultPossibleOnePresentParameter();
  //andOr();
}
