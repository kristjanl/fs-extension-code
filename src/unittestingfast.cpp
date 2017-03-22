#include "unittesting.h"

using namespace std;

vector<string> ugetVNames(int n) {
  vector<string> vars;
	vars.push_back("t");
	char name[10];
	for(int i=0; i<n; ++i) {
		sprintf(name, "%s%d", "a", i+1);
    vars.push_back(name);
	}
  return vars;
}

void self() {
  vector<int> d;
  d.push_back(2);
	Interval I(1);
	Monomial m(I, d);
  cout << m.equals(m) << endl;
}

void eq() {
  vector<int> d1;
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(1);
	Interval I2(1);
	Monomial m2(I2, d2);
  cout << m1.equals(m2) << endl;
}

void difs() {
  vector<int> d1;
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(1);
  d2.push_back(1);
	Interval I2(1);
	Monomial m2(I2, d2);
  int r = m1.equals(m2);
  
  cout << r << endl;
}

void difd() {
  vector<int> d1;
  d1.push_back(1);
  d1.push_back(1);
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(1);
  d2.push_back(1);
  d2.push_back(2);
	Interval I2(1);
	Monomial m2(I2, d2);
  int r = m1.equals(m2);
  
  cout << r << endl;
}

void transf() {
  vector<int> d1;
  d1.push_back(0);
  d1.push_back(0);
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(0);
  d2.push_back(1);
	Interval I2(1);
	Monomial m2(I2, d2);
  
  
  vector<int> map;
  map.push_back(1);
  Monomial m3 = m1.transform(map);  
  int r = m2.equals(m3);
  
  
  cout << r << endl;
  
}


void transf_padding() {
  vector<int> d1;
  d1.push_back(1);
  d1.push_back(2);
  d1.push_back(3);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(1);
  d2.push_back(0);
  d2.push_back(2);
  d2.push_back(0);
  d2.push_back(3);
  d2.push_back(0);
	Interval I2(1);
	Monomial m2(I2, d2);
  
  
  vector<int> map;
  map.push_back(PADDING_VARIABLE);
  map.push_back(0);
  map.push_back(PADDING_VARIABLE);
  map.push_back(1);
  map.push_back(PADDING_VARIABLE);
  Monomial m3 = m1.transform(map);  
  int r = m2.equals(m3);
  
  //cout << m1.toString(ugetVNames(10)) << endl;
  //cout << m3.toString(ugetVNames(10)) << endl;
  //cout << m2.toString(ugetVNames(10)) << endl;
  cout << r << endl;
}

void transfe() {
  vector<int> d1;
  d1.push_back(1);
  d1.push_back(2);
  d1.push_back(3);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(0);
  d2.push_back(1);
	Interval I2(1);
	Monomial m2(I2, d2);
  
  vector<int> map;
  map.push_back(0);
  
  try {
    Monomial m3 = m1.transform(map);
    cout << "didn't fail" << endl;
  } catch (const std::invalid_argument& e) {
    cout << "caught" << endl;
  }
}


void poly() {
  vector<int> d1;
  d1.push_back(0);
  d1.push_back(0);
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(0);
  d2.push_back(1);
  d2.push_back(0);
	Interval I2(2);
	Monomial m2(I2, d2);
  
  list<Monomial> ms;
  ms.push_back(m1);
  ms.push_back(m2);
  Polynomial p(ms);
  cout << p.toString(ugetVNames(3)) << endl;
}



void component() {
  MyComponent c1 = MyComponent();
  c1.varIndexes.push_back(0);
  c1.varIndexes.push_back(3);
  c1.varIndexes.push_back(5);
  MyComponent c2 = MyComponent();
  c2.varIndexes.push_back(2);
  c2.varIndexes.push_back(4);
  c2.varIndexes.push_back(5);
  c2.varIndexes.push_back(6);
  MyComponent c = MyComponent();
  c.previous.push_back(c1);
  c.previous.push_back(c2);
  c1.log();
  c2.log();
  c.log();
  vector< vector<int> > mappers = c.previousMappers();
  vector< vector<int> > expected;
  vector<int> e1;
  vector<int> e2;
  e1.push_back(0);
  e1.push_back(-2);
  e1.push_back(1);
  e1.push_back(-2);
  e1.push_back(2);
  e1.push_back(-2);
  
  e2.push_back(-2);
  e2.push_back(0);
  e2.push_back(-2);
  e2.push_back(1);
  e2.push_back(2);
  e2.push_back(3);
  expected.push_back(e1);
  expected.push_back(e2);
  
  for(int i = 0; i < expected.size(); i++) {
    vector<int> v1 = expected.at(i);
    vector<int> v2 = mappers.at(i);
    //logger.logVi("v1", v1);
    //logger.logVi("m2", v2);
    //logger.log(v1.size());
    //logger.log(v2.size());
    for(int j = 0; j < v1.size(); j++) {
      logger.log(sbuilder() << "v1(" << j << ") = " << v1.at(j) << ", m2(" << j << ") = " << v2.at(j));
    }
    if(equal(v1.begin(), v1.end(), v2.begin())) {
      cout << "eq" << endl;
      cout << "eq: " << equal(v1.begin(), v1.end(), v2.begin()) << endl;
    }
    else {
      cout << "ne" << endl;
      cout << "ne: " << equal(v1.begin(), v1.end(), v2.begin()) << endl;
    }
    cout << equal(v1.begin(), v1.end(), v2.begin()) << endl;
  }
  cout << "here" << endl;
}

void sw1() {
  MyComponent component;
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  vector<Interval> step_exp_table;
  vector<Interval> step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 0.1, 2*2);
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parse("my models {2 + 4*a + 0.5 * a^2 + [-0.2,0.2],1 + 3*b + a * b + [-0.1,0.1]}");
  component.pipes.push_back(parseResult.tmv);
	
  double q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  
  cout << q << endl;
}

void sw2() {
  MyComponent component;
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  vector<Interval> step_exp_table;
  vector<Interval> step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 0.1, 2*2);
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parse("my models {a + [-1,1],2*b + [-2,2]}");
	logger.logTMV("res", parseResult.tmv);
	component.pipes.push_back(parseResult.tmv);
	
  double q;
  q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  cout << q << endl;
}

void swSet() {
  MyComponent component;
  vector<Interval> domain;
	component.pipes.push_back(parseResult.tmv);
	component.compVars.push_back(0);
	component.allTMParams.push_back(0);
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parse("my models {3 + a + [-1,1]}");
	component.swInput = parseResult.tmv;
	logger.logTMV("tmv1", component.initSet);
	smallComp::shrinkWrapSet(component, &component, 2.0, domain);
  TaylorModelVec expected = parseTMV("my models {3 + 2 * a}");
	
	logger.logTMV("tmv", expected);
	logger.logTMV("set",  component.initSet);
	logger.log(sbuilder() << "is close: " << expected.isClose(component.initSet, 0));
}
void comp() {
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parseSetting.addVar("c");
	parseSetting.addVar("d");
  TaylorModelVec tmv = parseTMV("my models {a, b, c,d}");
  vector<HornerForm> hfs = parseHFFromPoly("my hfs {b+a,b,c,d+c}");
  
  MyComponent c1;
  c1.addVar(0);
  c1.addVar(1);
  MyComponent c2;
  c2.addVar(2);
  MyComponent c3;
  c3.addVar(3);
  c3.addDependency(2, &c2);
    
  vector<MyComponent *> comps;
	comps.push_back(&c1);
	comps.push_back(&c2);
	comps.push_back(&c3);
	
	prepareComponents(comps, tmv, hfs, domain);
	
	
  //c1.log();
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
  TaylorModelVec tmv1 = parseTMV("my models {a,b}");
  vector<HornerForm> hfs1 = parseHFFromPoly("my hfs {a+b,b}");
	logger.log(c1.initSet.isClose(tmv1, 1e-20));
	logger.log(c1.odes.size() == hfs1.size());
	for(int i = 0; i < c1.odes.size(); i++) {
  	logger.log(c1.odes[0].isClose(hfs1[0], 1e-20));
  }
  logger.log(c1.varIndexes == parseiVec("my iv <0,1>"));
  logger.log(c1.solveIndexes == parseiVec("my iv <0,1>"));
  logger.log(c1.tpIndexes == parseiVec("my iv <0,1>"));
  logger.log(c1.allTMParams == parseiVec("my iv <0,1>"));
  
  logger.log(c1.compMappers.size() == 0);
  
  c3.log();
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("c");
	parseSetting.addVar("d");
  TaylorModelVec tmv3 = parseTMV("my models {d, c}");
  vector<HornerForm> hfs3 = parseHFFromPoly("my hfs {c + d, 0}");
	logger.log(c3.initSet.isClose(tmv3, 1e-20));
	logger.log(c3.odes.size() == hfs3.size());
	for(int i = 0; i < c3.odes.size(); i++) {
  	logger.log(c3.odes[i].isClose(hfs3[i], 1e-20));
  }
  logger.log(c3.varIndexes == parseiVec("my iv <3>"));
  logger.log(c3.solveIndexes == parseiVec("my iv <0>"));
  logger.log(c3.tpIndexes == parseiVec("my iv <3>"));
  logger.log(c3.allTMParams == parseiVec("my iv <2,3>"));
  logger.log(c3.allTMParams == parseiVec("my iv <2,3>"));
  
  logger.log(c3.compMappers.at(0) == parseiVec("my iv <0,-2>"));
  
  exit(0);
}



void applySw() {
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  
  vector<Interval> step_exp_table;
  vector<Interval> step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 0.1, 2*2);
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	
	OutputWriter writer("dummy", 0, 1);
	
  TaylorModelVec parsed = parseTMV("my models {1+2*a + [-1,2],2+b}");
  vector<HornerForm> hfs = parseHFFromPoly("my hfs {1,2}");
		
  MyComponent c1;
  c1.addVar(0);
  
  MyComponent c2;
  c2.addVar(1);
  
	vector<MyComponent *> comps;
	comps.push_back(&c1);
	comps.push_back(&c2);
		
	prepareComponents(comps, parsed, hfs, domain);
	
	//pipes will be populated in actual integration
	c1.pipes.push_back(c1.initSet);
	c2.pipes.push_back(c2.initSet);
	
  MyComponent all = getSystemComponent(comps, parsed, hfs, domain);
  
  smallComp::applyShrinkWrapping(all, domain, step_end_exp_table, 
      comps, writer);
      
  logger.logTMV("c1", c1.initSet);
  logger.logTMV("c2", c2.initSet);
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec e1 = parseTMV("my models {1+4*a}");
  TaylorModelVec e2 = parseTMV("my models {2+2*a}");
	logger.log(sbuilder() << "is close: " << e1.isClose(c1.initSet, 0));
	logger.log(sbuilder() << "is close: " << e2.isClose(c2.initSet, 0));
}
void sw3() {
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  
  vector<Interval> step_exp_table;
  vector<Interval> step_end_exp_table;
	construct_step_exp_table(step_exp_table, step_end_exp_table, 0.1, 2*2);
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parseSetting.addVar("c");
	parseSetting.addVar("d");
	parseSetting.addVar("e");
	
	OutputWriter writer("dummy", 0, 1);
	
  TaylorModelVec parsed = parseTMV("my models {1+a,2+[-1,2],3+c,4+d,5+e + [-2,2]}");
  vector<HornerForm> hfs = parseHFFromPoly("my hfs {1,2,3,4,5}");
		
  MyComponent c1;
  c1.addVar(0);
  c1.addVar(1);
  c1.addVar(3);
  
  MyComponent c2;
  c2.addVar(2);
  
  MyComponent c3;
  c3.addVar(4);
  c3.addDependency(0, &c1);
  
	
	vector<MyComponent *> comps;
	comps.push_back(&c1);
	comps.push_back(&c2);
	comps.push_back(&c3);
	
	c1.retainEmptyParams = true;
	c2.retainEmptyParams = true;
	c3.retainEmptyParams = true;
	
	prepareComponents(comps, parsed, hfs, domain);
	
	//pipes will be populated in actual integration
	c1.pipes.push_back(c1.initSet);
	c2.pipes.push_back(c2.initSet);
	c3.pipes.push_back(c3.initSet);
  
  MyComponent all = getSystemComponent(comps, parsed, hfs, domain);
  
  /*
  logger.logTMV("c1", c1.initSet);
  logger.logTMV("c2", c2.initSet);
  logger.logTMV("c3", c3.initSet);
  */
  smallComp::applyShrinkWrapping(all, domain, step_end_exp_table, 
      comps, writer);
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	TaylorModelVec e2 = parseTMV("my models{3 + 3*a}");
	parseSetting.addVar("b");
	parseSetting.addVar("c");
	TaylorModelVec e1 = parseTMV("my models{1 + 3*a,2.5 + 4.5*b,4 + 3*c}");
	parseSetting.addVar("d");
	TaylorModelVec e3 = parseTMV("my models{5 + 3*d, 1 + 3*a}");
	
	logger.log(e1.isClose(c1.initSet, 1e-14));
	logger.log(e2.isClose(c2.initSet, 1e-14));
	logger.log(e3.isClose(c3.initSet, 1e-14));
	
}


int main() {
  /*
  self();
  eq();
  difd();
  difs();
  transf();
  transfe();
  poly();
  transf_padding();
  component();
  sw1();
  sw2();
  swSet();
  comp();
  applySw();
  */
  sw3();
}
