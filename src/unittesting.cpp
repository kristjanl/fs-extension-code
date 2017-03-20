#define BOOST_TEST_MODULE stringtest
#include <boost/test/included/unit_test.hpp>
#include "unittesting.h"

BOOST_AUTO_TEST_SUITE (equality_test)

BOOST_AUTO_TEST_CASE (monomial_self_equality) {
  vector<int> d;
  d.push_back(2);
	Interval I(1);
	Monomial m(I, d);
  BOOST_REQUIRE_EQUAL (m.equals(m), 0);
}

BOOST_AUTO_TEST_CASE (monomial_equality) {
  vector<int> d1;
  d1.push_back(1);
	Interval I1(1);
	Monomial m1(I1, d1);
  
  vector<int> d2;
  d2.push_back(1);
	Interval I2(1);
	Monomial m2(I2, d2);
  BOOST_REQUIRE_EQUAL (m1.equals(m2), 0);
}

BOOST_AUTO_TEST_CASE (monomial_inequality_sizes) {
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
  
  BOOST_CHECK_PREDICATE( std::not_equal_to<int>(), (r)(0) ); 
}

BOOST_AUTO_TEST_CASE (monomial_inequality_degrees) {
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
  
  BOOST_CHECK_PREDICATE( std::not_equal_to<int>(), (r)(0) ); 
}


BOOST_AUTO_TEST_CASE (monomial_transform) {
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
  BOOST_REQUIRE_EQUAL (r, 0);
}


BOOST_AUTO_TEST_CASE (monomial_padding_transform) {
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
  BOOST_REQUIRE_EQUAL (r, 0);
}

BOOST_AUTO_TEST_CASE (monomial_transform_exception) {
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
  map.push_back(0);
  
  try {
    Monomial m3 = m1.transform(map);
    BOOST_REQUIRE_EQUAL ("no exception", "");
  } catch (const std::invalid_argument& e) {
  }
}

BOOST_AUTO_TEST_CASE (subcomponent_transformation) {
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
  
  BOOST_REQUIRE_EQUAL(expected.size(), mappers.size());
  for(int i = 0; i < expected.size(); i++) {
    vector<int> v1 = expected.at(i);
    vector<int> v2 = mappers.at(i);
    BOOST_REQUIRE_EQUAL(v1.size(), v2.size());
    BOOST_REQUIRE_EQUAL(true, equal(v1.begin(), v1.end(), v2.begin()));
  }
}

BOOST_AUTO_TEST_CASE (tmv_equal_1) {
  //equal models
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {3 + 2 * a + [-1,2],3}");
  TaylorModelVec tmv2 = parseTMV("my models {3 + 2 * a + [-1,2],3}");
  BOOST_CHECK( tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_2) {
  //first model not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2],3}");
  TaylorModelVec tmv2 = parseTMV("my models {3 + 2 * a + [-1,2],3}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_3) {
  //second model not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2],3}");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1,2],2}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_4) {
  //different dimensions
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2],3}");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1,2],2,2}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_5) {
  //different different param dims
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2],3}");
	parseSetting.addVar("b");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1,2],2}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_6) {
  //constant part not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2]}");
  TaylorModelVec tmv2 = parseTMV("my models {3 + 2 * a + [-1,2]}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_7) {
  //polynomial part not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2.2 * a + [-1,2]}");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1,2]}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_8) {
  //remainder sup not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2]}");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1,2.2]}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_9) {
  //remainder sup not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2 + 2 * a + [-1,2]}");
  TaylorModelVec tmv2 = parseTMV("my models {2 + 2 * a + [-1.1,2]}");
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-20) );
}
BOOST_AUTO_TEST_CASE (tmv_equal_10) {
  //remainder sup not equal
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
  TaylorModelVec tmv1 = parseTMV("my models {2}");
  TaylorModelVec tmv2 = parseTMV("my models {2.1}");
  BOOST_CHECK(tmv1.isClose(tmv2, 1.1e-1) );
  BOOST_CHECK(false == tmv1.isClose(tmv2, 1e-1) );
}
BOOST_AUTO_TEST_CASE (sw_set) {
  MyComponent component;
  vector<Interval> domain;
	component.pipes.push_back(parseResult.tmv);
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parse("my models {3 + a + [-1,1]}");
	component.initSet = parseResult.tmv;
	
	smallComp::shrinkWrapSet(&component, 2.0, domain);
  TaylorModelVec expected = parseTMV("my models {3 + 2 * a}");
	
  BOOST_CHECK( expected.isClose(component.initSet, 1e-20) );
}

BOOST_AUTO_TEST_CASE (sw_paper) {
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
	
  component.swInput = parseTMV(
      "my models {2 + 4*a + 0.5 * a^2 + [-0.2,0.2],1 + 3*b + a * b + [-0.1,0.1]}");
  double q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  BOOST_CHECK_CLOSE(1.1125, q, 1e-10);
}
BOOST_AUTO_TEST_CASE (sw_linear_coef) {
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
	component.swInput = parseTMV("my models {2 + 0.5*a + [-1,1]}");
  double q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  BOOST_CHECK_CLOSE(3, q, 1e-10);
}
BOOST_AUTO_TEST_CASE (sw_linear_remainder) {
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
	component.swInput = parseTMV("my models {2 + 1*a + [-3,3]}");
  double q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  BOOST_CHECK_CLOSE(4, q, 1e-10);
}
BOOST_AUTO_TEST_CASE (sw_nonlinear) {
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
	component.swInput = parseTMV("my models {2 + a + 0.25 * a^2 + [-1,1]}");
  double q = smallComp::shrinkWrap(component, domain, step_end_exp_table);
  BOOST_CHECK_CLOSE(2 + 1.0/3, q, 1e-10);
}


BOOST_AUTO_TEST_CASE (initial_component_preparation) {
  vector<Interval> domain;
  domain.push_back(Interval(0,0.1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  domain.push_back(Interval(-1,1));
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parseSetting.addVar("c");
  TaylorModelVec tmv = parseTMV("my models {a, b, c}");
  vector<HornerForm> hfs = parseHFFromPoly("my hfs {b+a,b,c}");
  
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
	parseSetting.addVar("c");
  TaylorModelVec tmv2 = parseTMV("my models {a, b, c}");
  vector<HornerForm> hfs2 = parseHFFromPoly("my hfs {a + b,b,c}");
  
  
  MyComponent c1;
  c1.addVar(0);
  c1.addVar(1);
  c1.addVar(2);
    
  vector<MyComponent *> comps;
	comps.push_back(&c1);
	
	prepareComponents(comps, tmv, hfs, domain);
	
  BOOST_CHECK(c1.initSet.isClose(tmv2, 1e-20));
	BOOST_CHECK(c1.odes.size() == hfs2.size());
	for(int i = 0; i < c1.odes.size(); i++) {
  	BOOST_CHECK(c1.odes[i].isClose(hfs2[i], 1e-20));
  }
  BOOST_CHECK(c1.varIndexes == parseiVec("my iv <0,1,2>"));
  BOOST_CHECK(c1.solveIndexes == parseiVec("my iv <0,1,2>"));
  BOOST_CHECK(c1.tpIndexes == parseiVec("my iv <0,1,2>"));
  BOOST_CHECK(c1.allTMParams == parseiVec("my iv <0,1,2>"));
  
  
  BOOST_CHECK(c1.compMappers.size() == 0);
}

BOOST_AUTO_TEST_CASE (multiple_and_dependent_components) {
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
	
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a");
	parseSetting.addVar("b");
  TaylorModelVec tmv1 = parseTMV("my models {a,b}");
  vector<HornerForm> hfs1 = parseHFFromPoly("my hfs {a+b,b}");
	BOOST_CHECK(c1.initSet.isClose(tmv1, 1e-20));
	BOOST_CHECK(c1.odes.size() == hfs1.size());
	for(int i = 0; i < c1.odes.size(); i++) {
  	BOOST_CHECK(c1.odes[i].isClose(hfs1[i], 1e-20));
  }
  BOOST_CHECK(c1.varIndexes == parseiVec("my iv <0,1>"));
  BOOST_CHECK(c1.solveIndexes == parseiVec("my iv <0,1>"));
  BOOST_CHECK(c1.tpIndexes == parseiVec("my iv <0,1>"));
  BOOST_CHECK(c1.allTMParams == parseiVec("my iv <0,1>"));
  
  BOOST_CHECK(c1.compMappers.size() == 0);
  
	parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("c");
	parseSetting.addVar("d");
  TaylorModelVec tmv3 = parseTMV("my models {d, c}");
  vector<HornerForm> hfs3 = parseHFFromPoly("my hfs {c + d, 0}");
	BOOST_CHECK(c3.initSet.isClose(tmv3, 1e-20));
	BOOST_CHECK(c3.odes.size() == hfs3.size());
	for(int i = 0; i < c3.odes.size(); i++) {
  	BOOST_CHECK(c3.odes[i].isClose(hfs3[i], 1e-20));
  }
  BOOST_CHECK(c3.varIndexes == parseiVec("my iv <3>"));
  BOOST_CHECK(c3.solveIndexes == parseiVec("my iv <0>"));
  BOOST_CHECK(c3.tpIndexes == parseiVec("my iv <3>"));
  BOOST_CHECK(c3.allTMParams == parseiVec("my iv <2,3>"));
  BOOST_CHECK(c3.allTMParams == parseiVec("my iv <2,3>"));
  
  BOOST_CHECK(c3.compMappers.at(0) == parseiVec("my iv <0,-2>"));
}

BOOST_AUTO_TEST_SUITE_END( )
