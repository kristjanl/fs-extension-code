#include "MyLogger.h"

using namespace std;

/*
template <typename T>
void mylogger2::force(T t) {
  int temp = disabled;
  disabled = 0;
  log(t);
  disabled = temp;
}

template <typename T>
void mylogger2::log(T t) {
	if(disabled > 0)
		return;
	cout << string(ltab, ' ');
	cout << t << endl;
}*/

void mylogger2::force(string s) {
  int temp = disabled;
  disabled = 0;
  log(s);
  disabled = temp;
}

void mylogger2::log(string name, vector<Interval> v) {
	if(disabled > 0)
		return;
	for (unsigned i=0; i<v.size(); i++) {
		logger.log(sbuilder() << name <<"[" << i << "] = " << v.at(i).toString());
	}
  if(v.size() == 0)
    logger.log(sbuilder() << name << " is empty");
}

void mylogger2::log(string name, vector<string> v) {
  if(disabled > 0)
		return;
  if(v.size() == 0) {
    logger.log(sbuilder() << name <<" = [ ]");
    return;
  }
  string s;
  
  vector<string>::iterator it = v.begin();
  s.append(sbuilder() << (*it));
  
  for(it++; it < v.end(); it++) {
    s.append(sbuilder() << ", " << (*it));
  }
  logger.log(sbuilder() << name <<" = [" << s << "]");
}


void mylogger2::log(string name, vector<int> v) {
	if(disabled > 0)
		return;
  if(v.size() == 0) {
    logger.log(sbuilder() << name <<" = [ ]");
    return;
  }
  string s;
  
  vector<int>::iterator it = v.begin();
  s.append(sbuilder() << (*it));
  
  for(it++; it < v.end(); it++) {
    s.append(sbuilder() << ", " << (*it));
  }
  logger.log(sbuilder() << name <<" = [" << s << "]");
}

void mylogger2::log(string name, TaylorModelVec tmv) {
  if(disabled > 0)
		return;
  int dim = tmv.tms.size();
  
  if(tmv.tms.size() == 0)
    logger.log(sbuilder() << name << " is empty");
  
  //logger.log(sbuilder() << "dim: " << dim);
  //logger.log(sbuilder() << "pdim: " << paramDim);
	for (unsigned i=0; i<dim; i++) {
    vector<string> vNames;
    try{
      int paramDim = tmv.tms.at(i).getIgnoringParamCount();
      vNames = getVNames(paramDim);
    } catch(ArgumentException& ignore) { }
    
    //logger.log(paramDim);
    logger.log(sbuilder() << name << "[" << i << "] = " << 
        tmv.tms[i].toString(vNames));
	}
}
void mylogger2::log(string name, TaylorModel tm) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << name << " = " << 
      tm.toString(getVNames(tm.getIgnoringParamCount()))
  );
}

void mylogger2::log(string name, RangeTree *tree) {
	if(disabled > 0)
		return;
  if(tree == NULL) {
    logger.log("NULL");
    return;
  }
  logger.log(sbuilder() << name << " " << tree->ranges.size() << ", " << tree->children.size());
  logger.inc();
   list<Monomial> retained;
  list<Interval>::iterator rIt = tree->ranges.begin();
  list<RangeTree *>::iterator cIt = tree->children.begin();
  for(; rIt != tree->ranges.end(); rIt++) {
    logger.log(rIt->toString());
  }
  for(; cIt != tree->children.end(); cIt++) {
    logger.log(name, *cIt);
  }
  logger.dec();
  /*
	list<Interval> ranges;
	list<RangeTree *> children;*/
}

void mylogger2::log(string name, vector<RangeTree *> trees) {
	if(disabled > 0)
		return;
	for(int i = 0; i < trees.size(); i++) {
  	logger.log(sbuilder() << name << i, trees[i]);
	}
}

void mylogger2::log(string name, vector<HornerForm> hfs) {
  if(disabled > 0)
    return;
  for(int i = 0; i < hfs.size(); i++) {
    logger.log(sbuilder() << name << "[" << i << "]: " << hfs[i].toString());
  }
}

void mylogger2::log(string name, HornerForm hf) {
  if(disabled > 0)
    return;
  logger.log(sbuilder() << name << ": " << hf.toString());
}

void mylogger2::log(Monomial m) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << m.toString(getVNames(10)));
}
void mylogger2::log(string name, Matrix m) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << "matrix " << name);
  
  for(int i = 0; i < m.rows(); i++) {
    sbuilder s;
    for(int j = 0; j < m.cols(); j++) {
      s << m.get(i,j) << "\t";
    }
    logger.log(s);
  }
  
}

void mylogger2::log(string name, set<int> myset) {
  if(disabled > 0)
		return;
  if(myset.size() == 0) {
    logger.log(sbuilder() << name <<" = { }");
    return;
  }
  
  string s;
  set<int>::iterator it = myset.begin();
  s.append(sbuilder() << (*it));
  
  for(it++; it != myset.end(); it++) {
    s.append(sbuilder() << ", " << (*it));
  }
    
  logger.log(sbuilder() << name <<" = {" << s << "}");
}

void mylogger2::log(string s) {
	if(disabled > 0)
		return;
	cout << string(ltab, ' ');
	cout << s << endl;
}

void mylogger2::log(int i) {
	log(sbuilder() << i);
}

void mylogger2::log(Polynomial *p) {
	if(disabled)
		return;
	
	vector<string> vars; //TODO
	vars.push_back("t");
	vars.push_back("x1");
	vars.push_back("x2");
	vars.push_back("x3");
	string s;
	p->toString(s, vars);
	
	cout << string(ltab, ' ');
	cout << s << endl;
}
void mylogger2::inc() {
	ltab += 2;
}
void mylogger2::dec() {
	ltab -= 2;
}
void mylogger2::disable() {
	disabled++;
}
void mylogger2::enable() {
	disabled--;
}
int mylogger2::reset() {
  int old = disabled;
	disabled = 0;
	return old;
}
void mylogger2::restore(int old) {
  disabled = old;
}

mylogger2::mylogger2() {
	this->disabled = 0;
}

mylogger2::mylogger2(bool disabled):disabled(disabled) {
}

mylogger2 logger;
mylogger2 logger2(false);

vector<string> getVNames(int n) {
  vector<string> vars;
	vars.push_back("t");
	char name[10];
	for(int i=0; i<n; ++i) {
		sprintf(name, "%s%d", "a", i+1);
    vars.push_back(name);
	}
  return vars;
}

void mylogger2::log(string name, Polynomial poly) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << name << " = " << 
      poly.toString(getVNames(poly.getVariableCount()))
  );
}
void mylogger2::log(string name, Monomial m) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << name << " = " << 
      m.toString(getVNames(m.getVariableCount())));
}
