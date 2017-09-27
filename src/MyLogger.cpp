#include "MyLogger.h"

void mylogger2::force(string s) {
  int temp = disabled;
  disabled = 0;
  log(s);
  disabled = temp;
}

void mylogger2::logVI(vector<Interval> v) {
	if(disabled > 0)
		return;
	logger.log("vec interval <");
	logger.inc();
	for (unsigned i=0; i<v.size(); i++) {
		string s;
		v.at(i).toString(s);
		logger.log(sbuilder() << "v[" << i << "] = " << s);
	}
	logger.dec();
	logger.log("vec interval >");
}


void mylogger2::logVI(string name, vector<Interval> v) {
	if(disabled > 0)
		return;
	for (unsigned i=0; i<v.size(); i++) {
		logger.log(sbuilder() << name <<"[" << i << "] = " << v.at(i).toString());
	}
  if(v.size() == 0)
    logger.log(sbuilder() << name << " is empty");
}

void mylogger2::logVS(string name, vector<string> v) {
	if(disabled > 0)
		return;
	for (unsigned i=0; i<v.size(); i++) {
		logger.log(sbuilder() << name <<"[" << i << "] = " << v[i]);
	}
  if(v.size() == 0)
    logger.log(sbuilder() << name << " is empty");
}

void mylogger2::logVi(string name, vector<int> v) {
	if(disabled > 0)
		return;
	for (unsigned i=0; i<v.size(); i++) {
		logger.log(sbuilder() << name <<"[" << i << "] = " << v.at(i));
	}
  if(v.size() == 0)
    logger.log(sbuilder() << name << " is empty");
}



void mylogger2::listVi(string name, vector<int> v) {
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

void mylogger2::logTMV(string name, TaylorModelVec tmv) {
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
      int paramDim = tmv.tms.at(i).getParamCount();
      vNames = getVNames(paramDim);
    } catch(ArgumentException& ignore) { }
    
    //logger.log(paramDim);
    logger.log(sbuilder() << name << "[" << i << "] = " << 
        tmv.tms.at(i).toString(vNames));
	}
}
void mylogger2::logPM(string name, PrecondModel *pm) {
	if(disabled > 0)
		return;
  logger.force("a");
  //logger.logTMV("left", pm->left);
  throw invalid_argument("don't call logPM");
}


void mylogger2::logRT(string name, RangeTree *tree) {
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
    logger.logRT(name, *cIt);
  }
  logger.dec();
  /*
	list<Interval> ranges;
	list<RangeTree *> children;*/
}

void mylogger2::logVRT(string name, vector<RangeTree *> trees) {
	if(disabled > 0)
		return;
	for(int i = 0; i < trees.size(); i++) {
  	logger.logRT(sbuilder() << name << i, trees[i]);
	}
}


void mylogger2::logTMVRem(string name, TaylorModelVec tmv) {
	if(disabled > 0)
		return;
  
  if(tmv.tms.size() == 0)
    logger.log(sbuilder() << name << " is empty");
  
	for (unsigned i=0; i<tmv.tms.size(); i++) {
    logger.log(sbuilder() << name << "[" << i << "] = " << 
        tmv.tms[i].remainder.toString());
	}
}

void mylogger2::logVHF(string name, vector<HornerForm> hfs) {
  if(disabled > 0)
    return;
  for(int i = 0; i < hfs.size(); i++) {
    logger.log(sbuilder() << name << "[" << i << "]: " << hfs[i].toString());
  }
}


void mylogger2::logM(Monomial m) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << m.toString(getVNames(10)));
}
void mylogger2::logMatrix(string name, Matrix m) {
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
void mylogger2::logTM(string name, TaylorModel tm) {
	if(disabled > 0)
		return;
  logger.log(sbuilder() << name << " = " << tm.toString(getVNames(10)));
}

void mylogger2::log() {
  log("");
}

void mylogger2::log(string s) {
	if(disabled > 0)
		return;
	cout << string(ltab, ' ');
	cout << s << endl;
}

void mylogger2::printLevel() {
	cout << string(ltab, ' ');
	cout << disabled << endl;
}
void mylogger2::log(int i) {
	log(sbuilder() << i);
}
void mylogger2::logd(double d) {
	log(sbuilder() << d);
}

void mylogger2::logPoly(Polynomial *p) {
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
