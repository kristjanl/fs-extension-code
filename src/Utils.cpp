#include "Utils.h"


ShrinkWrappingCondition::ShrinkWrappingCondition(int steps): 
      steps(steps) {
  useSteps = true;
  useRemainder = false;
  cycleSteps = 0;
  count = 0;
}
//constructor for using remainder
ShrinkWrappingCondition::ShrinkWrappingCondition() {
  useRemainder = true;
  useSteps = false;
  count = 0;
}

void ShrinkWrappingCondition::log() const {
  mlog1(sbuilder() << "useRemainder: " << useRemainder);
  mlog1(sbuilder() << "useSteps: " << useSteps);
  mlog1(sbuilder() << "steps: " << steps);
  mlog1(sbuilder() << "cycleSteps: " << cycleSteps);
  
}

bool ShrinkWrappingCondition::checkApplicability(vector<MyComponent *> comps, 
      const vector<Interval> & estimation) {
  //log("");
  if(useSteps) {
    cycleSteps++;
    if(cycleSteps == steps)
      cycleSteps = 0;
    if(cycleSteps == 0) {
      count++;
      return true;
    }
    return false;
  }
  if(useRemainder) {
    double maxEstimation = -1; //using constant estimations
    for(int i = 0; i < estimation.size(); i++) {
      if(maxEstimation < estimation.at(i).sup())
        maxEstimation = estimation.at(i).sup();
    }
    for(int i = 0; i < comps.size(); i++) {
      for(int j = 0; j < comps.at(i)->initSet.tms.size(); j++) {
        Interval & interval = comps[i]->initSet.tms[j].remainder;
        if (interval.width() > 10*maxEstimation) {
          count++;
          return true;
        }
      }
    }
    return false;
  }
  throw std::runtime_error("should never get here");
}

int ShrinkWrappingCondition::getCount() const {
  return count;
}




PrecondModel::PrecondModel(TaylorModelVec left, TaylorModelVec right) : 
      left(left), right(right) {
}


TaylorModelVec PrecondModel::composed(MySettings *settings) {
  TaylorModelVec ret;
	//mforce3(tt1, "comp_left", left);
	//mforce3(tt2, "comp_right", right);
  
  vector<Interval> rightRange;
	//right.polyRange(rightRange, settings->domain);
  right.polyRangeNormal(rightRange, settings->step_exp_table);
	
	//left.insert_ctrunc(ret, right, rightRange, settings->domain, settings->order, 0);
	left.insert_ctrunc_normal(ret, right, rightRange, settings->step_exp_table,
	     settings->domain.size(), settings->order, settings->cutoff);
	
	//logger.log("cl", left.tms[0]);
	
	//pSerializer->add(right, "comp_right");
	//pSerializer->add(left, "comp_left");
	//pSerializer->add(ret, "comp_comp");
	return ret;
}


namespace utilsprivate {
  vector<string> createParams(const TaylorModelVec & tmv) {
    //mlog1("creating");
    
    vector<string> params;
    
    // time will always be a parameter
    params.push_back("t");
    
    //add other parameters
    for(int i = 1; i < tmv.getParamCount(); i++) {
      string param = sbuilder() << "a" << i;
      params.push_back(param);
    }
    return params;
  }
  
  vector<string> createParams(const vector<TaylorModelVec> & tmvs) {
    //mlog1("creating");
    
    vector<string> params;
    
    // time will always be a parameter
    params.push_back("t");
    
    int maxParams = 0;
    for(int i = 0; i < tmvs.size(); i++) {
      if(maxParams < tmvs[i].getParamCount())
        maxParams = tmvs[i].getParamCount();
    }
    
    //add other parameters
    for(int i = 1; i < maxParams; i++) {
      string param = sbuilder() << "a" << i;
      params.push_back(param);
    }
    return params;
  }
  
  void writeParamsToFile(const vector<string> & params, FILE *fp) {
    fprintf(fp, "vars{t"); 
    for(int i = 1; i < params.size(); i++) {
      string s = sbuilder() << ", " << params[i];
      fprintf(fp, s.c_str());
    }
    fprintf(fp, "}\n");
  }
}


void serializeTMV(TaylorModelVec & tmv, string filename) {
  mreset(old);
  mforce(sbuilder() << "serializing tmv to '" << filename << "'");
  mlog("tmv", tmv);
  FILE *fp = fopen(filename.c_str(), "w");
  vector<string> params = utilsprivate::createParams(tmv);
  utilsprivate::writeParamsToFile(params, fp);
  
  fprintf(fp, "flowpipes{\n");
  //tmv.pushConstantToRemainder();
  tmv.serialize(fp,  params);
  fprintf(fp, "}\n");
  fclose(fp);
  exit(0);
}

void serializeFlows(MyComponent *comp, string filename) {
  throw runtime_error("shoudln't call this anymore");
  mlog1("serializing");
  mlog1(sbuilder() << "writing flows to " << filename);
  FILE *fp = fopen(filename.c_str(), "w");
  vector<string> params = utilsprivate::createParams(comp->pipes[0]);
  utilsprivate::writeParamsToFile(params, fp);
  
  fprintf(fp, "flowpipes{\n");
  bool first = true;
  int n = 1;
  for(vector<TaylorModelVec>::iterator pipeIt = comp->pipes.begin();
        pipeIt < comp->pipes.end(); pipeIt++) {
    if(first == false)
      fprintf(fp, ",\n");
    first = false;
    
    //remove
    //pipeIt->tms.erase(pipeIt->tms.begin() + 1);
    //pipeIt->tms[0].remainder = ZERO_INTERVAL;
    //pipeIt->tms[0].expansion.filter(parseiVec("my iv <3,0,0>"));
    
    //pipeIt->pushConstantToRemainder();
    
    pipeIt->serialize(fp,  params);
  }
  fprintf(fp, "}\n");
  fclose(fp);
  exit(0);
}

vector<TaylorModelVec> & deserializeFlows(string filename) {
  mlog1(sbuilder() << "reading flows from " << filename);
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  mlog1(sbuilder() << "parsed pipes size: " << parseResult.pipes.size());
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    //mlog("p", *pipeIt);
  }
  return v;
}


vector<TaylorModelVec *> pDeserializeFlows(string filename) {
  mlog1(sbuilder() << "reading flows from " << filename);
  parseResult.pipes.clear();
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  mlog1(sbuilder() << "parsed *pipes size: " << parseResult.pipes.size());
  vector<TaylorModelVec *> ret;
  
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    //pipeIt should be memory leak
    TaylorModelVec *tmv = new TaylorModelVec(*pipeIt);
    ret.push_back(tmv);
  }
  return ret;
}


double compareIntervalVecs(vector<Interval> & f, vector<Interval> & s) {
  double sumOfMags = 0;
  if(f.size() != s.size()) {
    mforce("different vInt sizes");
    exit(0);
    return 0;
  }
  for(int i = 0; i < f.size(); i++) {
    Interval fSup, sSup, fInf, sInf;
    
    f[i].sup(fSup);
    s[i].sup(sSup);
    
    f[i].inf(fInf);
    s[i].inf(sInf);
    
    //mlog1(fSup.toString());
    //mlog1(sSup.toString());
    //mlog1(fInf.toString());
    //mlog1(sInf.toString());
    
    Interval supDif = fSup - sSup;
    Interval infDif = fInf - sInf;
    
    Interval infMag;
    infDif.mag(infMag);
    Interval supMag;
    supDif.mag(supMag);
    
    Interval totalDif = infDif + supDif;
    mlog1(sbuilder() << "totalDif: " << totalDif.toString());
    sumOfMags += totalDif.mag();
  }
  return sumOfMags;
  //exit(0);
}


double sumVImag(vector<Interval> vec) {
  Interval ret;
  for(int i = 0; i < vec.size(); i++) {
    ret += vec[i];
  }
  return ret.mag();
}


void compareFlows(vector<TaylorModelVec> & first, 
    vector<TaylorModelVec> & second) {
  vector<TaylorModelVec>::iterator fi = first.begin();
  vector<TaylorModelVec>::iterator si = second.begin();
  vector<TaylorModelVec *> fp;
  vector<TaylorModelVec *> sp;
  for(; fi < first.end(); fi++, si++) {
    fp.push_back(&*fi);
    sp.push_back(&*si);
  }
  compareFlows(fp, sp);
}

void compareFlows(vector<TaylorModelVec *> & first, 
    vector<TaylorModelVec *> & second) {
  mlog1(sbuilder() << "comparing, sizes: " 
      << first.size() << ", " << second.size());
  
  if(first.size() != second.size()) {
    mforce(sbuilder() << "compareFlows unequal size (" << 
        first.size() << ", " << second.size() << ")");
  }
  
  
  vector<TaylorModelVec *>::iterator fIt = first.begin();
  vector<TaylorModelVec *>::iterator sIt = second.begin();
  
  
  vector<Interval> domain;
  for(int i = 0; i < (*fIt)->tms[0].getParamCount(); i++) {
    domain.push_back(Interval(-1,1));  //don't really care about time step size
  }
  double sumOfMags = 0;
  
  for(; fIt < first.end() && sIt < second.end(); fIt++, sIt++) {
    TaylorModelVec dif;
    TaylorModelVec *f = *fIt;
    TaylorModelVec *s = *sIt;
    mlog("f", *f);
    mlog("s", *s);
    
    TaylorModelVec dis = f->distance(*s);
    vector<Interval> totalDistance;
    dis.intEval(totalDistance, domain);
    
    
    mlog("dis", dis);
    //mlog("totalDistance", totalDistance);
    
    double sumOfMag = sumVImag(totalDistance);
    sumOfMags += sumOfMag;
    mlog1(sbuilder() << "sumOfMag: " << sumOfMag);
    exit(0);
    
  }
  mlog1(sbuilder() << "sum of magnitudes: " << sumOfMags);
}


void printTMVFiles(string file1, string file2, string name, 
    int index1, int index2) {
  cout << "inspeciting files '" << file1 << "' and '" << file2 << "'" << endl;
  vector<NamedTMV> p1 = pDeserializeNamedFlows(file1);
  vector<NamedTMV> p2 = pDeserializeNamedFlows(file2);
  vector<TaylorModelVec> v1;
  vector<TaylorModelVec> v2;  
  for(int i = 0; i < p1.size(); i++) {
    if(p1[i].name == name) {
      v1.push_back(p1[i].tmv);
    }
  }
  for(int i = 0; i < p2.size(); i++) {
    if(p2[i].name == name) {
      v2.push_back(p2[i].tmv);
    }
  }
  mlog1(sbuilder() << "name: " << name <<
      ", index1: " << index1 << ", index2: " << index2);
  mlog1(v1.size());
  mlog1(v2.size());
  
  if(index1 >= v1.size() && index2 >= v2.size()) {
    throw std::invalid_argument(
        sbuilder() << "index was bigger than tmvs count in both files");
  }
  if(index1 >= v1.size()) {
    throw std::invalid_argument(
        sbuilder() << "index was bigger than tmvs count in " << file1);
  }
  if(index2 >= v2.size()) {
    throw std::invalid_argument(
        sbuilder() << "index was bigger than tmvs count in " << file2);
  }
  TaylorModelVec & first = v1[index1];
  TaylorModelVec & second = v2[index2];
  
  logger.log("");
  logger.log("f[0]", first.tms[0]);
  //logger.log("f", first);
  logger.log("");
  logger.log("s[0]", second.tms[0]);
  //logger.log("s", second);
  logger.log("");
  
  TaylorModelVec dis = first.distance(second);
  
  mlog("dis", dis);
  //mlog1(sbuilder() << "rem" << dis.tms[0].remainder.toString(50));
  
  vector<Interval> domain = getUnitBox(first.getParamCount());
  domain[0] = Interval(0, 0.1);
  
  vector<Interval> totalDistance;
  dis.intEval(totalDistance, domain);
  
  
  for(int i = 0; i < totalDistance.size(); i++) {
    Interval mag;
    totalDistance[i].mag(mag);
    string s = sbuilder() << "totalDistance[" << i << "] = " 
        << mag.toString(30);
    cout << s << endl;
    break;
  }
  
  /*
  cout << "----------------------";
  for(int i = 1; i < 31; i++) {
    cout << (i%10);
  }
  cout << endl;
  */
  /*
  mlog("totalDistance", totalDistance);
  
  double sumOfMag = sumVImag(totalDistance);
  mlog1(sbuilder() << "sumOfMag: " << sumOfMag);
  */
  
}


void toMathematica(string file) {
  stringstream ss;
  ss << "mathematica_" << file;
  mforce(sbuilder() << "writing" << ss.str());
  FILE *fp = fopen(ss.str().c_str(), "w");
  
  fprintf(fp, "{");
  
  vector<NamedTMV> named = pDeserializeNamedFlows(file);
  
  vector<string> allowed;
  allowed.push_back("leftStar");
  allowed.push_back("left_after_precond");
  //allowed.push_back("comp_left");
  //allowed.push_back("comp_right");
  allowed.push_back("right_after_precond");
  //allowed.push_back("composed_after_precond");
  //allowed.push_back("composed_before_precond");
  
  
  std::map<string,int> lookup;
  
  bool notFirst = false;
  for(int i = 0; i < named.size(); i++) {
    //mlog(named[i].name, named[i].tmv);
    
    //mlog1(named[i].name);
    
    if(find(allowed.begin(), allowed.end(), named[i].name) == allowed.end()) {
      continue;
    }
    
    lookup[named[i].name] = lookup[named[i].name] + 1;
    //if((lookup[named[i].name]-1)%4 != 0)
    //  continue;
    
    
    //mlog1(sbuilder() << "--" << named[i].name);
    mlog(named[i].name, named[i].tmv);
    
    //mlog(named[i].name, named[i].tmv);
    //mlog1(named[i].tmv.toMathematicaString());

    if(notFirst)
      fprintf(fp, ",");
    notFirst = true;
    
    fprintf(fp, "{");
    fprintf(fp, "\"%s\", ", named[i].name.c_str());
    fprintf(fp, named[i].tmv.toMathematicaString().c_str());
    fprintf(fp, "}");
  }
  
  fprintf(fp, "}\n");
  fclose(fp);
}

vector<Interval> getUnitBox(int n) {
  vector<Interval> ret;
  for(int i = 0; i < n; i++) {
    if(i == 0) {
      ret.push_back(Interval(INT_MIN, INT_MAX)); //fail with time variable
      continue;
    }
    ret.push_back(Interval(-1,1));
  }
  return ret;
}

TMVSerializer::TMVSerializer(string filename) : 
    filename(filename), maxSize(INT_MAX), active(true) {
}
TMVSerializer::TMVSerializer(string filename, int maxSize) : 
    filename(filename), maxSize(maxSize), active(true) {
}
TMVSerializer::TMVSerializer(string filename, int maxSize, bool active) : 
    filename(filename), maxSize(maxSize), active(active) {
}
void TMVSerializer::add(const TaylorModelVec & tmv) {
  add(tmv, "NULL");
}
void TMVSerializer::add(const TaylorModelVec & tmv, string name) {
  //mforce("adding");
  if(active == false)
    return;
  //mforce(sbuilder() << tmvs.size());
  //mforce(sbuilder() << "adding '" << name << "'");
  
  TaylorModelVec temp(tmv);
  tmvs.push_back(temp);
  names.push_back(name);
  if(tmvs.size() >= maxSize) {
    serialize();
  }
}

void TMVSerializer::serialize() {
  mreset(old);
  mlog1("");
  mforce(sbuilder() << "serializing tmvs to '" << filename << "'");
  mlog1(sbuilder() << "size: " << tmvs.size());
  FILE *fp = fopen(filename.c_str(), "w");
  
  
  vector<string> params = utilsprivate::createParams(tmvs);
  utilsprivate::writeParamsToFile(params, fp);
  fprintf(fp, "flowpipes{\n");
  for(int i = 0; i < tmvs.size(); i++) {
    //mlog1(sbuilder() << "i: " << i << " out of " << tmvs.size());
    //mlog(sbuilder() << names[i] << "[" << i << "]", tmvs[i]);
    if(i != 0) {
      fprintf(fp, ",\n");
    }
    fprintf(fp, "%s:", names[i].c_str());
    tmvs[i].serialize(fp,  params);
  }
  fprintf(fp, "}\n");
  
  fclose(fp);
  exit(0);
}

void TMVSerializer::activate() {
  active = true;
}

vector<NamedTMV> pDeserializeNamedFlows(string filename) {
  mlog1(sbuilder() << "reading flows from " << filename);
  parseResult.pipes.clear();
  parseResult.names.clear();
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  mlog1(sbuilder() << "parsed *pipes size: " << parseResult.pipes.size());
  vector<TaylorModelVec *> ret;
  vector<NamedTMV> ret2;
  
  mlog1(parseResult.names.size());
  
  for(int i = 0; i < parseResult.pipes.size(); i++) {
    //mlog(parseResult.names[i], parseResult.pipes[i]);
    //pipeIt should be memory leak
    //TaylorModelVec *tmv = new TaylorModelVec(parseResult.pipes[i]);
    ret2.push_back(NamedTMV(parseResult.names[i], parseResult.pipes[i]));
  }
  return ret2;
}

NamedTMV::NamedTMV(string name, TaylorModelVec tmv) : name(name), tmv(tmv) { }

TMVSerializer *pSerializer;

map<string, double> timeLookup;


void printTimes(string prefix) {
  for(map<string,double>::iterator iter = timeLookup.begin(); 
      iter != timeLookup.end(); ++iter) {
    //only print the Transformer clocks
    if(strncmp(iter->first.c_str(), prefix.c_str(), prefix.length()) == 0) {
      cout << iter->first << ": " << iter->second << endl;
    }
  }
}

void addTimeToInfo(string name, string clockName, vector<string> & infos) {
  infos.push_back(sbuilder() << name << ": " << timeLookup[clockName]);
}



int findPos(int value, vector<int> *v) {
  return find(v->begin(), v->end(), value) - v->begin();
}

int isIn(int value, vector<int> *v) {
  return find(v->begin(), v->end(), value) != v->end();
}

