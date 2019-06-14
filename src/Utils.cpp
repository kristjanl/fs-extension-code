#include <sys/time.h>
#include <sys/resource.h>

#include "Utils.h"

#include "TaylorModel.h"
#include "MyComponent.h"
#ifdef DInt
  #include "DoubleInterval.h"
#else
  #include "Interval.h"
#endif
#include "OutputWriter.h"
#include "Transformer.h"
#include "PreconditionedTMV.h"


using namespace std;

//in .cpp since then header doesn't have to depend on TaylorModelVec.h
class NamedTMV {
  public:
    string name;
    TaylorModelVec tmv;
    NamedTMV(string name, TaylorModelVec tmv);
};
vector<NamedTMV> pDeserializeNamedFlows(string filename);


MySettings::MySettings() : useFlow(false), discardEmptyParams(false), 
      autoComponents(false), only1Var(false) {
}
MySettings::MySettings(OutputWriter *writer, int order, 
      double step, double time, vector<Interval> estimation, 
      vector<Interval> step_exp_table, 
      vector<Interval> step_end_exp_table, 
      vector<Interval> domain, const Interval *cutoff)
      : writer(writer), maxOrder(order), step(step), time(time), 
      estimation(estimation), step_exp_table(step_exp_table), 
      step_end_exp_table(step_end_exp_table), domain(domain), cutoff(cutoff), 
      discardEmptyParams(false), autoComponents(false), only1Var(false) {
}

void MySettings::log() {
  mreset(old);
  mlog1("setting2 <");
  minc();
  mlog1(sbuilder() << "autoComponents: " << autoComponents);
  mlog1(sbuilder() << "order: " << maxOrder);
  mlog1(sbuilder() << "step: " << step);
  mlog1(sbuilder() << "time: " << time);
  mlog1(sbuilder() << "useFlow: " << useFlow);
  mlog1(sbuilder() << "discardEmptyParams: " << discardEmptyParams);
  mlog("estimation", estimation);
  if(step_exp_table.size() < 2) {
    mlog1("step_exp_table is empty");
  } else {
    mlog1(sbuilder() << "step_exp_table[1]: " << step_exp_table[1].toString());
  }
  if(step_end_exp_table.size() < 2) {
    mlog1("step_end_exp_table is empty");
  } else {
    mlog1(sbuilder() << "step_end_exp_table[1]: " << step_end_exp_table[1].toString());
  }
  mlog("domain", domain);
  mlog1(sbuilder() << "cutoff: " << cutoff->toString(5));
  mlog("varNames", varNames);
  mlog1(sbuilder() << "number of components: " << intComponents.size());
  for(int i = 0; i < intComponents.size(); i++) {
    mlog(sbuilder() << "comp[" << i << "]", intComponents[i]);
  }
  mdec();
  mlog1("setting2 >");
  mrestore(old);
}


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
      if(maxParams < tmvs[i].getIgnoringParamCount())
        maxParams = tmvs[i].getIgnoringParamCount();
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
      fprintf(fp, "%s", s.c_str());
    }
    fprintf(fp, "}\n");
  }
}


void serializeTMV(TaylorModelVec & tmv, string filename) {
  mreset(old);
  mforce1(sbuilder() << "serializing tmv to '" << filename << "'");
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
    mforce1("different vInt sizes");
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
    mforce1(sbuilder() << "compareFlows unequal size (" << 
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


void printNames(string file1, string file2, 
    vector<NamedTMV> & p1, vector<NamedTMV> & p2) {
  cout << "here" << endl;
  
  //dont want to import set
  map<string, int> m1;
  for(int i = 0; i < p1.size(); i++) {
    m1[p1[i].name] = 0;
  }
  map<string, int> m2;
  for(int i = 0; i < p2.size(); i++) {
    m2[p2[i].name] = 0;
  }
  mforce1(sbuilder() << "tmvs in '" << file1 << "'");
  for(map<string, int>::iterator it = m1.begin(); it != m1.end(); it++) {
    mforce1(sbuilder() << "  " << it->first);
  }
  mforce1(sbuilder() << "tmvs in '" << file2 << "'");
  for(map<string, int>::iterator it = m2.begin(); it != m2.end(); it++) {
    mforce1(sbuilder() << "  " << it->first);
  }
  exit(0);
}


void printTMVFiles(string file1, string file2, string name, 
    int index1, int index2) {
  cout << "inspeciting files '" << file1 << "' and '" << file2 << "'" << endl;
  vector<NamedTMV> p1 = pDeserializeNamedFlows(file1);
  vector<NamedTMV> p2 = pDeserializeNamedFlows(file2);
  
  if(name == "") {
    printNames(file1, file2, p1, p2);
  }
  
  
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
  //logger.log("f[0]", first.tms[0]);
  logger.log("f", first);
  logger.log("");
  //logger.log("s[0]", second.tms[0]);
  logger.log("s", second);
  logger.log("");
  
  TaylorModelVec dis = first.distance(second);
  
  mlog("dis", dis);
  //mlog1(sbuilder() << "rem" << dis.tms[0].remainder.toString(50));

  int dim = first.getIgnoringParamCount();
  dim = dim == 0 ? 1 : dim; //check that there are parameters, if not use only t
  
  vector<Interval> domain = getUnitBox(dim);
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
  
  
}


void toMathematica(string file) {
  stringstream ss;
  ss << "mathematica_" << file;
  mforce1(sbuilder() << "writing" << ss.str());
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
    if((lookup[named[i].name]-1)%5 != 0) 
      continue;
    
    //mlog1(sbuilder() << "--" << named[i].name);
    mlog(named[i].name, named[i].tmv);
    
    //mlog(named[i].name, named[i].tmv);
    //mlog1(named[i].tmv.toMathematicaString());

    if(notFirst)
      fprintf(fp, ",");
    notFirst = true;
    
    fprintf(fp, "{");
    fprintf(fp, "\"%s\", ", named[i].name.c_str());
    fprintf(fp, "%s", named[i].tmv.toMathematicaString().c_str());
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
  //mforce1("adding");
  if(active == false)
    return;
  //mforce1(sbuilder() << tmvs.size());
  //mforce1(sbuilder() << "adding '" << name << "'");
  
  TaylorModelVec temp(tmv);
  tmvs.push_back(temp);
  names.push_back(name);
  if(tmvs.size() >= maxSize) {
    serialize();
  }
}

void TMVSerializer::serialize() {
  //cout << "silently stopping serializer\n";
  return;
  mreset(old);
  mlog1("");
  mforce1(sbuilder() << "serializing tmvs to '" << filename << "'");
  //mlog1(sbuilder() << "size: " << tmvs.size());
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
void TMVSerializer::deactivate() {
  active = false;
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

//TODO remove after refactoring
void printComponents(MySettings *settings) {
  //printing out
  mreset(old);
  string s;
  for(int i = 0; i < settings->intComponents.size(); i++) {
    vector<int> comp = settings->intComponents[i];
    s.append(sbuilder() << "[" << settings->varNames[comp[0]]);
    for(int j = 1; j < comp.size(); j++) {
      s.append(sbuilder() << "," << settings->varNames[comp[j]]);
    }
    s.append(sbuilder() << "]");
  }
  cout << "components: " << s << endl;
  mrestore(old);
}

int findPos(int value, const vector<int> *v) {
  return find(v->begin(), v->end(), value) - v->begin();
}

int isIn(int value, const vector<int> *v) {
  return find(v->begin(), v->end(), value) != v->end();
}


void addMyInfo(vector<string> & info) {
  for(map<string,double>::iterator it = timeLookup.begin(); 
      it != timeLookup.end(); it++) {
    //cout << it->first << "\t" << it->second << endl;
    info.push_back(sbuilder() << it->first << ": " << it->second);
  }
  return;


  taddToInfo("remap 1", tr_remap1, info);
  taddToInfo("eval t", tr_eval, info);
  taddToInfo("precond time", tr_precond, info);
  taddToInfo("remap 2", tr_remap2, info);
  taddToInfo("int time", sc_integrate, info);
  
  
  taddToInfo("picard poly", sc_int_poly, info);
  taddToInfo("picard remainder", sc_int_rem, info);
  taddToInfo("picard decreasing", sc_int_find_dec, info);
  taddToInfo("picard refining", sc_int_refine, info);
  
  mlog1(sbuilder() << info.size());
  
  //cout << timeLookup["tr_eval"] << endl;
  taddToInfo("pstart", pre_start, info);
  taddToInfo("pm", pre_matrix, info);
  taddToInfo("pltr", pre_ltr, info);
  taddToInfo("prrange", pre_r_range, info);
  taddToInfo("pinsert", pre_insert, info);
  taddToInfo("pscaling", pre_scaling_m, info);
  taddToInfo("plintra", pre_lin_trans, info);
  taddToInfo("pleft", pre_left, info);
  taddToInfo("pend", pre_end, info);
  
  taddToInfo("refstart", ref_start, info);
  taddToInfo("reffirst", ref_first_picard, info);
  taddToInfo("refrem", ref_rem, info);
  taddToInfo("refsub", ref_subset, info);
}

void addFlowInfo(vector<string> & info) {
  taddToInfo("eval t", fl_eval, info);
  taddToInfo("precond time", fl_precond, info);
  taddToInfo("int time", fl_integrate, info);
  
  taddToInfo("picard poly", fl_int_poly, info);
  taddToInfo("picard remainder", fl_int_rem, info);
  taddToInfo("picard decreasing", sc_int_rem_setup, info);
  taddToInfo("picard refining", fl_int_refine, info);
  
  taddToInfo("reffirst", fl_ref_first_picard, info);
  
  
  
  taddToInfo("pstart", fl_pre_start, info);
  taddToInfo("pm", fl_pre_matrix, info);
  taddToInfo("pltr", fl_pre_ltr, info);
  taddToInfo("prrange", fl_pre_r_range, info);
  taddToInfo("pinsert", fl_pre_insert, info);
  taddToInfo("pscaling", fl_pre_scaling_m, info);
  taddToInfo("plintra", fl_pre_lin_trans, info);
  taddToInfo("pleft", fl_pre_left, info);
  taddToInfo("pend", fl_pre_end, info);
  
}

TaylorModelVec getUnitTmv(int varCount) {
  vector<TaylorModel> tms;
  tms.reserve(varCount);
  for(int i = 0; i < varCount; i++) {
    vector<Interval> temp;
    temp.push_back(Interval(0));
    for(int j = 0; j < varCount; j++) {
      if(i == j) {
        temp.push_back(Interval(1));
        continue;
      }
      temp.push_back(Interval(0));
    }
	  tms.push_back(TaylorModel(Polynomial(temp), Interval(0)));
  }
  TaylorModelVec ret(tms);
  return ret;
}
TaylorModelVec getNVarMParam(int varCount, int paramCount) {
  if(varCount == 0)
    return TaylorModelVec();
  int shift = paramCount - varCount;
  vector<TaylorModel> tms;
  Matrix coefs(varCount, paramCount + 1);
  for(int i = 0; i < varCount; i++) {
    for(int j = 0; j < paramCount; j++) {
      if(i == j) {
        coefs.set(1, i, j + 1 + shift);
      }
    }
  }
  TaylorModelVec ret(coefs);
  return ret;
}
TaylorModelVec getNVarMParam(int varCount, vector<int> params) {
  if(varCount == 0)
    return TaylorModelVec();
  vector<TaylorModel> tms;
  for(int i = 0; i < varCount; i++) {
    vector<Interval> temp;
    temp.push_back(Interval(0));
    for(int j = 0; j < varCount; j++) {
      if(isIn(j, &params) == false) {
        continue;
      }
      if(i == j) {
        temp.push_back(Interval(1));
        continue;
      }
      temp.push_back(Interval(0));
      //temp.push_back(Interval(1));
    }
    tms.push_back(TaylorModel(Polynomial(temp), Interval(0)));
  }
  TaylorModelVec ret(tms);
  return ret;
}


void createFullyCompositionalOutput(vector<MyComponent *> comps, 
    MyComponent & all, Transformer *transformer, MySettings *settings) {
  mreset(old);
  mdisable();
  mlog1("fully comp");
  
  mlog("allVars", all.allVars);
  mlog("varIndexes", all.varIndexes);
  mlog1(sbuilder() << all.dependencies.size());
  
  int pipes = all.dependencies[0]->pComp->pipePairs.size();
  all.pipePairs.reserve(pipes);
  TaylorModel left, right;
  
  for(int step = 0; step < pipes; step++) {
    mlog1(sbuilder() << "step: " << step);
    
    TaylorModelVec left, right;
    for(int var = 0; var < all.allVars.size(); var++) {
      mlog1(sbuilder() << "var: " << var);
      
      
      vector<int> lMapper;
      vector<int> rMapper;
      MyComponent *curComp = &all;
      
      
      while(true) {
        mlog1("in while");
        int pos = findPos(var, &curComp->varIndexes);
        mlog1(sbuilder() << "pos: " << pos);
        
        if(pos < curComp->varIndexes.size()) {
          mlog1("can return straight");
          
          TaylorModel singleLeft = curComp->pipePairs[step]->left.tms[pos];
          TaylorModel singleRight = curComp->pipePairs[step]->right.tms[pos];
          
          mlog("l", singleLeft);
          mlog("r", singleRight);
          
          left.tms.push_back(singleLeft.transform(lMapper));
          right.tms.push_back(singleRight.transform(rMapper));
          break;
        }
        
        
        bool foundDep = false;
        for(int j = 0; j < curComp->dependencies.size(); j++) {
          CompDependency *dep = curComp->dependencies[j];
          mlog("depAll", dep->pComp->allVars);
          
          //is variable in dependency implicit variables
          bool in = isIn(var, &dep->pComp->allVars);
          
          //if variable is not in dependency skip dependency
          if(in == false) {
            continue;
          }
          foundDep = true;
          curComp = dep->pComp;
          
          lMapper = concateMapper(dep->leftMapper, lMapper);
          rMapper = concateMapper(dep->rightMapper, rMapper);
          mlog("lm", lMapper);
          mlog("rm", rMapper);
        }
        if(foundDep == false) {
          stringstream ss;
          ss << "should never get to this point, ";
          ss << "either variable is in componenent or one of the dependencies";
          throw std::runtime_error(ss.str());
        }
      }
    }//end of var
    
    mlog("left", left);
    mlog("right", right);
    all.pipePairs.push_back(new PrecondModel(left, right));
  }//end of step
  mrestore(old);
}


void createOutput(vector<MyComponent *> comps, MyComponent & all, 
      Transformer *transformer, MySettings *settings) {
  tstart(sc_post_composing);
  mreset(old);
  //mdisable();
  if(transformer->isPreconditioned == false) {
    //last one manually since transforming is before integration
    //cout << "size: " << comps[0]->pipes.size() << endl;
    //comps[0]->pipes.push_back(comps[0]->timeStepPipe);
    return;
  }
  mlog1("making system flowpipes");  
  //if fully compositional  
  tstart(tr_remap3);
  
  if(transformer->getType() == TR_SINGLE_COMP) {
    //mforce1("SINGLE COMP");
    //transformer preconditions single component at a time
    //need to add last flowpipe for all components
    for(int i = 0; i < comps.size(); i++) {
      comps[i]->pipePairs.push_back(
          new PrecondModel(comps[i]->timeStepPipe, comps[i]->unpairedRight));
    }
    createFullyCompositionalOutput(comps, all, transformer, settings);
  } else if(transformer->getType() == TR_ALL_COMP) {
    //mforce1("ALL COMP");
    //tranformer maps everything to system, then precondtions
    //need to remap last integration result, add last flowpipe for system component
    all.remapTimeStepPipe();
    all.pipePairs.push_back(new PrecondModel(all.timeStepPipe, all.unpairedRight));
    //mforce1(all.getVarName(settings));
    //mforce("all.tsp", all.timeStepPipe);
    //mforce("all.tsp", all.unpairedRight);
    /*
    mforce1(sbuilder() << all.dependencies[0]->linkVar);
    mforce1(sbuilder() << all.dependencies[1]->linkVar);
    mforce1(sbuilder() << all.dependencies[2]->linkVar);
    mforce1(sbuilder() << all.dependencies[3]->linkVar);
    mforce1(sbuilder() << all.dependencies[4]->linkVar);
    mforce1(sbuilder() << all.dependencies[5]->linkVar);
    mforce1("end");
    all.log();
    */
    //mforce3(old3, "all.right3", all.unpairedRight);
  } else {
    //old code might not be compatible, look into it when problems arise
    throw std::runtime_error("not updated");
  }
  tend(tr_remap3);
  mlog1("composing flowpipes");
  all.output.reserve(all.pipePairs.size());
  for(int i = 0; i < all.pipePairs.size(); i++) {
    cerr << ".";
    //pSerializer->add(all.pipePairs[i]->left, "comp_left");
    //pSerializer->add(all.pipePairs[i]->right, "comp_right");
    TaylorModelVec composed = all.pipePairs[i]->composed(settings, &all);
    
    //pSerializer->add(composed, "composed");
    all.output.push_back(composed);
    //pSerializer->add(composed, "composed");
  }
  cout << all.output.size();
  cout << endl;
  tend(sc_post_composing);
  mrestore(old);
}



long getTime() {
  struct rusage usage1;
  int res1 = getrusage(RUSAGE_SELF, &usage1);
  return  (long)(usage1.ru_utime.tv_sec + usage1.ru_stime.tv_sec) * 1000000
      + usage1.ru_utime.tv_usec + usage1.ru_stime.tv_usec;
}

string getDiff(long t1, long t2, bool div=true) {
  stringstream ss;
  if(div)
    ss << (t2 - t1) / 1e6;
  else
    ss << (t2 - t1);
  return ss.str();
}
