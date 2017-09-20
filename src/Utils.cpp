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
  logger.log(sbuilder() << "useRemainder: " << useRemainder);
  logger.log(sbuilder() << "useSteps: " << useSteps);
  logger.log(sbuilder() << "steps: " << steps);
  logger.log(sbuilder() << "cycleSteps: " << cycleSteps);
  
}

bool ShrinkWrappingCondition::checkApplicability(vector<MyComponent *> comps, 
      const vector<Interval> & estimation) {
  //log();
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
  
  vector<Interval> rightRange;
	right.polyRange(rightRange, settings->domain);
	
	left.insert_ctrunc(ret, right, rightRange, settings->domain, settings->order, 0);
	
	return ret;
}


namespace utilsprivate {
  vector<string> createParams(const TaylorModelVec & tmv) {
    //logger.log("creating");
    
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
  logger.reset();
  logger.force(sbuilder() << "serializing tmv to '" << filename << "'");
  logger.logTMV("tmv", tmv);
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
  logger.log("serializing");
  logger.log(sbuilder() << "writing flows to " << filename);
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
  logger.log(sbuilder() << "reading flows from " << filename);
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  logger.log(sbuilder() << "parsed pipes size: " << parseResult.pipes.size());
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    //logger.logTMV("p", *pipeIt);
  }
  return v;
}


vector<TaylorModelVec *> pDeserializeFlows(string filename) {
  logger.log(sbuilder() << "reading flows from " << filename);
  parseResult.pipes.clear();
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  logger.log(sbuilder() << "parsed *pipes size: " << parseResult.pipes.size());
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
    logger.force("different vInt sizes");
    exit(0);
    return 0;
  }
  for(int i = 0; i < f.size(); i++) {
    Interval fSup, sSup, fInf, sInf;
    
    f[i].sup(fSup);
    s[i].sup(sSup);
    
    f[i].inf(fInf);
    s[i].inf(sInf);
    
    //logger.log(fSup.toString());
    //logger.log(sSup.toString());
    //logger.log(fInf.toString());
    //logger.log(sInf.toString());
    
    Interval supDif = fSup - sSup;
    Interval infDif = fInf - sInf;
    
    Interval infMag;
    infDif.mag(infMag);
    Interval supMag;
    supDif.mag(supMag);
    
    Interval totalDif = infDif + supDif;
    logger.log(sbuilder() << "totalDif: " << totalDif.toString());
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
  logger.log(sbuilder() << "comparing, sizes: " 
      << first.size() << ", " << second.size());
  
  if(first.size() != second.size()) {
    logger.force(sbuilder() << "compareFlows unequal size (" << 
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
    logger.logTMV("f", *f);
    logger.logTMV("s", *s);
    
    TaylorModelVec dis = f->distance(*s);
    vector<Interval> totalDistance;
    dis.intEval(totalDistance, domain);
    
    
    logger.logTMV("dis", dis);
    //logger.logVI("totalDistance", totalDistance);
    
    double sumOfMag = sumVImag(totalDistance);
    sumOfMags += sumOfMag;
    logger.log(sbuilder() << "sumOfMag: " << sumOfMag);
    exit(0);
    
    
    /*
    f->sub(dif, *s);
    //logger.logTMV("dif", dif);
    logger.logTMV("dif", dif);
    vector<Interval> difBound;
    dif.polyRange(difBound, domain);
    //logger.log(sbuilder() << "difBound[0]: " << difBound[0].toString());
    logger.logVI("difBound", difBound);
    
    for(int i = 0; i < difBound.size(); i++)
      sumOfMags += difBound[i].mag();
    
    vector<Interval> fRems = f->getRemainders();
    vector<Interval> sRems = s->getRemainders();
    
    //logger.logVI("fRems", fRems);
    //logger.logVI("sRems", sRems);
    sumOfMags += compareIntervalVecs(fRems, sRems);
    */
  }
  logger.log(sbuilder() << "sum of magnitudes: " << sumOfMags);
}


void printTMVFiles(string file1, string file2, int index1, int index2) {
  logger.log("printing");
  vector<TaylorModelVec *> v1 = pDeserializeFlows(file1);
  vector<TaylorModelVec *> v2 = pDeserializeFlows(file2);
  
  
  //why is -1 >= v1.size() true???
  if(index1 != -1 && index1 >= v1.size() || 
      index2 != -1 && index2 >= v2.size()) {
    logger.force(
        sbuilder() << "different sizes: " << v1.size() << ", " << v2.size());
    throw std::invalid_argument("index was bigger than tmvs count");
  }
  
  //set the indexes to last if they were -1
  if(index1 == -1) {
    index1 = v1.size() - 1;
  }
  if(index2 == -1) {  
    index2 = v2.size() - 1;
  }
  logger.log(sbuilder() << "indexes: " << index1 << ", " << index2);
  TaylorModelVec & first = *v1[index1];
  TaylorModelVec & second = *v2[index2];
  
  logger.log();
  logger.logTMV("f", first);
  logger.log();
  logger.logTMV("s", second);
  logger.log();
  
  TaylorModelVec dis = first.distance(second);
  logger.logTMV("dis", dis);
  vector<Interval> totalDistance;
  dis.intEval(totalDistance, getUnitBox(first.getParamCount()));
  logger.logVI("totalDistance", totalDistance);
  
  double sumOfMag = sumVImag(totalDistance);
  logger.log(sbuilder() << "sumOfMag: " << sumOfMag);
}


vector<Interval> getUnitBox(int n) {
  vector<Interval> ret;
  for(int i = 0; i < n; i++) {
    ret.push_back(Interval(-1,1));
  }
  return ret;
}

TMVSerializer::TMVSerializer(string filename) {
  TMVSerializer(filename, INT_MAX);
}
TMVSerializer::TMVSerializer(string filename, int maxSize) : 
    filename(filename), maxSize(maxSize) {
}
void TMVSerializer::add(const TaylorModelVec & tmv) {
  logger.force("adding");
  logger.force(sbuilder() << tmvs.size());
  logger.force(sbuilder() << (tmvs.size() >= maxSize));
  tmvs.push_back(tmv);
  if(tmvs.size() >= maxSize) {
    logger.force("forced serialization");
    serialize();
  }
}

void TMVSerializer::serialize() {
  logger.reset();
  logger.log();
  logger.force(sbuilder() << "serializing tmvs to '" << filename << "'");
  logger.log(sbuilder() << "size: " << tmvs.size());
  FILE *fp = fopen(filename.c_str(), "w");
  
  
  vector<string> params = utilsprivate::createParams(tmvs[0]);
  utilsprivate::writeParamsToFile(params, fp);
  
  fprintf(fp, "flowpipes{\n");
  for(int i = 0; i < tmvs.size(); i++) {
    logger.logTMV(sbuilder() << "tmv[" << i << "]", tmvs[i]);
    if(i != 0) {
      fprintf(fp, ",\n");
    }
    tmvs[i].serialize(fp,  params);
  }
  
  fprintf(fp, "}\n");
  
  fclose(fp);
  exit(0);
}

