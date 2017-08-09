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


void serializeFlows(MyComponent *comp, string filename) {
  logger.log("serializing");
  logger.log(sbuilder() << "writing flows to " << filename);
  FILE *fp = fopen(filename.c_str(), "w");
  vector<string> params;
  
  // time will always be a parameter
  params.push_back("t");
  fprintf(fp, "vars{t"); 
  
  //add other parameters
  for(int i = 0; i < comp->allTMParams.size(); i++) {
    string param = sbuilder() << "a" << (i+1);
    string s = sbuilder() << ", " << param;
    
    params.push_back(param);
    fprintf(fp, s.c_str());
  }
  fprintf(fp, "}\n");
  
  fprintf(fp, "flowpipes{\n");
  bool first = true;
  int n = 1;
  for(vector<TaylorModelVec>::iterator pipeIt = comp->pipes.begin();
        pipeIt < comp->pipes.end(); pipeIt++) {
    if(first == false)
      fprintf(fp, ",\n");
    first = false;
    
    //logger.logTMV("pip", *pipeIt);
    pipeIt->serialize(fp,  params);
    //parseTMV(sbuilder() << "my models{" << n++ << " * a1 + [-0.001,0.001]}").serialize(fpDumping, varNames, params);
    //break;
  }
  fprintf(fp, "}\n");

  fclose(fp);
}

vector<TaylorModelVec> & deserializeFlows(string filename) {
  logger.log(sbuilder() << "reading flows from " << filename);
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  logger.log(parseResult.pipes.size());
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    //logger.logTMV("p", *pipeIt);
  }
  return v;
}


vector<TaylorModelVec *> pDeserializeFlows(string filename) {
  logger.log(sbuilder() << "reading flows from " << filename);
  parseFile(filename);
  vector<TaylorModelVec> & v = parseResult.pipes;
  logger.log(parseResult.pipes.size());
  vector<TaylorModelVec *> ret;
  
  for(vector<TaylorModelVec>::iterator pipeIt = v.begin();
      pipeIt < v.end(); pipeIt++) {
    //pipeIt should be memory leak
    TaylorModelVec *tmv = new TaylorModelVec(*pipeIt);
    ret.push_back(tmv);
  }
  return ret;
}

void compareIntervalVecs(vector<Interval> & f, vector<Interval> & s) {
  if(f.size() != s.size()) {
    logger.force("different vInt sizes");
    return;
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
  }
  //exit(0);
}


void compareFlows(vector<TaylorModelVec> & first, 
    vector<TaylorModelVec> & second) {
  logger.log("here");
  if(first.size() != second.size()) {
    logger.force("compareFlows unequal size (1)");
    exit(0);  
  }
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
  logger.log("comparing");
  vector<TaylorModelVec *>::iterator fIt = first.begin();
  vector<TaylorModelVec *>::iterator sIt = second.begin();
  
  logger.log(first.size());
  logger.log(second.size());
  
  vector<Interval> domain;
  for(int i = 0; i < (*fIt)->tms[0].getParamCount(); i++) {
    domain.push_back(Interval(-1,1));  //don't really care about time step size
  }
  
  for(; fIt < first.end() && sIt < second.end(); fIt++, sIt++) {
    TaylorModelVec dif;
    TaylorModelVec *f = *fIt;
    TaylorModelVec *s = *sIt;
    //logger.logTMV("f", *f);
    //logger.logTMV("s", *s);
    f->sub(dif, *s);
    //logger.logTMV("dif", dif);
    vector<Interval> difBound;
    dif.polyRange(difBound, domain);
    //logger.log(sbuilder() << "difBound[0]: " << difBound[0].toString());
    logger.logVI("difBound", difBound);
    
    vector<Interval> fRems = f->getRemainders();
    vector<Interval> sRems = s->getRemainders();
    
    //logger.logVI("fRems", fRems);
    //logger.logVI("sRems", sRems);
    //compareIntervalVecs(fRems, sRems);
    exit(0);
  }
  logger.log(first.size());
  logger.log(second.size());
}

