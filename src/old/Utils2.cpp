#include "Utils2.h"
#include "Interval.h"
#include "OutputWriter.h"
#include "Transformer.h"


MySettings::MySettings() : useFlow(false), discardEmptyParams(false), 
      autoComponents(false) {
}
MySettings::MySettings(OutputWriter *writer, int order, 
      double step, double time, vector<Interval> estimation, 
      vector<Interval> step_exp_table, 
      vector<Interval> step_end_exp_table, 
      vector<Interval> domain, const Interval *cutoff)
      : writer(writer), order(order), step(step), time(time), 
      estimation(estimation), step_exp_table(step_exp_table), 
      step_end_exp_table(step_end_exp_table), domain(domain), cutoff(cutoff), 
      discardEmptyParams(false), autoComponents(false) {
}

void MySettings::log() {
  mreset(old);
  mlog1("setting2 <");
  minc();
  mlog1(sbuilder() << "autoComponents: " << autoComponents);
  mlog1(sbuilder() << "order: " << order);
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

/*
MySettings* MySettings2::toOld() {
  MySettings* ret = new MySettings;
  ret->writer = writer;
  ret->intComponents = intComponents;
  ret->autoComponents = autoComponents;
  ret->order = order;
  ret->step = step;
  ret->time = time;
  ret->estimation = estimation;
  ret->step_exp_table = step_exp_table;
  ret->step_end_exp_table = step_end_exp_table;
  ret->domain = domain;
  ret->cutoff = Interval(*cutoff);
  ret->useFlow = useFlow;
  ret->discardEmptyParams = discardEmptyParams;
  ret->varNames = varNames;
  ret->transformer = transformer;

  return ret;
}*/