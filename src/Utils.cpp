#include "Utils.h"

void foo() {
  cout << "aa" << endl;
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


