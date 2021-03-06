#ifndef COMPSOLVER_H_
#define COMPSOLVER_H_

#include <iostream>
#include <vector>


class MySettings;
//class MySettings2;
class ContinuousSystem;
class HornerForm;
class TaylorModelVec;
class Interval;
class MyComponent;

class IVP {
  public:
    IVP(ContinuousSystem & system);
	  TaylorModelVec *initSet;
	  std::vector<HornerForm> hfOde;
    std::vector<Interval> domain;
};

class Solver {
  private:
    void setUp(MySettings *settings, IVP & ivp);
    void post(MySettings *settings);
  public:
    MyComponent *all; //TODO free
    std::vector<MyComponent *> comps;
    void solveIVP(MySettings *settings, IVP ivp);
};

#endif /* COMPSOLVER_H_ */

