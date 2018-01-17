#ifndef TRANSFORMER_H_
#define TRANSFORMER_H_

#include "include.h"
#include "TaylorModel.h"
#include "MyLogger.h"
#include "MyComponent.h"
#include "Interval.h"
#include "OutputWriter.h"
#include "Continuous.h"
#include "Utils.h"

class MySettings;
class OutputWriter;

using namespace std;

class Transformer {
  public:
    Transformer(bool isPreconditioned, bool isWrapper);
    virtual void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings) = 0;
    void evaluateStepEnd(vector<MyComponent *> & comps, MySettings & settings, 
        bool fail);
    const bool isPreconditioned;
    const bool isWrapper;
    virtual void addInfo(vector<string> & info) = 0;
    virtual void setIntegrationMapper(vector<MyComponent *> comps) = 0;
};


class ShrinkWrapper: public Transformer {
  public:
    ShrinkWrapper(ShrinkWrappingCondition *swChecker);
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void addInfo(vector<string> & info);
    void setIntegrationMapper(vector<MyComponent *> comps);
  private:
    ShrinkWrappingCondition *swChecker;
};

class PreconditionedTransformer: public Transformer {
  public:
    PreconditionedTransformer();
    virtual void getA(Matrix & result, const TaylorModelVec & x0, 
        const int dim) = 0;
    virtual void getAInv(Matrix & result, const Matrix & A) = 0;
    void precond2(TaylorModelVec & badLeft, MySettings & settings, 
        MyComponent & all);
    void precond3(TaylorModelVec & leftStar, MySettings & settings, 
        MyComponent & all);
    void addInfo(vector<string> & info);
    void setIntegrationMapper(vector<MyComponent *> comps);
};

class QRTransformer: public Transformer {
  public:
    QRTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void addInfo(vector<string> & info); //remove after extending preconditioned
    void setIntegrationMapper(vector<MyComponent *> comps); //remove after estending preconditioned
};

class IdentityTransformer: public PreconditionedTransformer {
  public:
    IdentityTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    
    void getA(Matrix & result, const TaylorModelVec & x0, 
        const int dim);
    void getAInv(Matrix & result, const Matrix & A);
    
    void getA(Matrix & result, MyComponent *comp);
    
    void getScaling(Matrix & S, Matrix & SInv, vector<Interval> & rightRange);
    
    TaylorModelVec makeLeftFromA(Matrix & A, MyComponent *comp);
};

class SingleComponentIdentityTransformer: public IdentityTransformer {
  public:
    SingleComponentIdentityTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void preconditionSingleComponent(MyComponent *comp, MySettings & settings);
    void initialPrecondition(MyComponent *comp, MySettings & settings);
};

class NullTransformer: public Transformer {
  public:
    NullTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void addInfo(vector<string> & info);
    void setIntegrationMapper(vector<MyComponent *> comps);
};



#endif /* TRANSFORMER_H_ */
