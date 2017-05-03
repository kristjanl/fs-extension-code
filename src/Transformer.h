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
    void evaluateStepEnd(vector<MyComponent *> & comps, MySettings & settings);
    const bool isPreconditioned;
    const bool isWrapper;
};


class ShrinkWrapper: public Transformer {
  public:
    ShrinkWrapper(ShrinkWrappingCondition *swChecker);
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
  private:
    ShrinkWrappingCondition *swChecker;
};

class QRTransformer: public Transformer {
  public:
    QRTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
};

class IdentityTransformer: public Transformer {
  public:
    IdentityTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
};

class NullTransformer: public Transformer {
  public:
    NullTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
};


#endif /* TRANSFORMER_H_ */
