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

using namespace std;

class Transformer {
  public:
    virtual void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings) = 0;
    void evaluateStepEnd(vector<MyComponent *> & comps, MySettings & settings);
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
class NullTransformer: public Transformer {
  public:
    NullTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
};

//TODO ugly hack, ask about it
class Picker: public Transformer {
  public:
    Picker(int type, ShrinkWrapper & sw, QRTransformer & qr, 
      NullTransformer & null);
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    int type;
    ShrinkWrapper & sw;
    QRTransformer & qr;
    NullTransformer & null;
};
const int PICK_SW = 0;
const int PICK_QR = 1;
const int PICK_NULL = 2;


#endif /* TRANSFORMER_H_ */
