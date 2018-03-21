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


#define TR_UNKNOWN      0
#define TR_ALL_COMP			1
#define TR_SINGLE_COMP	2


using namespace std;

class Transformer {
  private:
    int transformerType;
  public:
    Transformer(bool isPreconditioned, bool isWrapper, int type, string name);
    virtual void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings) = 0;
    void makeNextInitSet(vector<MyComponent *> & comps, MySettings & settings, 
        bool fail);
    const bool isPreconditioned;
    const bool isWrapper;
    const string name;
    virtual void addInfo(vector<string> & info) = 0;
    virtual void setIntegrationMapper(vector<MyComponent *> comps) = 0;
    int getType();
    int count;
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
    PreconditionedTransformer(int type, string name);
    virtual void getMatrices(Matrix & a, Matrix & aInv, 
        const TaylorModelVec & x0) = 0;
    virtual TaylorModelVec getLeftToRight(TaylorModelVec & leftStar, 
        Matrix & invA) = 0;
    void precondition(TaylorModelVec & leftStar, MySettings & settings, 
        MyComponent & all);
    void transformFullSystem(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void addInfo(vector<string> & info);
    void setIntegrationMapper(vector<MyComponent *> comps);
    void getScaling(Matrix & S, Matrix & SInv, vector<Interval> & rightRange);
};

class QRTransformer: public PreconditionedTransformer {
  public:
    QRTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    virtual void getMatrices(Matrix & a, Matrix & aInv, 
        const TaylorModelVec & x0) = 0;
    void addInfo(vector<string> & info); //remove after extending preconditioned
    void setIntegrationMapper(vector<MyComponent *> comps); //remove after extending preconditioned
    TaylorModelVec getLeftToRight(TaylorModelVec & leftStar, Matrix & invA);
};


//using regular (column) based QR
class QRTransformerPlain: public QRTransformer {
  public:
    QRTransformerPlain();
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
};

//using row based QR
class QRTransformer1: public QRTransformer {
  public:
    QRTransformer1();
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
};

//using lin and lin^-1
class QRTransformer2: public QRTransformer {
  public:
    QRTransformer2();
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
};
//using R^T and R
class QRTransformer3: public QRTransformer {
  public:
    QRTransformer3();
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
};

class IdentityTransformer: public PreconditionedTransformer {
  public:
    IdentityTransformer();
    IdentityTransformer(int type, string name);
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
    TaylorModelVec getLeftToRight(TaylorModelVec & leftStar, Matrix & invA);
    
    TaylorModelVec makeLeftFromA(Matrix & A, MyComponent *comp);
};

class SingleComponentIdentityTransformer: public PreconditionedTransformer {
  public:
    SingleComponentIdentityTransformer();
    void transform(MyComponent & all, vector<MyComponent *> & comps, 
        MySettings & settings);
    void preconditionSingleComponent(MyComponent *comp, MySettings & settings);
    void initialPrecondition(MyComponent *comp, MySettings & settings);
    TaylorModelVec makeLeftFromA(Matrix & A, MyComponent *comp);
    
    void getA(Matrix & result, MyComponent *comp);
    void getAInv(Matrix & result, const Matrix & A);
    void getMatrices(Matrix & a, Matrix & aInv, const TaylorModelVec & x0);
    TaylorModelVec getLeftToRight(TaylorModelVec & leftStar, Matrix & invA);
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
