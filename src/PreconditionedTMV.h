#ifndef PRECONDITIONEDTMV_H_
#define PRECONDITIONEDTMV_H_

#include "TaylorModel.h"

class MySettings;
class MyComponent;

class PrecondModel {
  public:
    PrecondModel(TaylorModelVec left, TaylorModelVec right);
    
    TaylorModelVec left;
    TaylorModelVec right;
    
    TaylorModelVec composed(MySettings *settings, MyComponent *component);
};


#endif /* TAYLORMODEL_H_ */