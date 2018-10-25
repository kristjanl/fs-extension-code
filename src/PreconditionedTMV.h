#ifndef PRECONDITIONEDTMV_H_
#define PRECONDITIONEDTMV_H_

#include "TaylorModel.h"

class MySettings;

class PrecondModel {
  public:
    PrecondModel(TaylorModelVec left, TaylorModelVec right);
    
    TaylorModelVec left;
    TaylorModelVec right;
    
    TaylorModelVec composed(MySettings *settings);
};


#endif /* TAYLORMODEL_H_ */