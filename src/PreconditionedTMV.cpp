#include "PreconditionedTMV.h"
#include "Utils.h"
#include "MyComponent.h"

PrecondModel::PrecondModel(TaylorModelVec left, TaylorModelVec right) : 
      left(left), right(right) {
}

TaylorModelVec PrecondModel::composed(MySettings *settings, MyComponent *component) {
  TaylorModelVec ret;
	//mforce("comp_left", left);
	//mforce("comp_right", right);
  
  vector<Interval> rightRange;
	//right.polyRange(rightRange, settings->domain);
  right.polyRangeNormal(rightRange, settings->step_exp_table);
	
	//left.insert_ctrunc(ret, right, rightRange, settings->domain, settings->order, 0);
	//left.insert_ctrunc_normal(ret, right, rightRange, settings->step_exp_table,
	//     settings->domain.size(), settings->order, *settings->cutoff);
	left.insert_ctrunc_normal(ret, right, rightRange, settings, component,
	     settings->domain.size());

	//logger.log("cl", left.tms[0]);
	
	//pSerializer->add(right, "comp_right");
	//pSerializer->add(left, "comp_left");
	//pSerializer->add(ret, "comp_comp");
	//mforce("ret", ret);
	return ret;
}