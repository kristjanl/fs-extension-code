#include "Transformer.h"

void Transformer::evaluateStepEnd(vector<MyComponent *> & comps, 
      MySettings & settings) {
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    MyComponent & comp = **it;
    comp.isSolved = false;
    int lastIndex = comp.pipes.size()-1;
    
    //set the next initSet
    comp.pipes.at(lastIndex).evaluate_t(comp.initSet,
        settings.step_end_exp_table);
        
    //logger.logTMV("last", (*it)->pipes.at(lastIndex));
    //logger.listVi("vars", (*it)->varIndexes);
    //logger.logTMV("evaluated", (*it)->initSet);
  }
}


TaylorModelVec shrinkWrapTMV(TaylorModelVec tmv, double factor) {
  TaylorModelVec noConst = TaylorModelVec(tmv);
  noConst.rmConstant(); 
  
  TaylorModelVec constantPart;
  noConst.mul(constantPart, -1);
  constantPart.add_assign(tmv);
  for(int i = 0; i < constantPart.tms.size(); i++) {
    constantPart.tms.at(i).remainder = Interval(0); //nullify the remainder
  }
  
  vector<TaylorModel> qTM;
  for(int i = 0; i < noConst.tms.size(); i++) {
    TaylorModel tm = noConst.tms.at(i);
    TaylorModel mult;
    tm.mul(mult, factor);
    //take only the polynomial part
    qTM.push_back(TaylorModel(mult.expansion, Interval(0)));
  }
  
  TaylorModelVec wrappedSet = TaylorModelVec(qTM);
  wrappedSet.add_assign(constantPart);
  return wrappedSet;
}

//wraps all.swInput
void shrinkWrapSet(MyComponent & all, MyComponent * comp, double factor, 
      vector<Interval> domain) {
  int old = logger.reset();
  logger.disable();
  logger.log("sw set <");
  logger.inc();
  TaylorModelVec orig = comp->initSet;
  //comp->initSet = shrinkWrapTMV(comp->initSet, factor);
  
  logger.logTMV("input", all.swInput);
  TaylorModelVec wrapped = shrinkWrapTMV(all.swInput, factor);
  logger.log("");
  logger.logTMV("wrapped", wrapped);
  
  logger.listVi("compVars", all.compVars);
  
  TaylorModelVec newCompSet;
  
  for(int i = 0; i < comp->compVars.size(); i++) {
    //variable in component
    int varIndex = comp->compVars.at(i);
    //index of that variable in call component
    int indexInAll = find(all.compVars.begin(), all.compVars.end(), varIndex) - 
        all.compVars.begin();
    logger.log(sbuilder() << "index in all: " << indexInAll);
    logger.log(sbuilder() << "init var: " << varIndex);
    //using compVars, cause of the assumption that initial conditions 
    logger.listVi("tm paras", comp->allTMParams);
    logger.logTM("c1", wrapped.tms.at(indexInAll));
    TaylorModel tm = wrapped.tms.at(indexInAll).transform(comp->allTMParams);
    logger.logTM("c2", tm);
    logger.log(sbuilder() << "tm paramCount1: " << wrapped.tms[0].getParamCount());
    logger.log(sbuilder() << "tm paramCount2: " << tm.getParamCount());
    
    newCompSet.tms.push_back(tm);
  }
  
  
  //comp->initSet is the old set, newCompSet is the shrink wrapped one
  if(true) { // TODO add compiler flag
    if(comp->initSet.compare(newCompSet, domain) == false) {
      //throw std::runtime_error("something decreased after shrink wrapping");
    }
  }
  
  
  logger.logTMV("old", comp->initSet);
  comp->initSet = newCompSet;
  logger.logTMV("new", newCompSet);
  
  logger.dec();
  logger.log("sw set >");
  logger.restore(old);
}

void my_handler (const char * reason, const char * file, int line, 
      int gsl_errno) {
  throw IntegrationException(sbuilder() << reason);
}
  
double shrinkWrap(MyComponent & component, vector<Interval> domain, 
     vector<Interval> step_end_exp_table) {
  int old = logger.reset();
  logger.disable();
  logger.log("sw <");
  logger.inc();
  double s = -1, t = -1, d = -1, q;
  
  //TaylorModelVec last = component.pipes.at(component.pipes.size() - 1);
  //TaylorModelVec next;
  //last.evaluate_t(next, step_end_exp_table);
  
  TaylorModelVec next = component.swInput;
  
  //logger.logTMV("next", next);
  TaylorModelVec orig(next);
  logger.logTMV("orig", orig);
  
  int varSize = orig.tms.size();
  int paramCount = next.tms.at(0).getParamCount();
  
  
  //remove the constant part
  next.rmConstant(); 
  
  //logger.logTMV("noconst", next);
  
  
  //logger.log(sbuilder() << d);
  
  vector<int> vars;
  for(int i = 0; i < varSize; i++)
    vars.push_back(i);
  Matrix m(varSize);
  
  //set the elements of matrix m as the linear coefficients
  {
    for(int row = 0; row < varSize; row++) {
      Interval linearPart[varSize];
      next.tms.at(row).getLinearPart(linearPart, vars);
      for(int col = 0; col < varSize; col++) {
        m.set(linearPart[col].sup(), row, col); //TODO taking sup is not correct
      }
    }
  }
  logger.log("matrix M");
  logger.logMatrix(m);
  
  
  //find the inverse of m
  Matrix inv(varSize);
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler (&my_handler);
  try{
    m.inverse(inv);
  }catch(IntegrationException& e) {
    logger.reset();
    logger.logMatrix(m);
    throw e;
  }
  gsl_set_error_handler (old_handler);
  
  logger.log("M^-1");
  logger.logMatrix(inv);
  
  
  //find the nonlinear part of m^-1 * TM
  vector<TaylorModel> linUnitTMs;
  for(int row = 0; row < varSize; row++) {
    //logger.log(sbuilder() << "row:" << row);
    
    //^m-1 * TM part
    TaylorModel result;
    
    for(int col = 0; col < varSize; col++) {
      //logger.log(sbuilder() << "(" << row << "," << col <<")");
      //logger.logd(m.get(row, col));
      TaylorModel temp;
      next.tms.at(col).mul(temp, inv.get(row,col));
      result.add_assign(temp);
    }
    logger.logTM("m^-1 * TM", result);
    
    linUnitTMs.push_back(result);
  }
  logger.logTMV("linear unit TMs", linUnitTMs);
  
  
  
  vector<TaylorModel> nlTMs;
  for(int row = 0; row < varSize; row++) {
    TaylorModel result = linUnitTMs.at(row);
    //substract the linear part (by adding -1*linear part)
    vector<Interval> coefs;
    coefs.push_back(0); //time
    for(int i = 0; i < paramCount - 1; i++) { //-1 because time is counted too
      //logger.log(i);
      if(i == row) // TODO assuming that varIndex == paramIndex
        coefs.push_back(-1);
      else
        coefs.push_back(0);
    }
    logger.logVI("coefs", coefs);
    TaylorModel tm(Polynomial(coefs), Interval(0));
    logger.logTM("-1*lin", tm);
    result.add_assign(tm);
    //logger.logTM("res2", result);
    nlTMs.push_back(result);
  }
  
  
  //nullify remainder (so you get only nonlinear part)
  for(int i = 0; i < nlTMs.size(); i++) {
    nlTMs.at(i).remainder = Interval(0);
  }
  
  TaylorModelVec nlPart(nlTMs);
  
  
  logger.logTMV("nlPart", nlPart);
  
  //logger.logTMV("C^-1 * ", tms);
  //finding d
  for(vector<TaylorModel>::iterator it = linUnitTMs.begin(); 
      it < linUnitTMs.end(); it++) {
    //logger.log(it->remainder.toString());
    //logger.log(sbuilder() << it->remainder.mag());
    if(it->remainder.mag() > d)
      d = it->remainder.mag();
  }
  logger.log(sbuilder() << "d: " << d);
  
  
  
  //finding s
  //logger.logTMV("nl", nlPart);
  //bound the nonlinear part
  vector<Interval> nlRange;
  
  nlPart.polyRange(nlRange, domain);
  
  //logger.logVI("nl", nlRange);
  for(vector<Interval>::iterator it = nlRange.begin();
      it < nlRange.end(); it++) {
    double mag = it->mag();
    s = max(s,mag);
  }
  logger.log(sbuilder() << "s: " << s);
  
  //logger.listVi("params", component.allTMParams);
  
  //find t
  for(int param = 1; param < paramCount; param++) {
    TaylorModelVec pDer;
    nlPart.derivative(pDer, param); //derivative wrt to parameter
    //logger.log(param);
    //logger.logTMV("pDer", pDer);
    
    //bound of the derivative
    vector<Interval> pDerRange;
    pDer.polyRange(pDerRange, domain);
    //logger.logVI("pd", pDerRange);
    for(vector<Interval>::iterator it2 = pDerRange.begin();
        it2 < pDerRange.end(); it2++) {
      double mag = it2->mag();
      //logger.log(sbuilder() << mag);
      t = max(t,mag);
    }
  }
  logger.log(sbuilder() << "t: " << t);
  
  int n = paramCount - 1;
  q = 1 + d*(1/((1-(n-1)*t)*(1-s))); //TODO not safe floating point rounding
  //logger.force(sbuilder() << "f(t): " << (1-(n-1)*t));
  //logger.force(sbuilder() << "f(s): " << (1-s));
  //logger.force(sbuilder() << (1/((1-s))));
  if((1 - s) <= 0) {
    throw IntegrationException(sbuilder() << 
        "Map is not shrinkable (s = " << s << ")");
  }
  if((1 - n*t) <= 0) {
    throw IntegrationException(sbuilder() << 
        "Map is not shrinkable (n = " << n << ", t = " << t << ")");
  }
  
  logger.log(sbuilder() << "1 - s: " << (1-s));
  logger.log(sbuilder() << "1 - nt: " << (1-n*t));
  logger.log(sbuilder() << "q:" << q);
  logger.dec();
  logger.log("sw >");
  logger.restore(old);
  return q;
}


void parametrizeVars(TaylorModelVec & tmv, vector<int> varsToIntroduce, 
     int paramCount) {
  int old = logger.reset();
  logger.disable();
  logger.log("parametrizing <");
  logger.inc();
  //int paramCount = tmv.tms[0].getParamCount();
  for(int i = 0; i < varsToIntroduce.size(); i++) {
    int var = varsToIntroduce[i];
    //logger.log(sbuilder() << "var: " << var);
    TaylorModel & tm = tmv.tms[var];
    
    Interval & rem = tm.remainder;
    
    Interval inf;
    rem.inf(inf);
    Interval sup;
    rem.sup(sup);
    
    Interval shift = inf + sup;
    shift.div_assign(2);
    
    //logger.log(sbuilder() << "shift: " << shift.toString());
    //logger.log(sbuilder() << "inf: " << inf.toString());
    //logger.log(sbuilder() << "sup: " << sup.toString());
    //logger.log(sbuilder() << "tm: " << tm.remainder.toString());
    
    
    tm.expansion.add_assign(Monomial(shift, paramCount));
    Interval coef = (rem - shift).sup();
    
    //logger.log(sbuilder() << "coef: " << coef.toString());
    
    vector<int> degs;
    //logger.log(tm.getParamCount());
    for(int j = 0; j < paramCount; j++) {
      if(j == var + 1)
        degs.push_back(1);
      else
        degs.push_back(0);
    }
    
    tm.expansion.add_assign(Monomial(coef, degs));
    tm.remainder = Interval();
  }
  logger.dec();
  logger.log("parametrizing >");
  logger.restore(old);
}

void introduceParam(MyComponent & comp, vector<MyComponent *> comps, 
      vector<Interval> & domain) {
  int old = logger.reset();
  logger.disable();
  logger.log("introducing <");
  logger.inc();
  TaylorModelVec tmv = comp.swInput;
  
  
  //logger.logTMV("last", last);
  logger.logTMV("original", tmv);
  
  vector<int> varsToIntroduce;
  for(int i = 0; i < comps.size(); i++) {
    vector<int> & compVars = comps[i]->varsToBeIntroduced;
    varsToIntroduce.insert(varsToIntroduce.end(), 
        compVars.begin(), compVars.end());
    logger.listVi("vs: ", comps[i]->varsToBeIntroduced);
    comps[i]->varsToBeIntroduced.clear();
  }
  logger.listVi("introducing", varsToIntroduce);
  
  //no need to add more parameters now
  //int oParamCount = tmv.tms[0].getParamCount();
  
  //TaylorModelVec padded = tmv.addNParams(varsToIntroduce.size());
  TaylorModelVec padded = TaylorModelVec(tmv);
  int paramCount = padded.tms[0].getParamCount();
  
  
  if(varsToIntroduce.size() != 0) {
    parametrizeVars(padded, varsToIntroduce, paramCount);
  }
  
  logger.logTMV("introduced", padded);
  comp.swInput = padded;
  
  
  
  if(true) { // TODO add compiler flag
    if(tmv.compare(comp.swInput, domain) == false) {
      throw std::runtime_error("something decreased after introducing param");
    }
  }
  logger.dec();
  logger.log("introducing >");
  logger.restore(old);
}

double applyShrinkWrapping(MyComponent & all, vector<Interval> domain, 
      vector<Interval> step_end_exp_table, vector<MyComponent *> comps,
      MySettings & settings) {
  int old = logger.reset();
  logger.disable();
  logger.log("applying sw <");
  logger.inc();
  
  clock_t start = clock();
  
  all.remapLastFlowpipe();
  
  TaylorModelVec last = all.pipes.at(all.pipes.size() - 1);
  TaylorModelVec tmv;
  last.evaluate_t(tmv, step_end_exp_table);
  all.swInput = tmv;
  
  logger.logTMV("swInput", all.swInput);
  //introduces a parameter to swInput
  
  
  introduceParam(all, comps, domain);
  logger.log(all.swInput.tms[0].getParamCount());
  /*
  logger.log(all.compMappers.size());
  logger.listVi("0", all.compMappers.at(0));
  logger.listVi("1", all.compMappers.at(1));
  logger.listVi("2", all.compMappers.at(2));
  logger.listVi("3", all.compMappers.at(3));
  logger.listVi("4", all.compMappers.at(4));
  */
  logger.logTMV("swI aft", all.swInput);
  
  
  //logger.logTMV("all", all.pipes.at(all.pipes.size() - 1));
  //vector<TaylorModelVec> p2 = all.dependencies.at(0)->pComp->pipes;
  //logger.logTMV("p2", p2.at(p2.size() - 1));
  
  //calculate the shrink wrapping factor for swInput
  double swQ = shrinkWrap(all, domain, step_end_exp_table);
  
  logger.log(sbuilder() << "swQ: " << swQ);
  
  //use the computed shrink wrapping factor to modify the initial sets
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    //sets comp.initSet from q and all.swInput
    shrinkWrapSet(all, *it, swQ, domain);
  }
  
  clock_t end = clock();
  double singleSW = double(end - start) / CLOCKS_PER_SEC;
  settings.writer->swTime += singleSW;
  //logger.force(sbuilder() << "single: " << singleSW);
  //logger.force(sbuilder() << (end - start));
  
  logger.dec();
  logger.log("applying sw >");
  logger.restore(old);
}




Transformer::Transformer(bool isPreconditioned, bool isWrapper) : 
      isPreconditioned(isPreconditioned), isWrapper(isWrapper) {
}

ShrinkWrapper::ShrinkWrapper(ShrinkWrappingCondition *swChecker) : 
      Transformer(false, true) {
  this->swChecker = swChecker;
}

QRTransformer::QRTransformer() : Transformer(true, false) {
}

NullTransformer::NullTransformer() : Transformer(false, false) {
}
IdentityTransformer::IdentityTransformer() : Transformer(false, false) {
  //isPreconditioned = true;
}

void ShrinkWrapper::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  logger.log("sw transforming");
  evaluateStepEnd(comps, settings);
  if(swChecker->checkApplicability(comps, settings.estimation)) {
    logger.force("wrapping");
    applyShrinkWrapping(all, settings.domain, settings.step_end_exp_table, comps,
        settings);
  }
}

TaylorModelVec getUnitTmv(int varCount) {
  vector<TaylorModel> tms;
  logger.log(varCount);
  for(int i = 0; i < varCount; i++) {
    vector<Interval> temp;
    temp.push_back(Interval(0));
    for(int j = 0; j < varCount; j++) {
      if(i == j) {
        temp.push_back(Interval(1));
        continue;
      }
      temp.push_back(Interval(0));
    }
	  tms.push_back(TaylorModel(Polynomial(temp), Interval(0)));
  }
  TaylorModelVec ret(tms);
  
	logger.logTMV("ret", ret);
  return ret;
}


void precond(TaylorModelVec & tmv, MySettings & settings, 
      MyComponent & all) {
  int old = logger.reset();
  logger.disable();
  logger.log("precond <");
  logger.inc();
  
  logger.log("precond");
  
  logger.logTMV("tmv", tmv);
  
  vector<Interval> previousRange;
  TaylorModelVec previousRight;
  
  //number of taylor model parameters
  int paramCount = tmv.tms[0].getParamCount();
  //number of system variables (in the component)
  int varCount = tmv.tms.size();
  
  //linear TMV where TM[i] = 1*var_i
  TaylorModelVec unitTmv = getUnitTmv(varCount);
  
  
  if(all.pipePairs.size() == 0) {
    previousRight = unitTmv;
  } else {
    previousRight = all.pipePairs[all.pipePairs.size() - 1]->right;
  }
	previousRight.polyRange(previousRange, settings.domain);
  logger.logVI("range", previousRange);
  
  logger.log(sbuilder() << "paramCount: " << paramCount);
  logger.log(sbuilder() << "varCount: " << varCount);
  Matrix A(varCount,varCount), invA(varCount, varCount);
  
  if(true) { //use QR preconditioning
  	preconditionQR(A, tmv, varCount, varCount+1);
    A.transpose(invA);
  } else {
    logger.force("not implementated");
    exit(0);
  }
	
  logger.logMatrix(A);
  logger.logMatrix(invA);
  
  vector<Interval> center;
  tmv.constant(center);
  logger.logVI("center", center);
  
  TaylorModelVec c(tmv), lin, nl, rem;
  tmv.getParts(c, lin, nl, rem);
  
  logger.logTMV("c", c);
  logger.logTMV("lin", lin);
  logger.logTMV("nl", nl);
  logger.logTMV("rem", rem);
  
  nl.linearTrans_assign(invA);
  rem.linearTrans_assign(invA);
  lin.linearTrans_assign(invA);
  
  TaylorModelVec rightPart(lin);
  rightPart.add_assign(nl);
  rightPart.add_assign(rem);
  
  
  //calculate right bound
  
  TaylorModelVec rightComp;
	rightPart.insert_ctrunc(rightComp, previousRight, previousRange, settings.domain, settings.order, 0);
  
  
  vector<Interval> rightBound;
	rightComp.intEvalNormal(rightBound, settings.step_end_exp_table);
  logger.logTMV("rightComp", rightComp);
  logger.logVI("rightBound", rightBound);
  
  //TODO shift this to be around 0?
  
  Matrix S(varCount, varCount);
	Matrix invS(varCount, varCount);
  for(int i = 0; i < varCount; i++) {
    Interval intMag;
    rightBound[i].mag(intMag);
    double dMag = intMag.sup();
    if(intMag.subseteq(Interval()) == false) {
		  S.set(dMag, i, i);
		  invS.set(1/dMag, i, i);
	  }
	  else {
		  S.set(0, i, i);
		  invS.set(0, i, i); //doesn't matter if 0 or 1
		}
  }
  
  logger.logMatrix(S);
  logger.logMatrix(invS);
  
  rightComp.linearTrans_assign(invS);
	rightComp.intEvalNormal(rightBound, settings.step_end_exp_table);
  logger.logTMV("rightComp", rightComp);
  logger.logVI("rightBound", rightBound);
  
	//left part is made of constant part and desired linear part
	TaylorModelVec left(c);
	TaylorModelVec desiredLinearPart;
	unitTmv.linearTrans(desiredLinearPart, A * S);
	left.add_assign(desiredLinearPart);
	
	PrecondModel *pre = new PrecondModel(left, rightComp);
	
	logger.logTMV("left", left);
  logger.logTMV("right", rightComp);
	
  all.pipePairs.push_back(pre);
  
  
  logger.dec();
  logger.log("precond >");
  logger.restore(old);
}


void QRTransformSet(MyComponent & all, MyComponent * comp) {
  int old = logger.reset();
  logger.disable();
  logger.log("QRTransformSet <");
  logger.inc();
      
  PrecondModel *pre = all.pipePairs[all.pipePairs.size() - 1];

  TaylorModelVec newCompSet;
  
  for(int i = 0; i < comp->compVars.size(); i++) {
    //variable in component
    int varIndex = comp->compVars.at(i);
    //index of that variable in call component
    int indexInAll = 
        find(all.compVars.begin(), all.compVars.end(), varIndex) - 
        all.compVars.begin();
    logger.log(sbuilder() << "comp var: " << varIndex);
    logger.log(sbuilder() << "index in all: " << indexInAll);
    //using compVars, cause of the assumption that initial conditions 
    logger.listVi("tm params", comp->allTMParams);
    logger.logTM("c1", pre->left.tms.at(indexInAll));
    TaylorModel tm = pre->left.tms.at(indexInAll).
        transform(comp->allTMParams);
    logger.logTM("c2", tm);
    logger.log(sbuilder() << "tm paramCount1: " << 
        pre->left.tms[0].getParamCount());
    logger.log(sbuilder() << "tm paramCount2: " << tm.getParamCount());
    
    newCompSet.tms.push_back(tm);
  }
  logger.logTMV("init1", comp->initSet);
  comp->initSet = newCompSet;
  logger.logTMV("init2", comp->initSet);
  
  logger.dec();
  logger.log("QRTransformSet >");
  logger.restore(old);
}


void QRTransformer::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  int old = logger.reset();
  logger.disable();
  logger.log("qr transforming <");
  logger.inc();
  evaluateStepEnd(comps, settings);
  
  //remaps all the last pipes to system flowpipes at all component
  all.remapLastFlowpipe();
  
  TaylorModelVec last = all.pipes.at(all.pipes.size() - 1);
  TaylorModelVec tmv;
  last.evaluate_t(tmv, settings.step_end_exp_table);
  
  
  
  precond(tmv, settings, all);
  logger.logTMV("tmv", tmv);
  
  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    QRTransformSet(all, *it);
  }
  
  
  logger.dec();
  logger.log("qr transforming >");
  logger.restore(old);
}

void NullTransformer::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  logger.log("null transforming");
  evaluateStepEnd(comps, settings);
}
void IdentityTransformer::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  logger.log("identity transforming");
  logger.force("no impl");
  exit(0);
  evaluateStepEnd(comps, settings);
}



