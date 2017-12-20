#include "Transformer.h"

void Transformer::evaluateStepEnd(vector<MyComponent *> & comps, 
      MySettings & settings, bool fail) {
  mreset(old);
  if(true || fail) {
    //signaling need to refactor
    throw invalid_argument("evaluate called with fail"); 
  }
  for(int i = 0; i < comps.size(); i++) {
    mlog1(comps.size());
    mlog1("here");
    MyComponent & comp = *(comps[i]);
    mlog("pipe", comp.timeStepPipe);
    continue;
    int lastIndex = comp.pipes.size()-1;
    
    //set the next initSet
    comp.pipes.at(lastIndex).evaluate_t(comp.initSet,
        settings.step_end_exp_table);
        
    //mlog("last", (*it)->pipes.at(lastIndex));
    //mlog("vars", (*it)->varIndexes);
    //mlog("evaluated", (*it)->initSet);
  }
  exit(0);
  mrestore(old);
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
  mreset(old);
  mdisable();
  mlog1("sw set <");
  minc();
  TaylorModelVec orig = comp->initSet;
  //comp->initSet = shrinkWrapTMV(comp->initSet, factor);
  
  mlog("input", all.swInput);
  TaylorModelVec wrapped = shrinkWrapTMV(all.swInput, factor);
  mlog1("");
  mlog("wrapped", wrapped);
  
  mlog("compVars", all.compVars);
  
  TaylorModelVec newCompSet;
  
  for(int i = 0; i < comp->compVars.size(); i++) {
    //variable in component
    int varIndex = comp->compVars.at(i);
    //index of that variable in call component
    int indexInAll = find(all.compVars.begin(), all.compVars.end(), varIndex) - 
        all.compVars.begin();
    mlog1(sbuilder() << "index in all: " << indexInAll);
    mlog1(sbuilder() << "init var: " << varIndex);
    //using compVars, cause of the assumption that initial conditions 
    mlog("tm paras", comp->allTMParams);
    mlog("c1", wrapped.tms.at(indexInAll));
    TaylorModel tm = wrapped.tms.at(indexInAll).transform(comp->allTMParams);
    mlog("c2", tm);
    mlog1(sbuilder() << "tm paramCount1: " << wrapped.tms[0].getParamCount());
    mlog1(sbuilder() << "tm paramCount2: " << tm.getParamCount());
    
    newCompSet.tms.push_back(tm);
  }
  
  
  //comp->initSet is the old set, newCompSet is the shrink wrapped one
  if(true) { // TODO add compiler flag
    if(comp->initSet.compare(newCompSet, domain) == false) {
      throw std::runtime_error("something decreased after shrink wrapping");
    }
  }
  
  
  mlog("old", comp->initSet);
  comp->initSet = newCompSet;
  mlog("new", newCompSet);
  
  mdec();
  mlog1("sw set >");
  mrestore(old);
}

void my_handler (const char * reason, const char * file, int line, 
      int gsl_errno) {
  throw IntegrationException(sbuilder() << reason);
}
  
double shrinkWrap(MyComponent & component, vector<Interval> domain, 
     vector<Interval> step_end_exp_table) {
  mreset(old);
  mdisable();
  mlog1("sw <");
  minc();
  double s = -1, t = -1, d = -1, q;
  
  //TaylorModelVec last = component.pipes.at(component.pipes.size() - 1);
  //TaylorModelVec next;
  //last.evaluate_t(next, step_end_exp_table);
  
  TaylorModelVec next = component.swInput;
  
  //mlog("next", next);
  TaylorModelVec orig(next);
  mlog("orig", orig);
  
  int varSize = orig.tms.size();
  int paramCount = next.tms.at(0).getParamCount();
  
  
  //remove the constant part
  next.rmConstant(); 
  
  //mlog("noconst", next);
  
  
  //mlog1(sbuilder() << d);
  
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
  mlog("M", m);
  
  
  //find the inverse of m
  Matrix inv(varSize);
  
  gsl_error_handler_t *old_handler = gsl_set_error_handler (&my_handler);
  try{
    m.inverse(inv);
  }catch(IntegrationException& e) {
    mreset(old2);
    mlog("M", m);
    throw e;
  }
  gsl_set_error_handler (old_handler);
  
  mlog1("M^-1");
  mlog("M inv", inv);
  
  
  //find the nonlinear part of m^-1 * TM
  vector<TaylorModel> linUnitTMs;
  for(int row = 0; row < varSize; row++) {
    //mlog1(sbuilder() << "row:" << row);
    
    //^m-1 * TM part
    TaylorModel result;
    
    for(int col = 0; col < varSize; col++) {
      //mlog1(sbuilder() << "(" << row << "," << col <<")");
      TaylorModel temp;
      next.tms.at(col).mul(temp, inv.get(row,col));
      result.add_assign(temp);
    }
    mlog("m^-1 * TM", result);
    
    linUnitTMs.push_back(result);
  }
  mlog("linear unit TMs", linUnitTMs);
  
  
  
  vector<TaylorModel> nlTMs;
  for(int row = 0; row < varSize; row++) {
    TaylorModel result = linUnitTMs.at(row);
    //substract the linear part (by adding -1*linear part)
    vector<Interval> coefs;
    coefs.push_back(0); //time
    for(int i = 0; i < paramCount - 1; i++) { //-1 because time is counted too
      //mlog1(i);
      if(i == row) // TODO assuming that varIndex == paramIndex
        coefs.push_back(-1);
      else
        coefs.push_back(0);
    }
    mlog("coefs", coefs);
    TaylorModel tm(Polynomial(coefs), Interval(0));
    mlog("-1*lin", tm);
    result.add_assign(tm);
    //mlog("res2", result);
    nlTMs.push_back(result);
  }
  
  
  //nullify remainder (so you get only nonlinear part)
  for(int i = 0; i < nlTMs.size(); i++) {
    nlTMs.at(i).remainder = Interval(0);
  }
  
  TaylorModelVec nlPart(nlTMs);
  
  
  mlog("nlPart", nlPart);
  
  //mlog("C^-1 * ", tms);
  //finding d
  for(vector<TaylorModel>::iterator it = linUnitTMs.begin(); 
      it < linUnitTMs.end(); it++) {
    //mlog1(it->remainder.toString());
    //mlog1(sbuilder() << it->remainder.mag());
    if(it->remainder.mag() > d)
      d = it->remainder.mag();
  }
  mlog1(sbuilder() << "d: " << d);
  
  
  
  //finding s
  //mlog("nl", nlPart);
  //bound the nonlinear part
  vector<Interval> nlRange;
  
  nlPart.polyRange(nlRange, domain);
  
  //mlog("nl", nlRange);
  for(vector<Interval>::iterator it = nlRange.begin();
      it < nlRange.end(); it++) {
    double mag = it->mag();
    s = max(s,mag);
  }
  mlog1(sbuilder() << "s: " << s);
  
  //mlog("params", component.allTMParams);
  
  //find t
  for(int param = 1; param < paramCount; param++) {
    TaylorModelVec pDer;
    nlPart.derivative(pDer, param); //derivative wrt to parameter
    //mlog1(param);
    //mlog("pDer", pDer);
    
    //bound of the derivative
    vector<Interval> pDerRange;
    pDer.polyRange(pDerRange, domain);
    //mlog("pd", pDerRange);
    for(vector<Interval>::iterator it2 = pDerRange.begin();
        it2 < pDerRange.end(); it2++) {
      double mag = it2->mag();
      //mlog1(sbuilder() << mag);
      t = max(t,mag);
    }
  }
  mlog1(sbuilder() << "t: " << t);
  
  int n = paramCount - 1;
  q = 1 + d*(1/((1-(n-1)*t)*(1-s))); //TODO not safe floating point rounding
  //mforce(sbuilder() << "f(t): " << (1-(n-1)*t));
  //mforce(sbuilder() << "f(s): " << (1-s));
  //mforce(sbuilder() << (1/((1-s))));
  if((1 - s) <= 0) {
    throw IntegrationException(sbuilder() << 
        "Map is not shrinkable (s = " << s << ")");
  }
  if((1 - n*t) <= 0) {
    throw IntegrationException(sbuilder() << 
        "Map is not shrinkable (n = " << n << ", t = " << t << ")");
  }
  
  mlog1(sbuilder() << "1 - s: " << (1-s));
  mlog1(sbuilder() << "1 - nt: " << (1-n*t));
  mlog1(sbuilder() << "q:" << q);
  mdec();
  mlog1("sw >");
  mrestore(old);
  return q;
}


void parametrizeVars(TaylorModelVec & tmv, vector<int> varsToIntroduce, 
     int paramCount) {
  mreset(old);
  mdisable();
  mlog1("parametrizing <");
  minc();
  //int paramCount = tmv.tms[0].getParamCount();
  for(int i = 0; i < varsToIntroduce.size(); i++) {
    int var = varsToIntroduce[i];
    //mlog1(sbuilder() << "var: " << var);
    TaylorModel & tm = tmv.tms[var];
    
    Interval & rem = tm.remainder;
    
    Interval inf;
    rem.inf(inf);
    Interval sup;
    rem.sup(sup);
    
    Interval shift = inf + sup;
    shift.div_assign(2);
    
    //mlog1(sbuilder() << "shift: " << shift.toString());
    //mlog1(sbuilder() << "inf: " << inf.toString());
    //mlog1(sbuilder() << "sup: " << sup.toString());
    //mlog1(sbuilder() << "tm: " << tm.remainder.toString());
    
    
    tm.expansion.add_assign(Monomial(shift, paramCount));
    Interval coef = (rem - shift).sup();
    
    //mlog1(sbuilder() << "coef: " << coef.toString());
    
    vector<int> degs;
    //mlog1(tm.getParamCount());
    for(int j = 0; j < paramCount; j++) {
      if(j == var + 1)
        degs.push_back(1);
      else
        degs.push_back(0);
    }
    
    tm.expansion.add_assign(Monomial(coef, degs));
    tm.remainder = Interval();
  }
  mdec();
  mlog1("parametrizing >");
  mrestore(old);
}

void introduceParam(MyComponent & comp, vector<MyComponent *> comps, 
      vector<Interval> & domain) {
  mreset(old);
  mdisable();
  mlog1("introducing <");
  minc();
  TaylorModelVec tmv = comp.swInput;
  
  
  //mlog("last", last);
  mlog("original", tmv);
  
  vector<int> varsToIntroduce;
  for(int i = 0; i < comps.size(); i++) {
    vector<int> & compVars = comps[i]->varsToBeIntroduced;
    varsToIntroduce.insert(varsToIntroduce.end(), 
        compVars.begin(), compVars.end());
    mlog("vs: ", comps[i]->varsToBeIntroduced);
    comps[i]->varsToBeIntroduced.clear();
  }
  mlog("introducing", varsToIntroduce);
  
  //no need to add more parameters now
  //int oParamCount = tmv.tms[0].getParamCount();
  
  //TaylorModelVec padded = tmv.addNParams(varsToIntroduce.size());
  TaylorModelVec padded = TaylorModelVec(tmv);
  int paramCount = padded.tms[0].getParamCount();
  
  
  if(varsToIntroduce.size() != 0) {
    parametrizeVars(padded, varsToIntroduce, paramCount);
  }
  
  mlog("introduced", padded);
  comp.swInput = padded;
  
  
  
  if(true) { // TODO add compiler flag
    if(tmv.compare(comp.swInput, domain) == false) {
      throw std::runtime_error("something decreased after introducing param");
    }
  }
  mdec();
  mlog1("introducing >");
  mrestore(old);
}

double applyShrinkWrapping(MyComponent & all, vector<Interval> domain, 
      vector<Interval> step_end_exp_table, vector<MyComponent *> comps,
      MySettings & settings) {
  mreset(old);
  mdisable();
  mlog1("applying sw <");
  minc();
  
  clock_t start = clock();
  
  all.remapLastFlowpipe();
  
  TaylorModelVec last = all.pipes.at(all.pipes.size() - 1);
  TaylorModelVec tmv;
  last.evaluate_t(tmv, step_end_exp_table);
  all.swInput = tmv;
  
  mlog("swInput", all.swInput);
  //introduces a parameter to swInput
  
  
  introduceParam(all, comps, domain);
  mlog1(all.swInput.tms[0].getParamCount());
  /*
  mlog1(all.compMappers.size());
  mlog("0", all.compMappers.at(0));
  mlog("1", all.compMappers.at(1));
  mlog("2", all.compMappers.at(2));
  mlog("3", all.compMappers.at(3));
  mlog("4", all.compMappers.at(4));
  */
  mlog("swI aft", all.swInput);
  
  
  //mlog("all", all.pipes.at(all.pipes.size() - 1));
  //vector<TaylorModelVec> p2 = all.dependencies.at(0)->pComp->pipes;
  //mlog("p2", p2.at(p2.size() - 1));
  
  //calculate the shrink wrapping factor for swInput
  double swQ = shrinkWrap(all, domain, step_end_exp_table);
  
  mlog1(sbuilder() << "swQ: " << swQ);
  
  //use the computed shrink wrapping factor to modify the initial sets
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    //sets comp.initSet from q and all.swInput
    shrinkWrapSet(all, *it, swQ, domain);
  }
  
  clock_t end = clock();
  double singleSW = double(end - start) / CLOCKS_PER_SEC;
  settings.writer->swTime += singleSW;
  //mforce(sbuilder() << "single: " << singleSW);
  //mforce(sbuilder() << (end - start));
  
  mdec();
  mlog1("applying sw >");
  mrestore(old);
}




Transformer::Transformer(bool isPreconditioned, bool isWrapper) : 
      isPreconditioned(isPreconditioned), isWrapper(isWrapper) {
}

ShrinkWrapper::ShrinkWrapper(ShrinkWrappingCondition *swChecker) : 
      Transformer(false, true) {
  this->swChecker = swChecker;
}

PreconditionedTransformer::PreconditionedTransformer() : Transformer(true, false) {
}

QRTransformer::QRTransformer() : 
    //is preconditioned, not shrink wrapping
    Transformer(true, false) {
}

NullTransformer::NullTransformer() : Transformer(false, false) {
}
IdentityTransformer::IdentityTransformer() : PreconditionedTransformer() {
  //isPreconditioned = true;
}

void IdentityTransformer::getA(Matrix & result, const TaylorModelVec & x0, 
    const int dim) {
  mlog1("getting A");
  for(int i=0; i<dim; ++i) {
    result.set(1, i, i);
	}
	//if(dim == 2)
  //	result.set(2, 0, 1);
}
void IdentityTransformer::getAInv(Matrix & result, const Matrix & A) {
  mlog1("getting A inverse");
  //inverse of identity is also identity, just copy
  for(int i = 0; i < A.rows(); i++) {
    result.set(1, i, i);
  }
}

void ShrinkWrapper::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  mlog1("sw transforming");
  exit(0);
  evaluateStepEnd(comps, settings, true);
  if(swChecker->checkApplicability(comps, settings.estimation)) {
    mforce("wrapping");
    applyShrinkWrapping(all, settings.domain, settings.step_end_exp_table, comps,
        settings);
  }
}

TaylorModelVec getUnitTmv(int varCount) {
  vector<TaylorModel> tms;
  mlog1(varCount);
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
  
	mlog("ret", ret);
  return ret;
}


void precond(TaylorModelVec & tmv, MySettings & settings, 
      MyComponent & all) {
  mreset(old);
  //mdisable();
  mlog1("precond <");
  minc();
  
  mlog1("precond");
  
  mlog("tmv", tmv);
  
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
  mlog("range", previousRange);
  
  mlog1(sbuilder() << "paramCount: " << paramCount);
  mlog1(sbuilder() << "varCount: " << varCount);
  Matrix A(varCount,varCount), invA(varCount, varCount);
  
  if(true) { //use QR preconditioning
  	preconditionQR(A, tmv, varCount, varCount+1);
    A.transpose(invA);
  } else {
    mforce("not implementated");
  }
  mlog("A", A);
  mlog("A inv", invA);
  
  vector<Interval> center;
  tmv.constant(center);
  mlog("center", center);
  
  TaylorModelVec c(tmv), lin, nl, rem;
  tmv.getParts(c, lin, nl, rem);
  
  mlog("c", c);
  mlog("lin", lin);
  mlog("nl", nl);
  mlog("rem", rem);
  
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
  mlog("rightComp", rightComp);
  mlog("rightBound", rightBound);
  
  //TODO shift this to be around 0?
  
  Matrix S(varCount, varCount);
	Matrix invS(varCount, varCount);
  for(int i = 0; i < varCount; i++) {
    Interval intMag;
    rightBound[i].mag(intMag);
    double dMag = intMag.sup();
    if(intMag.subseteq(ZERO_INTERVAL) == false) {
		  S.set(dMag, i, i);
		  invS.set(1/dMag, i, i);
	  }
	  else {
		  S.set(0, i, i);
		  invS.set(0, i, i); //doesn't matter if 0 or 1
		}
  }
  
  mlog("S", S);
  mlog("S inv", invS);
  
  rightComp.linearTrans_assign(invS);
	rightComp.intEvalNormal(rightBound, settings.step_end_exp_table);
  mlog("rightComp", rightComp);
  mlog("rightBound", rightBound);
  
	//left part is made of constant part and desired linear part
	TaylorModelVec left(c);
	TaylorModelVec desiredLinearPart;
	unitTmv.linearTrans(desiredLinearPart, A * S);
	left.add_assign(desiredLinearPart);
	
	PrecondModel *pre = new PrecondModel(left, rightComp);
	
	mlog("left", left);
  mlog("right", rightComp);
	
  all.pipePairs.push_back(pre);
  
  
  mdec();
  mlog1("precond >");
  mrestore(old);
}

void PreconditionedTransformer::precond2(TaylorModelVec & badLeft, 
    MySettings & settings, MyComponent & all) {
  mreset(old);
  mdisable();
  mlog1("precond2 <");
  minc();
  
  mlog("badLeft", badLeft);
  vector<Interval> previousRange;
  TaylorModelVec previousRight = all.unpairedRight;
  
  //number of taylor model parameters
  int paramCount = badLeft.tms[0].getParamCount();
  //number of system variables (in the component)
  int varCount = badLeft.tms.size();
  
  //linear TMV where TM[i] = 1*var_i
  TaylorModelVec unitTmv = getUnitTmv(varCount);
  
  tstart(tr_part_range1);
	previousRight.polyRange(previousRange, settings.domain);
	tend(tr_part_range1);
  mlog("range", previousRange);
  
  mlog1(sbuilder() << "paramCount: " << paramCount);
  mlog1(sbuilder() << "varCount: " << varCount);
  Matrix A(varCount,varCount), invA(varCount, varCount);
  
  tstart(tr_part_matrix);
  getA(A, badLeft, varCount);
  getAInv(invA, A);
  tend(tr_part_matrix);
  
  mlog("A", A);
  mlog("A inv", invA);
  
  
  //center the remainder and push shift to constant part of the polynomial
  badLeft.centerRemainder();
  
  tstart(tr_part_partitioning);
  //partition left model into parts: constant, linear, nonlinear and remainder
  TaylorModelVec c, lin, nl, rem;
  badLeft.getParts(c, lin, nl, rem);
  tend(tr_part_partitioning);
  
  mlog("c", c);
  mlog("lin", lin);
  mlog("nl", nl);
  mlog("rem", rem);
  
  
  tstart(tr_part_part1);
  //InvA o non-constant-part of left model
  nl.linearTrans_assign(invA);
  rem.linearTrans_assign(invA);
  lin.linearTrans_assign(invA);
  
  
  //gather parts that are moving from left to right
  TaylorModelVec leftToRight(lin);
  leftToRight.add_assign(nl);
  leftToRight.add_assign(rem);
  
  
  //compose new parts from left model with older right model
  TaylorModelVec rightComp;
	leftToRight.insert_ctrunc(rightComp, previousRight, previousRange, settings.domain, settings.order, 0);
  
  tend(tr_part_part1);
	
  //calculate the bounds for the new right model
  vector<Interval> rightBound;
	rightComp.intEvalNormal(rightBound, settings.step_end_exp_table);
  mlog("rightComp", rightComp);
  mlog("rightBound", rightBound);
  
  //TODO shift this to be around 0?
  
  tstart(tr_part_part2);
  //calculate scaling matrix S and it's inverse
  Matrix S(varCount, varCount);
	Matrix invS(varCount, varCount);
  for(int i = 0; i < varCount; i++) {
    Interval intMag;
    rightBound[i].mag(intMag);
    double dMag = intMag.sup();
    if(intMag.subseteq(ZERO_INTERVAL) == false) {
		  S.set(dMag, i, i);
		  invS.set(1/dMag, i, i);
	  } else {
		  S.set(0, i, i);
		  invS.set(0, i, i); //doesn't matter if 0 or 1
		}
  }
  
  mlog("S", S);
  mlog("S inv", invS);
  
  //apply inverse of S to right model
  rightComp.linearTrans_assign(invS);
	rightComp.intEvalNormal(rightBound, settings.step_end_exp_table);
  mlog("rightComp", rightComp);
  mlog("rightBound", rightBound);
  
  
	//left part is made of constant part and with scaling
	TaylorModelVec left = TaylorModelVec(c);
	TaylorModelVec desiredLinearPart;
	unitTmv.linearTrans(desiredLinearPart, A * S);
	left.add_assign(desiredLinearPart);
	
	//pSerializer->add(left, "left_after_precond");
	
	all.unpairedRight = rightComp;
  all.initSet = left;
	
  tend(tr_part_part2);
  mdec();
  mlog1("precond >");
  mrestore(old);
}


void PreconditionedTransformer::precond3(TaylorModelVec & leftStar, 
    MySettings & settings, MyComponent & all) {
  mreset(old);
  mdisable();
  mlog1("precond3 <");
  minc();
  
  tstart(tr_part_all);
  tstart(tr_part1);
  
  int paramCount = leftStar.tms[0].getParamCount();
  int varCount = leftStar.tms.size();
  
  int rangeDim = varCount;
  	
	//center remainder around zero and push the shift to constant
	leftStar.centerRemainder();
  
  vector<Interval> constantVec;
	leftStar.constant(constantVec);
	TaylorModelVec c0(constantVec, paramCount);

  //remove constant part from old part
	leftStar.rmConstant();
	
	mlog("left", leftStar);
	
	
	
	Matrix A(varCount,varCount), invA(varCount, varCount);
  tend(tr_part1);
  tstart(tr_part_matrix);
  getA(A, leftStar, varCount);
  getAInv(invA, A);
  
  
  // if we want left model to be c0 + A id
  // then we need to move A^-1 * <old left without constant part> to right model
  TaylorModelVec leftToRight;
  //leftStar.linearTrans(leftToRight, invA);
  leftToRight = leftStar; //TODO remove this
  
  tend(tr_part_matrix);
  
  tstart(tr_part2);
  
  mlog("leftToRight", leftToRight); 
    
  TaylorModelVec & prevRight = all.unpairedRight;
  
	vector<Interval> prevRightPolyRange;
	prevRight.polyRangeNormal(prevRightPolyRange,
	    settings.step_end_exp_table);
	mlog("range", prevRightPolyRange);
	
	
	//current right is leftToRight o previousRight
	TaylorModelVec currentRight;
	leftToRight.insert_ctrunc_normal(currentRight, prevRight,
	    prevRightPolyRange, settings.step_end_exp_table, paramCount, 
	    settings.order, settings.cutoff);
  mlog("currentRight", currentRight);
  
  vector<Interval> currentRightRange;
  currentRight.intEvalNormal(currentRightRange, settings.step_end_exp_table);
  mlog("currentRightRange", currentRightRange);
  
  //calculate scaling matrix from current right model
  Matrix S(varCount, varCount);
	Matrix invS(varCount, varCount);
	for(int i=0; i<varCount; i++) {
		Interval intSup;
		currentRightRange[i].mag(intSup);
		if(intSup.subseteq(ZERO_INTERVAL)) {
			S.set(0,i,i);
			invS.set(1,i,i);
		} else {
			double dSup = intSup.sup();
			S.set(dSup, i, i);
			invS.set(1/dSup, i, i);
			//assuming that this is overriden to be domain
			currentRightRange[i] = UNIT_INTERVAL;
		}
	}
  
  
  //apply scaling
  A = A * S;
	invA = invS * invA; //not needed actually
	

  //apply scaling to right model
	currentRight.linearTrans_assign(invS);
	currentRight.cutoff_normal(settings.step_end_exp_table, settings.cutoff);
	
	

  //construct left matrix using A coefficients and constant part
	Matrix leftLinCoefs(rangeDim, rangeDim+1);
	for(int i=0; i<rangeDim; i++) {
		for(int j=0; j<rangeDim; j++) {
			leftLinCoefs.set( A.get(i,j) , i, j+1);
		}
	}
	TaylorModelVec newLeft(leftLinCoefs);
  newLeft.add_assign(c0);
  
  
	mlog("newLeft", newLeft);
	mlog("right", currentRight);
	
	
	all.unpairedRight = currentRight;
  all.initSet = newLeft;
  
  tend(tr_part2);
  tend(tr_part_all);
  
  
  mdec();
  mlog1("precond3 >");
  mrestore(old);
}


void makeCompInitSets(MyComponent & all, MyComponent * comp) {
  mreset(old);
  mdisable();
  mlog1("makeCompInitSets <");
  minc();

  TaylorModelVec newCompSet;
  
  for(int i = 0; i < comp->compVars.size(); i++) {
    //variable in component
    int varIndex = comp->compVars.at(i);
    //index of that variable in call component
    int indexInAll = 
        find(all.compVars.begin(), all.compVars.end(), varIndex) - 
        all.compVars.begin();
    mlog1(sbuilder() << "comp var: " << varIndex);
    mlog1(sbuilder() << "index in all: " << indexInAll);
    //using compVars, cause of the assumption that initial conditions 
    mlog("tm params", comp->allTMParams);
    mlog("c1", all.initSet.tms.at(indexInAll));
    TaylorModel tm = all.initSet.tms.at(indexInAll).
        transform(comp->allTMParams);
    mlog("c2", tm);
    mlog1(sbuilder() << "tm paramCount1: " << 
        all.initSet.tms[0].getParamCount());
    mlog1(sbuilder() << "tm paramCount2: " << tm.getParamCount());
    
    newCompSet.tms.push_back(tm);
  }
  mlog("init1", comp->initSet);
  comp->initSet = newCompSet;
  mlog("init2", comp->initSet);
  
  mdec();
  mlog1("makeCompInitSets >");
  mrestore(old);
}


void QRTransformer::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  mreset(old);
  mdisable();
  mlog1("qr transforming <");
  minc();
  evaluateStepEnd(comps, settings, true);
  
  //remaps all the last pipes to system flowpipes at all component
  all.remapLastFlowpipe();
  
  TaylorModelVec last = all.pipes.at(all.pipes.size() - 1);
  TaylorModelVec tmv;
  last.evaluate_t(tmv, settings.step_end_exp_table);
  
  
  
  precond(tmv, settings, all);
  mlog("tmv", tmv);
  
  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    makeCompInitSets(all, *it);
  }
  
  
  mdec();
  mlog1("qr transforming >");
  mrestore(old);
}

void NullTransformer::transform(MyComponent & all, vector<MyComponent *> & comps, 
      MySettings & settings) {
  //mlog1("null transforming");
  evaluateStepEnd(comps, settings, true);
}
void IdentityTransformer::transformFullSystem(MyComponent & all, 
      vector<MyComponent *> & comps, MySettings & settings) {
  mreset(old);
  //mdisable();
  mlog1("id transforming <");
  minc();
  //evaluateStepEnd(comps, settings, false);
  
  //TODO evaluate first, then remap maybe
  tstart(tr_remap1);
  //remaps all the last pipes to system flowpipes at all component
  all.remapLastFlowpipe();
  tend(tr_remap1);
  
  tstart(tr_pipe);
  TaylorModelVec tsp = all.timeStepPipe;
  //mlog("tsp", tsp);
  //mlog("upright", all.unpairedRight); 
  
  mlog("unpairedRight", all.unpairedRight);
  
  all.pipePairs.push_back(new PrecondModel(tsp, all.unpairedRight));
  tend(tr_pipe);
  
  
  tstart(tr_eval);
  TaylorModelVec leftStar;
  tsp.evaluate_t(leftStar, settings.step_end_exp_table);
  tend(tr_eval);
  
  //pSerializer->add(leftStar, "leftStar");
  
  
  //pSerializer->serialize();
  
  mlog1(sbuilder() << "size: " << leftStar.tms.size());
  //mlog("leftStar", leftStar);
  
  tstart(tr_precond);
  precond3(leftStar, settings, all);
  tend(tr_precond);
  //mlog("upright", all.unpairedRight);
  
  
	//mlog("left", all.initSet);
  //mlog("right", all.unpairedRight);
  
  
  //mlog("pair.left", all.lastPre()->left);
  //mlog("pair.right", all.lastPre()->right);

  tstart(tr_remap2);  
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    makeCompInitSets(all, *it);
    //mlog("ci", (*it)->initSet);
  }
  tend(tr_remap2);
  
  mdec();
  mlog1("id transforming >");
  mrestore(old);
  exit(0);
  //evaluateStepEnd(comps, settings);
}

void IdentityTransformer::preconditionSingleComponent(MyComponent *comp,
      MySettings & settings) {
  if(comp->isPreconditioned)
    return;
  mlog1("single <");
  minc();
  for(int i = 0; i < comp->dependencies.size(); i++) {
    MyComponent *pComp = comp->dependencies[i]->pComp;
    preconditionSingleComponent(pComp, settings);
  }
  
  
  int varCount = comp->compVars.size();
  int paramCount = comp->allTMParams.size() + 1; //+1 for time
  
  
  mlist("varIndexes", comp->varIndexes);
  mlist("linkVars", comp->linkVars);
  mlist("compVars", comp->compVars);
  mlist("solveIndexes", comp->solveIndexes);
  mlist("allParams", comp->allTMParams);
  
  mlog1(varCount);
  
  TaylorModelVec leftStar;
  TaylorModelVec right;
  
  for(int i = 0; i < comp->compVars.size(); i++) {
    int variable = comp->compVars[i];
    bool isSolve = find(comp->solveIndexes.begin(), comp->solveIndexes.end(),
        i) != comp->solveIndexes.end();
    mlog1(sbuilder() << "var: " << variable << ", solve: " << isSolve);
    if(isSolve) {
      leftStar.tms.push_back(comp->initSet.tms[i]);
      right.tms.push_back(comp->unpairedRight.tms[i]);
    } else {
      mlog1("need to transform");
      
      leftStar.tms.push_back(comp->remapSingleInitTM(variable));
      right.tms.push_back(comp->remapSingleRightTM(variable));
    }
  }
  mlog("left*", leftStar);
  mlog("right", right);
  
  
  //might make sense to cut this out from big A
  Matrix A(varCount, varCount), invA(varCount, varCount);
  
  
  //center remainder around zero and push the shift to constant
	leftStar.centerRemainder();
  
  vector<Interval> constantVec;
	leftStar.constant(constantVec);
	TaylorModelVec c0(constantVec, paramCount);

  //remove constant part from old part
	leftStar.rmConstant();
	
  mlog("l*", leftStar);
  
  //TODO previous variables are unnecessary
  getA(A, leftStar, varCount);
  getAInv(invA, A);
  
  mlog("A", A);
  
  TaylorModelVec leftToRight;
  //leftStar.linearTrans(leftToRight2, invA);
  leftToRight = leftStar; //TODO remove this (when not using id preconditioning)
  
  mlog("ltr", leftToRight);
  
  TaylorModelVec & prevRight = right;
  mlog("prevRight", prevRight);
  
	vector<Interval> prevRightPolyRange;
	prevRight.polyRangeNormal(prevRightPolyRange,
	    settings.step_end_exp_table);
	mlog("range", prevRightPolyRange);
	
	//current right is leftToRight o previousRight
	TaylorModelVec currentRight;

	leftToRight.insert_ctrunc_normal(currentRight, prevRight,
	    prevRightPolyRange, settings.step_end_exp_table, paramCount, 
	    settings.order, settings.cutoff);
  mlog("currentRight", currentRight);
  mlog1(paramCount);
  
  vector<Interval> currentRightRange;
  currentRight.intEvalNormal(currentRightRange, settings.step_end_exp_table);
  mlog("currentRightRange", currentRightRange);
  
  for(int i = 0; i < comp->solveIndexes.size(); i++) {
    int varInComp = comp->solveIndexes[i];
    int actualVar = comp->varIndexes[varInComp];
    
    //pick the ith row from A and make Taylor Model
    RowVector rowVec(varCount + 1);
    Matrix rowMatrix(1, varCount);
	  for(int j=0; j < varCount; j++) {
      rowVec.set(A.get(varInComp, j), j + 1);
      rowMatrix.set(A.get(varInComp, j), 0, j);
	  }
	  mlog("rowVec", rowVec);
	  mlog("rowMatrix", rowMatrix);
	  TaylorModel newLeftRow(rowVec);
	  logger.log("newLeftRow", newLeftRow);
    
    mlog1(sbuilder() << "varInComp: " << varInComp
        << ", actualVar: " << actualVar);
  }
  
  if(comp->compVars.size() > 1)
    exit(0);
  
  comp->isPreconditioned = true; //TODO make sure to reset it
  mdec();
  mlog1("single >");
}

void IdentityTransformer::transform(MyComponent & all, 
      vector<MyComponent *> & comps, MySettings & settings) {
  mreset(old);
  //mdisable();
  mlog1("id transforming by component <");
  minc();
  
  for(int i = 0; i < comps.size(); i++) {
    MyComponent c1 = *comps[i];
    preconditionSingleComponent(comps[i], settings);
    mlog("unpaired r", c1.unpairedRight);
      
    if(c1.linkVars.size() > 10)
      exit(0);
  }
  //evaluate all components at t=t_end (left*)
  
  //make new flowpipe in pipe pairs from unpaired right
  
  //precondition all the components
  
  //set the initial sets for new component
  
  
  mdec();
  mlog1("id transforming by component>");
  mrestore(old);
  exit(0);
  //evaluateStepEnd(comps, settings);
}

void PreconditionedTransformer::addInfo(vector<string> & info) {
  int old = logger.reset();
  logger.disable();
  mlog1("adding");
  mlog1(sbuilder() << timeLookup["tr_precond"]);
  
  taddToInfo("remap 1", tr_remap1, info);
  taddToInfo("evaluate t", tr_eval, info);
  taddToInfo("precond time", tr_precond, info);
  taddToInfo("remap 2", tr_remap2, info);
  taddToInfo("int time", sc_integrate, info);
  
  logger.restore(old);
}

void QRTransformer::addInfo(vector<string> & info) {
  throw std::runtime_error("don't call add info on QR yet");
}
void NullTransformer::addInfo(vector<string> & info) {
  throw std::runtime_error("don't call add info on Null yet");
}
void ShrinkWrapper::addInfo(vector<string> & info) {
  throw std::runtime_error("don't call add info on SW yet");
}




