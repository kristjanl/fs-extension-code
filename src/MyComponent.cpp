#include "MyComponent.h"

using namespace std;


//TODO remove! (duplicate of Transformer function
namespace MyComponentRemove{
  TaylorModelVec getUnitTmv(int varCount) {
    vector<TaylorModel> tms;
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
    
	  //mlog("ret", ret);
    return ret;
  }
  
  //TODO make this smarter (to support any given initial TM)
  TaylorModelVec getNVarMParam(int varCount, int paramCount) {
    if(varCount == 0)
      return TaylorModelVec();
    int shift = paramCount - varCount;
    vector<TaylorModel> tms;
    Matrix coefs(varCount, paramCount + 1);
    for(int i = 0; i < varCount; i++) {
      for(int j = 0; j < paramCount; j++) {
        if(i == j) {
          coefs.set(1, i, j + 1 + shift);
        }
      }
    }
    TaylorModelVec ret(coefs);
    return ret;
  }
  
  //duplicated in Utils
  int findPos(int value, vector<int> *v) {
    return find(v->begin(), v->end(), value) - v->begin();
  }
  //duplicated in Utils
  int isIn(int value, vector<int> *v) {
    return find(v->begin(), v->end(), value) != v->end();
  }
}

MyComponent::MyComponent(vector<int> vs, vector<int> tps):varIndexes(vs),tpIndexes(tps) {
  isSolved = false;
  isPrepared = false;
  isPreconditioned = false;
  retainEmptyParams = false;
  firstPrecondition = true;
}
MyComponent::MyComponent() {
  isSolved = false;
  isPrepared = false;
  isPreconditioned = false;
  retainEmptyParams = false;
  firstPrecondition = true;
}

void MyComponent::log() {
  mlog1("component <");
  minc();
  stringstream varSs;
  mlog1(varIndexes.size());
  for(int i = 0; i < varIndexes.size(); i++) {
    varSs << varIndexes.at(i) << ", ";
  }
  mlog1(sbuilder() << "varIndexes: " << varSs.str());
  mlog1(sbuilder() << "odes size: " << odes.size());
  {
    int odeIndex = 0;
    for(vector<HornerForm>::iterator it = odes.begin(); it < odes.end(); it++) {
      mlog1(sbuilder() << odeIndex++ << ": " << (*it).toString());
    }
  }
  mlog("compDomain", dom);  
  mlog("compInit", initSet);
  
  for(int i = 0; i < initSet.tms.size(); i++) {
    mlog1(sbuilder() << "tm param size(" << i << "): "
        << initSet.tms.at(i).getParamCount());
  }
  mlog1(sbuilder() << "init set size: " << initSet.tms.size());
  mlog("to be solved variables", varIndexes);
  mlog("comp variables", compVars);
  mlog("all variables", allVars);
  mlog("solve indexes", solveIndexes);
  mlog("new params", tpIndexes);
  mlog("all params", allTMParams);
  mlog1(sbuilder() << "mappers size: " << compMappers.size());
  for(int i = 0; i < compMappers.size(); i++) {
    mlog(
        (sbuilder() << "m" << i << 
        " (link is " << dependencies.at(i)->linkVar <<")"),
        compMappers.at(i));
  }
  mdec();
  mlog1("component >");
  
}

/*
creates mappers for previous components
mappers are of the form [i1,i2,...,in]
where is are either indexes to parameters in the previous components
or -2 to indicate placeholders
[-2,-2,0,-2,1] maps a1*a2 to a3*a5
*/
void MyComponent::prepareMappers() {
  mreset(old);
  mdisable();
  mlog1("preparing mappers <");
  minc();
  vector< vector<int> > mappers;
  

  for(int i = 0; i < dependencies.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    //mlog1(sbuilder() << "link: " << dependencies.at(i)->linkVar);
    //previous component parameters
    vector<int> prevCompPm = pComp->allTMParams;
    //mlog("pm", prevCompPm);
    vector<int> mapper;
    //iterate all the variables in all the components
    for(int j = 0; j < allTMParams.size(); j++) {
      int var = allTMParams.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = 
          find(prevCompPm.begin(), prevCompPm.end(), var);
      
      if(it != prevCompPm.end()) {
        //find position if it was
        int pos = distance(prevCompPm.begin(), it);
        //map old position to new
        mapper.push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper.push_back(PADDING_VARIABLE);
      }
    }
    dependencies.at(i)->mapper = mapper;
    dependencies.at(i)->rightMapper = mapper;
    //mlog("m", mapper);
    mappers.push_back(mapper);
  }
  compMappers = mappers;
  
  //mapper for left TaylorModel
  for(int i = 0; i < dependencies.size(); i++) {
    MyComponent *pComp = dependencies.at(i)->pComp;
    
    //variables in previous component
    vector<int> prevCompVars = pComp->allVars;
    //mlog("pv", prevCompVars);
    vector<int> mapper;
    
    //loop over all current component variables
    for(int j = 0; j < allVars.size(); j++) {
      int var = allVars.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = 
          find(prevCompVars.begin(), prevCompVars.end(), var);
      
      if(it != prevCompVars.end()) {
        //find position if it was
        int pos = distance(prevCompVars.begin(), it);
        //map old position to new
        mapper.push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper.push_back(PADDING_VARIABLE);
      }
    }
    //mlog("mapper", mapper);
    leftMappers.push_back(mapper);
    dependencies.at(i)->leftMapper = mapper;
  }
  
  mdec();
  mlog1("preparing mappers >");
  mrestore(old);
}


vector< vector<int> > MyComponent::previousMappers2() {
  vector<int> allParams;
  //include current component parameters
  for(vector<int>::iterator it = tpIndexes.begin(); 
        it < tpIndexes.end(); it++) {
    allParams.push_back(*it);
  }
   
  mlog("allParams (current)", allParams);
  //gather all indexes of all previous components
  for(int i = 0; i < dependencies.size(); i++) {
    mlog1(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    
    mlog("vi", pComp->varIndexes);
    mlog("pi", pComp->tpIndexes);
    allParams.insert(allParams.end(), 
        pComp->allTMParams.begin(), pComp->allTMParams.end());
  }
  mlog("ai", allParams);
  //sort all indexes
  sort(allParams.begin(), allParams.end());
  //remove duplicate indexes
  allParams.erase( unique( allParams.begin(), allParams.end() ), allParams.end() );
  
  allTMParams = allParams;
  mlog("allParams (with dependencies)", allParams);
  
  
  mlog("ai", allParams);
  vector< vector<int> > mappers;
  
  for(int i = 0; i < dependencies.size(); i++) {
    mlog1(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    mlog1(sbuilder() << "link: " << dependencies.at(i)->linkVar);
    //previous component parameters
    vector<int> prevCompPm = pComp->allTMParams;
    mlog("pm", prevCompPm);
    vector<int> mapper;
    //iterate all the variables in all the components
    for(int j = 0; j < allParams.size(); j++) {
      int var = allParams.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = find(prevCompPm.begin(), prevCompPm.end(), var);
      if(it != prevCompPm.end()) {
        //find position if it was
        int pos = distance(prevCompPm.begin(), it);
        //map old position to new
        mapper.push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper.push_back(PADDING_VARIABLE);
      }
    }
    dependencies.at(i)->mapper = mapper;
    mlog("m", mapper);
    mappers.push_back(mapper);
  }
  //mlog("mapper0", mappers.at(0));
  //mlog("mapper1", mappers.at(1));
  
  //return *mappers;
  return mappers;
}


vector< vector<int> > MyComponent::previousMappers() {
  vector<int> allIndexes;
  
  //gather all indexes of all components
  for(int i = 0; i < previous.size(); i++) {
    allIndexes.insert(allIndexes.end(), 
        previous.at(i).varIndexes.begin(), previous.at(i).varIndexes.end());
  }
  //sort all indexes
  sort(allIndexes.begin(), allIndexes.end());
  //remove duplicate indexes
  allIndexes.erase( unique( allIndexes.begin(), allIndexes.end() ), allIndexes.end() );
  
  //mlog("allIndexes", allIndexes);
  
  vector< vector<int> > *mappers = new vector< vector<int> >();
  
  for(int i = 0; i < previous.size(); i++) {
    //previous component variables
    vector<int> prevCompInd = previous.at(i).varIndexes;
    vector<int> *mapper = new vector<int>();
    //iterate all the variables in all the components
    for(int j = 0; j < allIndexes.size(); j++) {
      int var = allIndexes.at(j);
      
      //find if the variable was in the previous component
      vector<int>::iterator it = find(prevCompInd.begin(), prevCompInd.end(), var);
      if(it != prevCompInd.end()) {
        //find position if it was
        int pos = distance(prevCompInd.begin(), it);
        //map old position to new
        mapper->push_back(pos);
      } else {
        //add a padding flag to variable if it's new to component
        mapper->push_back(PADDING_VARIABLE);
      }
    }
    mappers->push_back(*mapper);
  }
  //mlog("mapper0", mappers->at(0));
  //mlog("mapper1", mappers->at(1));
  
  return *mappers;
}


CompDependency::CompDependency(int link, MyComponent *pComp)
    :linkVar(link), pComp(pComp) {
  //mlog1("creating");
}


vector<int> concateMapper(vector<int> & smaller, vector<int> & bigger) {
  mlog1("concate");
  mlog("smaller", smaller);
  mlog("bigger", bigger);
  
  if(bigger.size() == 0)
    return smaller;
  
  vector<int> ret;
  //look into all the parameters in mapper
  for(int i = 0; i < bigger.size(); i++) {
    int param = bigger[i];
    
    //retain padding parameters
    if(param == PADDING_VARIABLE) {
      ret.push_back(PADDING_VARIABLE);
      continue;
    }
    //else use the dependency parameter
    ret.push_back(smaller[param]);
  }
  return ret;
}

void MyComponent::addDependency(int linkVar, MyComponent *pComp) {
  dependencies.push_back(new CompDependency(linkVar, pComp));
}

void MyComponent::addVar(int var) {
  //mlog1("adding");
  varIndexes.push_back(var);
  compVars.push_back(var);
}

void addEmptyTM(vector<TaylorModelVec> & pipes, int dim) {
  throw invalid_argument("don't use addEmptyTM"); 
  TaylorModelVec temp;
  TaylorModel t;
  for(int i = 0; i < dim; i++) {
    temp.tms.push_back(t);
  }
  pipes.push_back(temp);
}

TaylorModelVec getEmptyTMV(int dim) {  
  TaylorModelVec temp;
  TaylorModel t;
  for(int i = 0; i < dim; i++) {
    temp.tms.push_back(t);
  }
  return temp;
}

void MyComponent::remapIVP(TaylorModelVec tmv, const vector<HornerForm> & ode, 
      vector<Interval> domain) {
  mreset(old);
  mdisable();
  mlog1("remapping IV <");
  minc();

  
  mlog("all params", allTMParams);
  mlog("all variables", compVars);
  mlog("link variables", linkVars);
  
  vector<HornerForm> compOdes;
  for(int i = 0; i < compVars.size(); i++) {
    int varIndex = compVars.at(i);
    mlog1(sbuilder() << "v: " << varIndex);
    minc();
    HornerForm hf;
    if(linkVars.end() != find(linkVars.begin(), linkVars.end(), varIndex)) {
      // we don't care about the ode if the variable is a link
    } else {
      try{
        mlog("compVars", compVars);
        mlog1(sbuilder() << "ode: " << ode.at(varIndex).toString());
        hf = (ode.at(varIndex)).transform(compVars);
      }
      catch (exception& e) {
        mreset(old);
        mlog1(e.what());
        mlog("indexes", compVars);
        throw e;
      }
    }
    compOdes.push_back(hf);
    mlog1(hf.toString());
    mdec();
  }
  
  vector<Interval> compDomain;
  compDomain.push_back(domain.at(0)); //time domain
  for(int i = 0; i < allTMParams.size(); i++) {
    int paramIndex = allTMParams.at(i);
    compDomain.push_back(domain.at(paramIndex + 1)); //+1 for time
  }
  
  TaylorModelVec *compInit = new TaylorModelVec();
  for(int i = 0; i < compVars.size(); i++) {
    int varIndex = compVars.at(i);
    mlog1(sbuilder() << "init var: " << varIndex);
    //using compVars, cause of the assumption that initial conditions 
    mlog("tm paras", allTMParams);
    mlog("c1", tmv.tms.at(varIndex));
    TaylorModel tm = tmv.tms.at(varIndex).transform(allTMParams);
    mlog("c2", tm);
    (*compInit).tms.push_back(tm);
  }
  
  odes = compOdes;
  dom = compDomain;
  initSet = *compInit;
  
  mdec();
  mlog1("remapping IV >");
  mrestore(old);
}
  

//prepares all the dependent components
//prepares variables, mappers and intial set
void MyComponent::prepareComponent(TaylorModelVec init, 
    const vector<HornerForm> & ode, vector<Interval> domain) {
  //return if variables have already been prepared
  if(isPrepared) {
    return;
  }
  mreset(old);
  mdisable();
  mlog1("start");
  mlog1("preparing component <");
  minc();
  mlog("vars", varIndexes);
  mlog1(sbuilder() << "dependencies: " << dependencies.size());
  //prepare variables for all the dependencies
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
        it < dependencies.end(); it++) {
    mlog1(sbuilder() << "link: " << (*it)->linkVar);
    MyComponent *pComp = (*it)->pComp;
    pComp->prepareComponent(init, ode, domain);
  }
  prepareVariables(init, ode);
  prepareMappers();
  remapIVP(init, ode, domain);
  timeStepPipe = TaylorModelVec(this->initSet);
  
  //only create a part of right taylor model 
  //unpairedRight = MyComponentRemove::getUnitTmv(allTMParams.size());
  unpairedRight = MyComponentRemove::getUnitTmv(tpIndexes.size());
  unpairedRight = MyComponentRemove::getNVarMParam(tpIndexes.size(),
      allTMParams.size());
  //mforce3(aaaa, "in prepare2", unpairedRight);
  isPrepared = true;
  
  mdec();
  mlog1("preparing component >");
  mrestore(old);
}
  
void MyComponent::prepareVariables(TaylorModelVec init, 
    const vector<HornerForm> & ode) {
  mreset(old);
  mdisable();
  mlog1("preparing variables <");
  minc();
  
  mlog("vars", compVars);
  for(int i = 0; i < compVars.size(); i++) {
    allVars.push_back(compVars[i]);
  }
  for(int i = 0; i < dependencies.size(); i++) {
    MyComponent dep = *dependencies[i]->pComp;
    allVars.insert(allVars.end(), dep.allVars.begin(), dep.allVars.end());
    
    /*version that doesn't sort
    for(int j = 0; j < dep.allVars.size(); j++) {
      bool isNew = find(allVars.begin(), allVars.end(), dep.allVars[j])
          == allVars.end(); 
      mlog1(sbuilder() << "v: " << dep.allVars[j] << ", isNew: " << isNew);
      if(isNew) {
        allVars.push_back(dep.allVars[j]);
      }
    }*/
  }
  
  //need to sort to remove duplicates
  sort(allVars.begin(), allVars.end());
  allVars.erase( 
      unique(allVars.begin(), allVars.end()), allVars.end() );
  //*/
  mlog("all", allVars);
  //add all linking variables
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++) {
    compVars.push_back((*it)->linkVar);
    //TODO can remove and replace with unlinked
    linkVars.push_back((*it)->linkVar);
  }
  
  //find the indexes (relative in the component) of variables which 
  //need to be solved
  {
    vector<int> *v1 = &(varIndexes);
    vector<int> *v2 = &(compVars);
    for(vector<int>::iterator it = v1->begin(); it < v1->end(); it++) {
      //mlog1(*it);
      int pos = find(v2->begin(), v2->end(), *it) - v2->begin();
      //mlog1(sbuilder() << "pos: " << pos);
      solveIndexes.push_back(pos);
    }
  }
  
  
  
  //gather nonzero taylor model parameters from initial conditions
  for(vector<int>::iterator it = varIndexes.begin(); 
      it < varIndexes.end(); it++) {
    mlog1(sbuilder() << "var: " << *it);
    //mlog("tm", init.tms.at(*it));
    vector<int> tmParams = init.tms.at(*it).getParams();
    mlog("params", tmParams);
    
    
    //retain the parameter if started with point initial conditions
    if(tmParams.size() == 0 && retainEmptyParams) {
      varsToBeIntroduced.push_back(*it);
      tmParams.push_back(*it);
    }
    tpIndexes.insert(tpIndexes.end(), tmParams.begin(), tmParams.end());
  }
  sort(tpIndexes.begin(), tpIndexes.end());
  tpIndexes.erase( 
      unique(tpIndexes.begin(), tpIndexes.end()), tpIndexes.end() );
  
  
  vector<int> allParams;
      
  //include current component parameters
  for(vector<int>::iterator it = tpIndexes.begin(); 
        it < tpIndexes.end(); it++) {
    allParams.push_back(*it);
  }
  
  //mlog("allParams (current)", allParams);
  //gather all indexes of all previous components
  for(int i = 0; i < dependencies.size(); i++) {
    //mlog1(sbuilder() << "i: " << i);
    MyComponent *pComp = dependencies.at(i)->pComp;
    
    //mlog("vi", pComp->varIndexes);
    //mlog("pi", pComp->tpIndexes);
    allParams.insert(allParams.end(), 
        pComp->allTMParams.begin(), pComp->allTMParams.end());
  }
  
  //sort all params
  sort(allParams.begin(), allParams.end());
  //remove duplicate params
  allParams.erase( 
      unique( allParams.begin(), allParams.end() ), allParams.end() );
  
  allTMParams = allParams;
  //mlog("allParams (with dependencies)", allParams);
  
  mdec();
  mlog1("preparing variables >");
  mrestore(old);
}  

void MyComponent::remapFlowpipes() {
  vector< vector<int> > mappers = previousMappers2();
  mlog1(sbuilder() << "mappers size: " << mappers.size());
  

  //mapping of dependent component flowpipes to current component
  int depIndex = 0;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++, depIndex++) {
    mlog1(sbuilder() << "depIndex: " << depIndex);
    int link = (*it)->linkVar;
    mlog1(sbuilder() << "link: " << link);
    MyComponent *pComp = (*it)->pComp;
    mlog("dts", pComp->tpIndexes);
    mlog("mapper", (*it)->mapper);
    mlog1(pComp->pipes.size());
    mlog1("here");
    mlog("p", pComp->pipes.at(0));
    
    
    mlog1("here");
    mlog("dcv", pComp->compVars);
    
    
    //position of linking variable in precomputed component
    int dLinkPos = find(pComp->compVars.begin(), pComp->compVars.end(), link) - 
        pComp->compVars.begin();
    //position of linking variable in current component
    int linkPos = find(compVars.begin(), compVars.end(), link) - 
        compVars.begin();
    if(dLinkPos >= pComp->compVars.size()) {
      mlog1(sbuilder() << "dLinkPos: " << dLinkPos);
      mlog("compVars", compVars);
      throw std::invalid_argument("depend link wasn't in compVar");
    }
    if(linkPos >= compVars.size()) {
      throw std::invalid_argument("link wasn't in compVars");
    }
    mlog1(sbuilder() << "dLinkPos: " << dLinkPos);
    mlog1(sbuilder() << "linkPos: " << linkPos);
    
    
    mlog("vars", varIndexes);
    mlog("dvars", pComp->varIndexes); 
      
    //actual mapping of flowpipes
    int pipeIndex = 0;
    for(vector<TaylorModelVec>::iterator pipeIt = pComp->pipes.begin();
        pipeIt < pComp->pipes.end(); pipeIt++, pipeIndex++) {
      mlog1(sbuilder() << "pipeIndex: " << pipeIndex);
      mlog("tm", (*pipeIt).tms.at(dLinkPos));
      mlog("ma", mappers.at(depIndex));
      
      
      TaylorModel transformedLink = 
          (*pipeIt).tms.at(dLinkPos).transform(mappers.at(depIndex));
      pipes.at(pipeIndex).tms.at(linkPos) = transformedLink;
      //mlog("pipe", pipes.at(pipeIndex));
    }
    
  }
}

void MyComponent::remapLastFlowpipe() {
  mreset(old);
  mdisable();
  mlog1("remapping last <");
  minc();
  
  int varSize = varIndexes.size() + dependencies.size();

  //pipes.push_back(getEmptyTMV(varSize));
  timeStepPipe = getEmptyTMV(varSize);
  mlog1(dependencies.size());
  
  if(dependencies.size() == 0) {
    mdec();
    mrestore(old);
    return;
  }
  
  //mapping of dependent component flowpipes to current component
  int depIndex = 0;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++, depIndex++) {
    mlog1(sbuilder() << "depIndex: " << depIndex);
    int link = (*it)->linkVar;
    mlog1(sbuilder() << "link: " << link);
    MyComponent *pComp = (*it)->pComp;
    mlog("dts", pComp->tpIndexes);
    mlog("mapper", (*it)->mapper);
    mlog("rightMapper", (*it)->rightMapper);
    mlog("leftMapper", (*it)->leftMapper);
    
    
    //mlog1("here");
    mlog("dcv", pComp->compVars);
    
    
    //position of linking variable in precomputed component
    int dLinkPos = find(pComp->compVars.begin(), pComp->compVars.end(), link) - 
        pComp->compVars.begin();
    //position of linking variable in current component
    int linkPos = find(compVars.begin(), compVars.end(), link) - 
        compVars.begin();
        
        
    if(dLinkPos >= pComp->compVars.size()) {
      mlog1(sbuilder() << "dLinkPos: " << dLinkPos);
      mlog("compVars", compVars);
      throw std::invalid_argument("depend link wasn't in compVar");
    }
    if(linkPos >= compVars.size()) {
      throw std::invalid_argument("link wasn't in compVars");
    }
    mlog1(sbuilder() << "dLinkPos: " << dLinkPos);
    //mlog1(sbuilder() << "linkPos: " << linkPos);
    
    
    //mlog("vars", varIndexes);
    //mlog("dvars", pComp->varIndexes);
    
    mlog1(pComp->timeStepPipe.tms.size());
    
    TaylorModel transformedLink = pComp->timeStepPipe.tms
        .at(dLinkPos).transform(compMappers.at(depIndex));
    mlog1("between");
    timeStepPipe.tms.at(linkPos) = transformedLink;
    mlog("remapped", transformedLink);
  }
  
  //mlog("init", initSet);
  mdec();
  mlog1("remapping last >");
  mrestore(old);
}

void MyComponent::prepare(TaylorModelVec tmv, const vector<HornerForm> & ode, 
      vector<Interval> domain) {
  prepareVariables(tmv, ode);
  
  //add flowpipes if the component is not initial one
  //variables in component + link vars
  int varSize = varIndexes.size() + dependencies.size();
  
  if(dependencies.size() != 0)
    for(int i = 0; i < dependencies.at(0)->pComp->pipes.size(); i++) {
      addEmptyTM(pipes, varSize);
    }
  
  remapFlowpipes();
  mlog("all tm params", allTMParams);
  
  mlog("to solve variables", varIndexes);
  mlog("component variables", compVars);
  mlog("linking variables", linkVars);
  mlog("solve indexes", solveIndexes);
  mlog1(sbuilder() << "# of pipes: " << pipes.size());
  remapIVP(tmv, ode, domain);
}

vector<MyComponent *> createComponents(vector< vector<int> > compIndexes, 
    const vector<HornerForm> & ode) {
  mreset(old);
  mdisable();
  mlog1("creating <");
  minc();
  
  //number of variable is equal to number of ODEs
  int varSize = ode.size();
  
  vector<MyComponent *> components;
  MyComponent *lookup[varSize];
  
  //add variables to components
  //create a lookup from variables to components
  for(vector< vector <int> >::iterator cit = compIndexes.begin(); 
      cit < compIndexes.end(); cit++) {
    mlog("componenet vars", *cit);
    MyComponent *c = new MyComponent();
    for(vector<int>::iterator it = cit->begin(); it < cit->end(); it++) {
      c->addVar(*it);
      lookup[*it] = c;
    }
    components.push_back(c);
    
  }
  
  int vars[varSize];
  for(vector<MyComponent *>::iterator cit = components.begin(); 
      cit < components.end(); cit++) {
    
    //initial vars with 0s
    memset(vars, 0, varSize*sizeof(int) );
    
    //find all the variables that the component depends on
    for(vector<int>::iterator it = (*cit)->varIndexes.begin(); 
        it < (*cit)->varIndexes.end(); it++) {
      mlog1(sbuilder() << "var: " << *it);
      mlog1(ode.at(*it).toString());
      
      //populate vars with variables that exist in the ode
      bool b = ode.at(*it).getVars(vars);
    }
    
    //nullify variables that are going to be solved in component
    for(vector<int>::iterator it = (*cit)->varIndexes.begin(); 
        it < (*cit)->varIndexes.end(); it++) {
      vars[*it] = 0;
    }
    
    //add component dependency for variables that are not in the current
    //component, but who still exist in the ode.
    for(int i = 0; i < varSize; i++) {
      if(vars[i] > 0) {
        mlog1(sbuilder() << "i: " << i);
        (*cit)->addDependency(i, lookup[i]);
      }
    }
    vector<int> v(vars, vars + sizeof vars / sizeof vars[0]);
    mlog("v", v);
  }
  mdec();
  mlog1("creating >");
  mrestore(old);
  return components;
}



MyComponent getSystemComponent(vector<MyComponent *> comps,
    TaylorModelVec init, const vector<HornerForm> & ode,
    vector<Interval> domain) {
  MyComponent allVars;
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    for(vector<int>::iterator i2 = (*it)->varIndexes.begin();
        i2 < (*it)->varIndexes.end(); i2++) {
      allVars.addDependency(*i2, *it);
    }
  }
  allVars.prepareComponent(init, ode, domain);
  allVars.unpairedRight = MyComponentRemove::getUnitTmv(init.tms.size());
  return allVars;
}


void prepareComponents(vector<MyComponent *> & comps, TaylorModelVec init, 
    const vector<HornerForm> & ode, vector<Interval> domain) {
  mreset(old);
  mdisable();
  mlog1("preparing components <");
  minc();
  mforce("shouldn't be used");
  exit(0);
  
  mlog("init", init);
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    (*it)->prepareComponent(init, ode, domain);
    mlog("pre Init", (*it)->initSet);
  }

  mdec();
  mlog1("preparing components >");
  mrestore(old);
}

int MyComponent::nextFreeParam = -1;



PrecondModel *MyComponent::lastPre() {
  return pipePairs[pipePairs.size() - 1];
}

TaylorModelVec MyComponent::lastPipe() {
  throw invalid_argument("don't use");
  return pipes[pipes.size() - 1];
}

bool MyComponent::isSolveVar(int var) {
  vector<int>::iterator it = find(solveIndexes.begin(), solveIndexes.end(), var);
  return it != solveIndexes.end();
}

bool MyComponent::belongsToComp(int param) {
  vector<int>::iterator it = find(tpIndexes.begin(), tpIndexes.end(), param);
  return it != tpIndexes.end();
}

void MyComponent::serializeFlows() {
  mlog1("writing flows");
  //exit(0);
  
  mreset(old);
  mlog1("bar");
  FILE *fpDumping = fopen("temp.txt", "w");
  mlog("last", lastPipe().tms[0]);
  vector<string> names;
  names.push_back("t");
  names.push_back("a1");
  names.push_back("a2");
  
  vector<string> varNames;
  varNames.push_back("x1");
  varNames.push_back("x2");
  
  
  parseSetting.clear();
	parseSetting.addVar("t");
	parseSetting.addVar("a1");
	parseSetting.addVar("a2");
  
  
  fprintf(fpDumping, "flowpipes{\n");
  bool first = true;
  int n = 1;
  for(vector<TaylorModelVec>::iterator pipeIt = pipes.begin();
        pipeIt < pipes.end(); pipeIt++) {
    if(first == false)
      fprintf(fpDumping, ",\n");
    first = false;
    
    //mlog("pip", *pipeIt);
    pipeIt->serialize(fpDumping, names);
    //parseTMV(sbuilder() << "my models{" << n++ << " * a1 + [-0.001,0.001]}").serialize(fpDumping, varNames, names);
    //break;
  }
  fprintf(fpDumping, "}\n");

  fclose(fpDumping);
}


TaylorModel MyComponent::remapSingleInitTM(int variable) {
  mreset(old);
  mdisable();
  minc();
    
  mlog1("remapping single left");
  mlog1(sbuilder() << "variable: " << variable);
  
  int depPos, dLinkPos, linkPos;
  computeMappingPositions(variable, &depPos, &dLinkPos, &linkPos);
  
  MyComponent *pComp = dependencies[depPos]->pComp;
  
  mlog("pInit", pComp->initSet);
  mlog1(pComp->initSet.getParamCount());
  mlog1(initSet.getParamCount());
  mlog("ma", compMappers[depPos]);
  
  //actual mapping of flowpipes
  TaylorModel tm = pComp->initSet.tms.at(dLinkPos).transform(compMappers[depPos]);
  
  mlog("tr", tm);
  mlog1(tm.getParamCount());
  
  mdec();
  mrestore(old);
  
  return tm;
}

void MyComponent::computeMappingPositions(int variable, int *depPos, 
    int *dLinkPos, int *linkPos) {
  mlog1("computing positions");
  mlog1(sbuilder() << "variable: " << variable);
  
  //dependency position
  *depPos = MyComponentRemove::findPos(variable, &linkVars);
  
  if(*depPos >= linkVars.size()) {
    throw invalid_argument(sbuilder() << "not wasn't linkVars, pos: " << *depPos); 
  }
  
  MyComponent *pComp = dependencies[*depPos]->pComp;
  mlog1(sbuilder() << "*depPos: " << *depPos);
  
  //position of linking variable in precomputed component
  *dLinkPos = MyComponentRemove::findPos(variable, &(pComp->compVars));
  //position of linking variable in current component
  *linkPos = MyComponentRemove::findPos(variable, &compVars);
  if(*dLinkPos >= pComp->compVars.size()) {
    mlog1(sbuilder() << "*dLinkPos: " << *dLinkPos);
    mlog("compVars", compVars);
    throw std::invalid_argument("depend link wasn't in compVar");
  }
  if(*linkPos >= compVars.size()) {
    throw std::invalid_argument("link wasn't in compVars");
  }
  mlog1(sbuilder() << "*dLinkPos: " << *dLinkPos);
  mlog1(sbuilder() << "*linkPos: " << *linkPos);
  
  mlog("compVars", compVars);
  mlog("dvars", pComp->varIndexes); 
  
}

TaylorModel MyComponent::remapSingleRightTM(int variable) {
  mreset(old);
  mdisable();
  minc();
  
  mlog1("remapping single right");
  
  int depPos, dLinkPos, linkPos;
  computeMappingPositions(variable, &depPos, &dLinkPos, &linkPos);
  MyComponent *pComp = dependencies[depPos]->pComp;
  
  mlog("pInit", pComp->initSet);
  mlog1(pComp->initSet.getParamCount());
  mlog1(initSet.getParamCount());
  mlog("ma", compMappers[depPos]);
  
  //actual mapping of flowpipes
  TaylorModel tm = pComp->unpairedRight.tms.at(dLinkPos)
      .transform(compMappers[depPos]);
  
  mlog("tm", tm);
  mlog1(tm.getParamCount());
  
  
  mdec();
  mrestore(old);
  return tm;
}

void MyComponent::getIthPipePair(vector<int> lMapper, vector<int> rMapper, 
        TaylorModel & left, TaylorModel & right, int var, int i) {
  mforce("getting ith");
  
  int pos = MyComponentRemove::findPos(var, &varIndexes);
  if(pos < varIndexes.size()) {
    //TaylorModel left = pipePairs[i]->left.tms[pos];
    
    //native to this component
    if(lMapper.size() == 0) {
      //return tm;
    }
    
    //TaylorModel mapped = tm.transform(mapper);
    
    //mlog("tm", tm);
    //mlog("ma", mapped);
    //mlog1(tm.getParamCount());
    //mlog1(mapped.getParamCount());
  }
  
  
}



TaylorModel MyComponent::getRightModelForVar(vector<int> mapper, int var) {
  mreset(old);
  mdisable();
  mlog1(sbuilder() << "getting flowpipe for " << var);
  mlog("mapper", mapper);
  
  int pos = MyComponentRemove::findPos(var, &varIndexes);
  if(pos < varIndexes.size()) {
    mlog1(mapper.size());
    mlog1("can return straight");
  
    TaylorModel tm = unpairedRight.tms[pos];
    
    //native to this component
    if(mapper.size() == 0) {
      mrestore(old);
      return tm;
    }
    
    TaylorModel mapped = tm.transform(mapper);
    
    //mlog("tm", tm);
    //mlog("ma", mapped);
    //mlog1(tm.getParamCount());
    //mlog1(mapped.getParamCount());
    mrestore(old);
    return mapped;
  }
  
  for(int j = 0; j < dependencies.size(); j++) {
    CompDependency *dep = dependencies[j];
    mlog("depAll", dep->pComp->allVars);
    
    //is variable in dependency implicit variables
    bool in = MyComponentRemove::isIn(var, &dep->pComp->allVars);
    
    //if variable is not in dependency skip dependency
    if(in == false) {
      continue;
    }
    mlog("mapper", mapper);
    mlog("dmapper", dep->rightMapper);
    
    vector<int> concatenated = concateMapper(dep->rightMapper, mapper);
    
    mlog("con", concatenated);
    mrestore(old);
    return dep->pComp->getRightModelForVar(concatenated, var);
  }
  
  
  throw std::runtime_error("should never get to this point, either variable is in componenent or one of the dependencies");
}









