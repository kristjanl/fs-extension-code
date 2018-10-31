#include "MyComponent.h"

//todo ask if this is ok (due to invalid use of forward declaration from
//MySettings class (in Cont.h)
//https://stackoverflow.com/questions/6988924/invalid-use-of-incomplete-type-forward-declaration
#include "Continuous.h"
#include "Utils.h"
#include "PreconditionedTMV.h"
#include "Transformer.h"

using namespace std;

MyComponent::MyComponent(vector<int> vs, vector<int> tps):varIndexes(vs),tpIndexes(tps) {
  isSolved = false;
  isPrepared = false;
  isPreconditioned = false;
  firstPrecondition = true;
  usingPreconditioning = false;
}
MyComponent::MyComponent() {
  isSolved = false;
  isPrepared = false;
  isPreconditioned = false;
  firstPrecondition = true;
  usingPreconditioning = false;
}

void MyComponent::log() {
  mreset(old);
  //mdisable();
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
      mlog1(sbuilder() << odeIndex++ << ": " << it->toString());
      //mforce1(sbuilder() << "size: " << it->hornerForms.size());
    }
  }
  mlog("compDomain", dom);
  mlog("compInit", initSet);
  for(int i = 0; i < initSet.tms.size(); i++) {
    mlog1(sbuilder() << "tm param size(" << i << "): "
        << initSet.tms.at(i).getIgnoringParamCount());
  }
  mlog1(sbuilder() << "init set size: " << initSet.tms.size());
  mlog("to be solved variables", varIndexes);
  mlog("comp variables", compVars);
  mlog("all variables", allVars);
  mlog("solve indexes", solveIndexes);
  mlog("new params", tpIndexes);
  mlog("all params", allTMParams);
  mdec();
  mlog1("component >");
  mrestore(old);
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
    MyComponent *pComp = dependencies[i]->pComp;
    //mlog1(sbuilder() << "link: " << dependencies.at(i)->linkVar);
    //previous component parameters
    vector<int> prevCompPm = pComp->allTMParams;
    mlog("pm", prevCompPm);
    vector<int> mapper;
    //iterate all the variables in all the components
    for(int j = 0; j < allTMParams.size(); j++) {
      int var = allTMParams[j];
      
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
    dependencies[i]->mapper = mapper;
    dependencies[i]->rightMapper = mapper;
    mappers.push_back(mapper);
  }
  
  //mapper for left TaylorModel
  for(int i = 0; i < dependencies.size(); i++) {
    MyComponent *pComp = dependencies[i]->pComp;
    
    //variables in previous component
    vector<int> prevCompVars = pComp->allVars;
    //mlog("pv", prevCompVars);
    vector<int> mapper;
    
    //loop over all current component variables
    for(int j = 0; j < allVars.size(); j++) {
      int var = allVars[j];
      
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
    dependencies[i]->leftMapper = mapper;
  }
  
  for(int i = 0; i < dependencies.size(); i++) {
    mlog1(sbuilder() << dependencies[i]->linkVar);
    mlog("left", dependencies[i]->leftMapper);
    mlog("right", dependencies[i]->rightMapper);
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
  mlog("compVars", compVars);
  mlog("link variables", linkVars);
  mlog("tmv", tmv);
  
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
  
  mlog("compInit", initSet);
  
  mdec();
  mlog1("remapping IV >");
  mrestore(old);
}
  

//prepares all the dependent components
//prepares variables, mappers and intial set
//TODO could refactor to use IVP
void MyComponent::prepareComponent(TaylorModelVec init, 
    const vector<HornerForm> & ode, MySettings *settings) {
  //return if variables have already been prepared
  if(isPrepared) {
    return;
  }
  mreset(old);
  mdisable();
  mlog1("preparing component <");

  minc();
  mlog("vars", varIndexes);
  mlog1(sbuilder() << "dependencies: " << dependencies.size());
  //prepare variables for all the dependencies
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
        it < dependencies.end(); it++) {
    mlog1(sbuilder() << "link: " << (*it)->linkVar);
    MyComponent *pComp = (*it)->pComp;
    pComp->prepareComponent(init, ode, settings);
  }
  prepareVariables(init, ode, settings->discardEmptyParams);
  prepareMappers();
  remapIVP(init, ode, settings->domain);
  timeStepPipe = TaylorModelVec(this->initSet);
  
  //only create a part of right taylor model 
  //unpairedRight = getUnitTmv(allTMParams.size());
  unpairedRight = getUnitTmv(tpIndexes.size());
  unpairedRight = getNVarMParam(tpIndexes.size(),
      allTMParams.size());
  //mforce3(aaaa, "in prepare2", unpairedRight);


  //denote that using preconditioned transformer
  usingPreconditioning = settings->transformer->isPreconditioned;

  for(int j = 0; j < dependencies.size(); j++) {
    CompDependency dep = *dependencies[j];
    //use left mapper when preconditioned
    if(settings->transformer->isPreconditioned) {
      dep.mapper = dep.leftMapper;
    } else {
      dep.mapper = dep.rightMapper;
    }
  }

  isPrepared = true;
  
  mdec();
  mlog1("preparing component >");
  mrestore(old);
}
  
void MyComponent::prepareVariables(TaylorModelVec init, 
    const vector<HornerForm> & ode, bool discardEmptyParams) {
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
    
    //*don't retain anymore
    //retain the parameter if started with point initial conditions
    if(tmParams.size() == 0 && discardEmptyParams == false) {
      varsToBeIntroduced.push_back(*it);
      tmParams.push_back(*it);
    }
    //*/
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
  mlog("allParams (with dependencies)", allParams);
  
  mdec();
  mlog1("preparing variables >");
  mrestore(old);
}  

//can only be used in transform
TaylorModelVec MyComponent::orderedTSPRemap(bool first) {
  mreset(old);
  mdisable();
  mlog1("ordered tsp remapping <");
  minc();
  
  int varSize = varIndexes.size() + dependencies.size();

  //pipes.push_back(getEmptyTMV(varSize));
  timeStepPipe = getEmptyTMV(varSize);
  mlog1(dependencies.size());
  
  if(dependencies.size() == 0) {
    throw std::invalid_argument("orderedTSPRemap can only be used with deps");
  }
  
  TaylorModelVec ret(compVars.size());
  TaylorModelVec ret2;
  
  map<int, CompDependency *>  lookup;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++) {
    lookup[(*it)->linkVar] = (*it);
  }
  
  for(int i = 0; i < allVars.size(); i++) {
    int var = allVars[i];
    mlog1(sbuilder() << "var: " << var);
    CompDependency *dep = lookup[var];
    MyComponent *pComp = dep->pComp;
    mlog("tmv", pComp->timeStepPipe);
    mlog("vs", pComp->allVars);
    
    int dLinkPos = find(pComp->compVars.begin(), pComp->compVars.end(), var)
        - pComp->compVars.begin();
    TaylorModel & source = pComp->timeStepPipe.tms[dLinkPos];
    mlog("source", source);
    
    mlog("vs", pComp->allVars);
    
    mlog("l", dep->leftMapper);
    mlog("r", dep->rightMapper);
    
    TaylorModel mapped;
    
    //first one uses TM parameters, other ones use parameters for variable subs
    if(first) {
      mapped = source.transform(dep->rightMapper);
    } else {
      mapped = source.transform(dep->leftMapper);
    }
    mlog("mapped", mapped);
    ret2.tms.push_back(mapped);
    //mlog1("between");
  }
  mdec();
  mlog1("orderd tsp remapping >");
  mrestore(old);
  return ret2;
}



void MyComponent::remapTimeStepPipe() {
  mreset(old);
  mdisable();
  mlog1("remapping tsp <");
  minc();
  
  int varSize = varIndexes.size() + dependencies.size();

  //pipes.push_back(getEmptyTMV(varSize));
  timeStepPipe = getEmptyTMV(varSize);
  mlog1(sbuilder() << "dep size: " << dependencies.size());
  
  if(dependencies.size() == 0) {
    mdec();
    mlog1("remapping tsp >");
    mrestore(old);
    return;
  }
  
  //mapping of dependent component flowpipes to current component
  int depIndex = 0;
  for(vector<CompDependency *>::iterator it = dependencies.begin(); 
      it < dependencies.end(); it++, depIndex++) {
    mlog1(sbuilder() << "depIndex: " << depIndex);
    minc();
    int link = (*it)->linkVar;
    mlog1(sbuilder() << "link: " << link);
    MyComponent *pComp = (*it)->pComp;
    //mlog("dts", pComp->tpIndexes);
    //mlog("mapper", (*it)->mapper);
    //mlog("rightMapper", (*it)->rightMapper);
    mlog("leftMapper", (*it)->leftMapper);
    
    vector<int> leftMapper = (*it)->leftMapper;
    
    
    //mlog1("here");
    mlog("dav", pComp->allVars);
    
    
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
    
    mlog("stsp", pComp->timeStepPipe);
    
    mlog("source", pComp->timeStepPipe.tms[dLinkPos]);
    
    mlog1(sbuilder() << "source params: " << pComp->timeStepPipe.tms[dLinkPos].
        getIgnoringParamCount());
    
    TaylorModel transformedLink = pComp->timeStepPipe.tms[dLinkPos].
        transform(leftMapper);
    timeStepPipe.tms[linkPos] = transformedLink;
    mlog("mapped", transformedLink);
    
    mlog1(sbuilder() << "mapped params: " << transformedLink.
        getIgnoringParamCount());
        
    mdec();
  }
  
  //mlog("init", initSet);
  mdec();
  mlog1("remapping tsp >");
  mrestore(old);
}

#include <set>

//TODO remove after refactoring is complete
void makeCompIndexes(MySettings *settings, const vector<HornerForm> & ode) {
  mreset(old);
  mdisable();
  mlog1("making indexes <");
  minc();
  
  //number of variables
  int varSize = ode.size();
  mlog1(sbuilder() << "varSize: " << varSize);
  
  int vars[varSize];
  
  //array of influencers for each variable
  set<int> sets[varSize];
  
  for(int var = 0; var < varSize; var++) {
     mlog1(sbuilder() << "var: " << var);
     //clear
     memset(vars, 0, varSize*sizeof(int) );
     
     //array of whether a variable is present in the derivative or not
     //(signalled by >0 int)
     ode[var].getVars(vars);
     
     //add variable itself to its influencers
     sets[var].insert(var);
     for(int var2 = 0; var2 < varSize; var2++) {
        //add variable to influencers if it was present in the derivative
        if(vars[var2] != 0) {
          sets[var].insert(var2);
        }
     }
     mlog("set", sets[var]);
  }
  
  bool repeat = true;
  while(repeat) {
    //reset repeating
    repeat = false;
    for(int var = 0; var < varSize; var++) {
      //currrent influencer of variable
      set<int> & infls = sets[var];
      int oldSize = infls.size();
  
      //could optimize by only inserting the ones that changed, but likely this
      //shouldn't matter much
      //loop over all influencers
      for(set<int>::iterator it = infls.begin(); it != infls.end(); it++) {
        //add the influencers of the influencer
        infls.insert(sets[*it].begin(), sets[*it].end());
      }
      int newSize = infls.size();
      
      //if new influencer is found, then need to repeat
      if(oldSize != newSize) {
        repeat = true;
      }
    }
  }
  
  for(int var = 0; var < varSize; var++) {
    mlog(sbuilder() << "s" << var, sets[var]);
  }
  
  //variable already in some component
  set<int> seenVars;
  for(int var = 0; var < varSize; var++) {
    //skip variable if it is already in component
    if(seenVars.count(var) != 0)
      continue;
    set<int> & infls = sets[var];
    
    //variable in the same component as var
    vector<int> comp;
    for(set<int>::iterator it = infls.begin(); it != infls.end(); it++) {
      int infl = *it;
      //add variable to this component if var is influencer of infl
      //(var itself is always influencer of itself, so this is also added)
      if(sets[infl].count(var) != 0) {
        comp.push_back(infl);
        seenVars.insert(infl);
      }
    }
    //mlog("comp", comp);
    
    //add this component to settings
    settings->intComponents.push_back(comp);
  }
  
  mdec();
  mlog1("making indexes >");
  mrestore(old);
}

vector<MyComponent *> createComponents(MySettings *settings, 
    const vector<HornerForm> & ode) {
  mreset(old);
  mdisable();
  mlog1("creating <");
  minc();
  
  if(settings->autoComponents) {
    //settings->intComponents.clear(); //comment in if forcing auto comp
    makeCompIndexes(settings, ode);
  }
  mlog1(sbuilder() << "ci size: " << settings->intComponents.size());
  
  vector< vector<int> > compIndexes = settings->intComponents;
  //mlog("ci", compIndexes[0]);
  mlog1(sbuilder() << "ac:" << settings->autoComponents);
  
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
      ode.at(*it).getVars(vars);
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


MyComponent* pGetSystemComponent(vector<MyComponent *> comps,
    TaylorModelVec init, const vector<HornerForm> & ode,
    MySettings *settings) {
  mreset(old);
  mdisable();
  MyComponent *ret = new MyComponent;
  for(vector<MyComponent *>::iterator it = comps.begin(); 
      it < comps.end(); it++) {
    for(vector<int>::iterator i2 = (*it)->varIndexes.begin();
        i2 < (*it)->varIndexes.end(); i2++) {
      ret->addDependency(*i2, *it);
    }
  }
  ret->prepareComponent(init, ode, settings);
  
  //mlog("av", ret.allVars);
  //mlog("cv", ret.compVars);
  
  TaylorModelVec temp;  
  //need to reorder initset wrt allVars
  for(int i = 0; i < ret->allVars.size(); i++) {
    int var = ret->allVars[i];
    int place = findPos(var, &ret->compVars);
    mlog1(sbuilder() << "place: " << place);
    temp.tms.push_back(ret->initSet.tms[place]);
  }
  ret->initSet = temp;
  
  //ret->unpairedRight = getUnitTmv(init.tms.size());
  ret->unpairedRight = getNVarMParam(init.tms.size(), 
      ret->allTMParams);
  mrestore(old);
  
  ret->usingPreconditioning = comps[0]->usingPreconditioning;
  
  return ret;
}

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



void MyComponent::computeMappingPositions(int variable, int *depPos, 
    int *dLinkPos, int *linkPos) {
  mlog1("computing positions");
  mlog1(sbuilder() << "variable: " << variable);
  
  //dependency position
  *depPos = findPos(variable, &linkVars);
  
  if(*depPos >= linkVars.size()) {
    throw invalid_argument(sbuilder() << "not wasn't linkVars, pos: " << *depPos); 
  }
  
  MyComponent *pComp = dependencies[*depPos]->pComp;
  mlog1(sbuilder() << "*depPos: " << *depPos);
  
  //position of linking variable in precomputed component
  *dLinkPos = findPos(variable, &(pComp->compVars));
  //position of linking variable in current component
  *linkPos = findPos(variable, &compVars);
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

void MyComponent::getIthPipePair(vector<int> lMapper, vector<int> rMapper, 
        TaylorModel & left, TaylorModel & right, int var, int i) {
  mforce1("getting ith");
  
  int pos = findPos(var, &varIndexes);
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
  
  int pos = findPos(var, &varIndexes);
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
    bool in = isIn(var, &dep->pComp->allVars);
    
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



int MyComponent::getIntergrationParamCount() {
  if(usingPreconditioning) {
    return allVars.size() + 1;
  }
  return allTMParams.size() + 1;
}

string MyComponent::getVarName(MySettings *settings) {
  stringstream ss;
  vector<int> & v = varIndexes;
  vector<int>::iterator it = v.begin();
  if(it != v.end())
    ss << settings->varNames[*it++];
  for(; it != v.end(); it++) {
    ss << "," << settings->varNames[*it];
  }
  return ss.str();
}




