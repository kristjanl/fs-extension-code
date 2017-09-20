#!/usr/bin/env python

import sys
import os
import re

# directory of models
modelDir = os.path.join("..", "models", "composition")

# directory of the outputs
outDir = os.path.join("..", "models", "compositional", "id")



def change_method_to_id(m, f, changer):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), "identity precondition"))
def change_output(m, f, changer):
  #use previous level of tabbing and file prefix
  f.write("%s%s\n" %(m.group(1), changer.suffix))
  
def same_line(m, f, changer):
  f.write("%s\n" %m.group(0))
  
def change_no_comp(m, f, changer):
  #use previous level of tabbing
  f.write("%sno decomposition\n" %m.group(1))
    
def no_line(m, f, changer):
  #print "no line")
  pass
  
  
  
def change_modelname(name, suffix):
  #matches single identifier
  m = re.search('(.*)_(.*)\\.', name)
  return "%s_%s.model" %(m.group(1), suffix)
  
def dummy(m, f, changer):
  print "dummy"
  raise Exception('Dummy function needs to be swapped out')
  

class BaseChanger:
  #captures parameters for method, outputname, implementation, composition
  methodPat = '(.*?)((?:no precondition|shrink wrapping))'
  outputPat = '(.*?output .*_)(.*)'
  implPat = '.*?alg_small_comp flow impl'
  decompPat = '(.*?)(decomposition .*|no decomposition)'
  
  lookup = {
    methodPat:dummy, 
    outputPat:dummy, 
    implPat:dummy, 
    decompPat:dummy,
  }
  def __init__(self, suffix):
    self.suffix = suffix
    
  def changeLine(self, line, outFile):
    #print line
    for key, changerF in self.lookup.iteritems():
      m = re.search(key, line)
      if m != None:
        changerF(m, outFile, self)
        return
    outFile.write(line)

class IdFlowChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_flow"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.implPat] = no_line
    self.lookup[self.decompPat] = no_line
    
class IdSameImplCompChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_comp2"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.implPat] = same_line
    self.lookup[self.decompPat] = same_line
    
class IdSameImplNoCompChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_nocomp"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.implPat] = same_line
    self.lookup[self.decompPat] = change_no_comp

def foo():
  
  changer = IdFlowChanger()
  changer = IdSameImplCompChanger()
  #changer = IdSameImplNoCompChanger()

  models = []
  dirpath = ""
  
  # filter out some of the files
  for (dirpath, dirnames, filenames) in os.walk(modelDir):
    for (i, file) in enumerate(filenames):
      r1 = re.search('(?:_2_|_\d0_)', file)
      r2 = re.search('_comp', file)
      if r1 != None and r2 != None:
        models.append(file)
        
        
  #modify the remaining models with changer object
  for model in models:
    inFile = os.path.join(dirpath, model)
    
    outModel = change_modelname(model, changer.suffix)
    outName = os.path.join(outDir, outModel)
    outFile = open(outName, 'w')
    print outName
    
    with open(inFile) as f:
      for line in f:
        changer.changeLine(line, outFile)
    outFile.close()
foo()  

   
