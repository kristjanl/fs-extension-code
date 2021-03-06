#!/usr/bin/env python

import sys
import os
import re

# directory of models
#modelDir = os.path.join("..", "models", "composition")
modelDir = os.path.join("..", "models", "compositional", "id_comp_precond")

# directory of the outputs
#outDir = os.path.join("..", "models", "compositional", "temp")
outDir = os.path.join("..", "models", "compositional", "qr")


def change_init(m, f, changer, f_lower, f_upper):
  lower = float(m.group(4))
  upper = float(m.group(5))
  f.write("%s%s%s[%s,%s]\n"%(
      m.group(1),m.group(2), m.group(3), 
      f_lower(m.group(2), lower), f_upper(m.group(2), upper)))

def fixed_value(m, f, changer, value):
  f.write("%s%s\n" %(m.group(1), value))

def change_method_to_id(m, f, changer):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), "identity precondition"))
  
def change_method(m, f, changer, value):
  #use previous level of tabbing
  #print "changing method"
  f.write("%s%s\n" %(m.group(1), value))
  
def change_output(m, f, changer):
  #use previous level of tabbing and file prefix
  f.write("%s%s%s\n" %(m.group(1), m.group(2), changer.suffix))
  
def change_output2(m, f, changer):
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
  
def change_modelname2(name, suffix):
  #matches single identifier
  m = re.search('(.*)_(.*)_(.*)\\.', name)
  return "%s_%s.model" %(m.group(1), suffix)
  
def dummy(m, f, changer):
  print "dummy"
  print m.group(0)
  print changer
  raise Exception('Dummy function needs to be swapped out')
  

class BaseChanger:
  #captures parameters for method, outputname, implementation, composition
  methodPat = '(.*?)((?:no precondition|shrink wrapping|identity precondition))'
  outputPat = '(.*?output .*_)(.*_)(.*)'
  implPat = '.*?alg_small_comp flow impl'
  decompPat = '(.*?)(decomposition .*|no decomposition)'
  stepPat = '(.*?fixed steps )(.*)'
  timePat = '(.*?time )(.*)'
  initPat = '(.*?)(x\d*)( in )\[(.*?),(.*?)\]'
  
  lookup = {
    methodPat:dummy, 
    outputPat:dummy, 
    implPat:dummy, 
    decompPat:dummy, 
    stepPat:same_line, 
    timePat:same_line, 
    initPat:same_line,
  }
  def __init__(self, suffix):
    self.suffix = suffix
    
  def changeLine(self, line, outFile):
    #print line
    for key, changerF in self.lookup.iteritems():
      #print key
      m = re.search(key, line)
      if m != None:
        changerF(m, outFile, self)
        return
    outFile.write(line)

class SameChanger(BaseChanger):
  def __init__(self):
    self.suffix = None
    self.lookup[self.methodPat] = same_line
    self.lookup[self.outputPat] = same_line
    self.lookup[self.implPat] = same_line
    self.lookup[self.decompPat] = same_line
  

class IdFlowChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_flow"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.implPat] = no_line
    self.lookup[self.decompPat] = no_line
    
class IdSameImplCompChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_comp"
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
  
  #changer.lookup[changer.stepPat] = (lambda m, f, c: fixed_value(m, f, c, 1.5))
  #changer.lookup[changer.timePat] = (lambda m, f, c: fixed_value(m, f, c, 60))
  #changer.lookup[changer.initPat] = lambda m, f, c: change_init(m, f, c, 
  #    lambda x, v: 0.9, #lower limit function
  #    lambda x, v: 1) #upper limit function
  
  changer = SameChanger()
  changer.suffix = "qr_flow"
  changer.lookup[changer.outputPat] = change_output2
  changer.lookup[changer.methodPat] = (lambda m, f, c: change_method(m, f, c, "QR precondition"))
  
  models = []
  dirpath = ""
  
  # filter out some of the files
  for (dirpath, dirnames, filenames) in os.walk(modelDir):
    for (i, file) in enumerate(filenames):
      rs = []
      rs.append(re.search('(?:_2_|_\d0_)', file))
      #rs.append(re.search('(?:_2_)', file))
      rs.append(re.search('_id', file))
      rs.append(re.search('_flow', file))
      #rs.append(re.search('_flow', file))
      #rs.append(re.search('lin_dep_\d', file))
      
      if reduce(lambda l, r: l and r, map(lambda r: r != None, rs)):
        models.append(file)
        
  
  #modify the remaining models with changer object
  for model in models:
    print model
    #continue
    
    inFile = os.path.join(dirpath, model)
    
    #change_modelname changes last one part, change_modelname2 last 2
    outModel = change_modelname2(model, changer.suffix)
    outName = os.path.join(outDir, outModel)
    outFile = open(outName, 'w')

    with open(inFile) as f:
      for line in f:
        changer.changeLine(line, outFile)
    outFile.close()
foo()  

   
