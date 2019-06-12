#!/usr/bin/env python

import sys
import os
import re

# directory of models
modelDir = os.path.join(".", "comp_base")

# directory of the outputs
outDir = os.path.join(".", "comp")


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

def change_method_to_qr(m, f, changer):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), "QR precondition"))

def change_method_to_tqr(m, f, changer):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), "TQR precondition"))

def change_method_to_pa(m, f, changer):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), "parallelepiped precondition"))
  
def change_method(m, f, changer, value):
  #use previous level of tabbing
  f.write("%s%s\n" %(m.group(1), value))
  
def change_output(m, f, changer):
  #use previous level of tabbing and file prefix
  f.write("%s%s%s\n" %(m.group(1), m.group(2), changer.suffix))
  
def change_output2(m, f, changer):
  #use previous level of tabbing and file prefix
  #print m.group(0)
  #print("%s%s%s\n" %(m.group(1), m.group(2), changer.suffix))
  f.write("%s%s_%s\n" %(m.group(1), m.group(2), changer.suffix))

def change_output_strip(m, f, changer):
  #use previous level of tabbing and file prefix
  #print m.group(0)
  #print("%s%s%s\n" %(m.group(1), m.group(2), changer.suffix))
  output = stripName(m.group(2))
  f.write("%s%s_%s\n" %(m.group(1), output, changer.suffix))
  
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
  print name
  m = re.search('(.*)_(.*)_(.*)\\.', name)
  return "%s_%s.model" %(m.group(1), suffix)
  
def change_modelname3(name, suffix):
  #matches single identifier
  m = re.search('(.*)\.model', name)
  return "%s_%s.model" %(m.group(1), suffix)

def stripName(s):
  s = s.strip()
  parts = s.split("_")
  discardedValues = ["id"]
  parts = filter(lambda s: not s in discardedValues, parts)
  ret = "_".join(parts)
  return ret

def change_modelname_strip(name, suffix):
  m = re.search('(.*)\.model', name)
  return "%s_%s.model" %(stripName(m.group(1)), suffix)
  
def dummy(m, f, changer):
  print "dummy"
  print m.group(0)
  print changer
  raise Exception('Dummy function needs to be swapped out')
  

class BaseChanger:
  
  def __init__(self, suffix):
    self.suffix = suffix

    #captures parameters for method, outputname, implementation, composition
    self.methodPat = '(\s*)((?:shrink wrapping|.*?precondition.*))'
    self.outputPat = '(.*?output )(.*)'
    self.decompPat = '(.*?)(decomposition .*|no decomposition)'
    self.cflowPat = '(.*?)use cflow'
    self.compPat = '(.*?)(auto|no) components'
    self.precompTypePat = '(.*?)(left model|fully) compositional'
    self.stepPat = '(.*?fixed steps )(.*)'
    self.timePat = '(.*?time )(.*)'
    self.initPat = '(.*?)(x\d*)( in )\[(.*?),(.*?)\]'
    
    self.lookup = {
      self.methodPat:dummy, 
      self.outputPat:dummy, 
      self.decompPat:dummy,
      self.cflowPat:dummy,
      self.compPat:dummy,
      self.precompTypePat:dummy,
      self.stepPat:same_line, 
      self.timePat:same_line, 
      self.initPat:same_line,
    }
    
    
  def changeLine(self, line, outFile):
    #print line
    for key, changerF in self.lookup.iteritems():
      m = re.search(key, line)
      if m != None:
        """if "identity" in line:
          print line
          print m
          print changerF
          print dummy
        """  
        changerF(m, outFile, self)
        return
    outFile.write(line)

class SameChanger(BaseChanger):
  def __init__(self):
    raise Exception("Can't call")
    """
    self.suffix = None
    self.lookup = {}
    self.lookup[self.methodPat] = same_line
    self.lookup[self.outputPat] = same_line
    self.lookup[self.decompPat] = same_line
    """
  

class IdFlowChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_flow"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.decompPat] = no_line
    
class IdSameImplCompChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_comp"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.decompPat] = same_line
    
class IdSameImplNoCompChanger(BaseChanger):
  def __init__(self):
    self.suffix = "id_nocomp"
    self.lookup[self.methodPat] = change_method_to_id
    self.lookup[self.outputPat] = change_output
    self.lookup[self.decompPat] = change_no_comp

def foo():
  changerN = BaseChanger("nocomp")
  changerN.lookup[changerN.cflowPat] = same_line
  changerN.lookup[changerN.compPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "no components"))
  changerN.lookup[changerN.precompTypePat] = same_line
  changerN.lookup[changerN.methodPat] = same_line
  changerN.lookup[changerN.outputPat] = change_output2
  
  changerL = BaseChanger("lcomp")
  changerL.lookup[changerL.cflowPat] = same_line
  changerL.lookup[changerL.compPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "auto components"))
  changerL.lookup[changerL.precompTypePat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "left model compositional"))
  changerL.lookup[changerL.methodPat] = same_line
  changerL.lookup[changerL.outputPat] = change_output2

  changerF = BaseChanger("fcomp")
  changerF.lookup[changerF.cflowPat] = same_line
  changerF.lookup[changerF.compPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "auto components"))
  changerF.lookup[changerF.precompTypePat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "fully compositional"))
  changerF.lookup[changerF.methodPat] = same_line
  changerF.lookup[changerF.outputPat] = change_output2
  
  changers = [changerN, changerL, changerF]
  #changers = [changerFId]

  models = []
  dirpath = ""
  
  # filter out some of the files
  for (dirpath, dirnames, filenames) in os.walk(modelDir):
    for (i, file) in enumerate(filenames):
      rs = []
      #rs.append(re.search('moore', file))
      rs.append(re.search('\.model', file))
      #rs.append(re.search('and_or', file))
      if reduce(lambda l, r: l and r, map(lambda r: r != None, rs)):
        models.append(file)
  
  #modify the remaining models with changer object
  for model in models:
    print model
    #continue
    
    inFile = os.path.join(dirpath, model)
    
    for c in changers:
      #change_modelname changes last one part, change_modelname2 last 2
      #outModel = change_modelname3(model, c.suffix)
      outModel = change_modelname3(model, c.suffix)
      outName = os.path.join(outDir, outModel)
      outFile = open(outName, 'w')

      with open(inFile) as f:
        for line in f:
          c.changeLine(line, outFile)
      outFile.close()
foo()  

   
