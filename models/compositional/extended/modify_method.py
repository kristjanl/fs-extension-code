#!/usr/bin/env python

import sys
import os
import re

# directory of models
modelDir = os.path.join(".", "method_base")

# directory of the outputs
outDir = os.path.join(".", "method")


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
  changerFQr = BaseChanger("flowqr")
  changerFQr.lookup[changerFQr.cflowPat] = no_line
  changerFQr.lookup[changerFQr.compPat] = no_line
  changerFQr.lookup[changerFQr.precompTypePat] = no_line
  changerFQr.lookup[changerFQr.outputPat] = change_output_strip
  changerFQr.lookup[changerFQr.methodPat] = change_method_to_qr
  
  changerFId = BaseChanger("flowid")
  changerFId.lookup[changerFId.cflowPat] = no_line
  changerFId.lookup[changerFId.compPat] = no_line
  changerFId.lookup[changerFId.precompTypePat] = no_line
  changerFId.lookup[changerFId.outputPat] = change_output_strip
  changerFId.lookup[changerFId.methodPat] = change_method_to_id

  changerQr = BaseChanger("qr")
  changerQr.lookup[changerQr.cflowPat] = same_line
  changerQr.lookup[changerQr.compPat] = same_line
  changerQr.lookup[changerQr.precompTypePat] = no_line
  changerQr.lookup[changerQr.outputPat] = change_output_strip
  changerQr.lookup[changerQr.methodPat] = change_method_to_qr
  

  changerId = BaseChanger("id")
  changerId.lookup[changerId.cflowPat] = same_line
  changerId.lookup[changerId.compPat] = same_line
  changerId.lookup[changerId.precompTypePat] = same_line
  changerId.lookup[changerId.outputPat] = change_output_strip
  changerId.lookup[changerId.methodPat] = change_method_to_id

  changerPa = BaseChanger("pa")
  changerPa.lookup[changerPa.cflowPat] = same_line
  changerPa.lookup[changerPa.compPat] = same_line
  changerPa.lookup[changerPa.precompTypePat] = same_line
  changerPa.lookup[changerPa.outputPat] = change_output_strip
  changerPa.lookup[changerPa.methodPat] = change_method_to_pa

  changerTqr = BaseChanger("tqr")
  changerTqr.lookup[changerTqr.cflowPat] = same_line
  changerTqr.lookup[changerTqr.compPat] = same_line
  changerTqr.lookup[changerTqr.precompTypePat] = no_line
  changerTqr.lookup[changerTqr.outputPat] = change_output_strip
  changerTqr.lookup[changerTqr.methodPat] = change_method_to_tqr

  changerNop = BaseChanger("nop")
  changerNop.lookup[changerNop.cflowPat] = same_line
  changerNop.lookup[changerNop.compPat] = same_line
  changerNop.lookup[changerNop.precompTypePat] = no_line
  changerNop.lookup[changerNop.outputPat] = change_output_strip
  changerNop.lookup[changerNop.methodPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "no processing"))

  changerSw1 = BaseChanger("sw1")
  changerSw1.lookup[changerSw1.cflowPat] = same_line
  changerSw1.lookup[changerSw1.compPat] = same_line
  changerSw1.lookup[changerSw1.precompTypePat] = no_line
  changerSw1.lookup[changerSw1.outputPat] = change_output_strip
  changerSw1.lookup[changerSw1.methodPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "shrink wrapping 1"))

  changerSw2 = BaseChanger("sw2")
  changerSw2.lookup[changerSw2.cflowPat] = same_line
  changerSw2.lookup[changerSw2.compPat] = same_line
  changerSw2.lookup[changerSw2.precompTypePat] = no_line
  changerSw2.lookup[changerSw2.outputPat] = change_output_strip
  changerSw2.lookup[changerSw2.methodPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "shrink wrapping 2"))

  changerSw5 = BaseChanger("sw5")
  changerSw5.lookup[changerSw5.cflowPat] = same_line
  changerSw5.lookup[changerSw5.compPat] = same_line
  changerSw5.lookup[changerSw5.precompTypePat] = no_line
  changerSw5.lookup[changerSw5.outputPat] = change_output_strip
  changerSw5.lookup[changerSw5.methodPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "shrink wrapping 5"))

  changerSw10 = BaseChanger("sw10")
  changerSw10.lookup[changerSw10.cflowPat] = same_line
  changerSw10.lookup[changerSw10.compPat] = same_line
  changerSw10.lookup[changerSw10.precompTypePat] = no_line
  changerSw10.lookup[changerSw10.outputPat] = change_output_strip
  changerSw10.lookup[changerSw5.methodPat] = lambda m, f, c: f.write("%s%s\n" %(m.group(1), "shrink wrapping 10"))
  
  changers = [changerFQr, changerFId, changerQr, changerId, changerPa, changerTqr, changerNop, changerSw1, \
		changerSw2, changerSw5, changerSw10]
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
      outModel = change_modelname_strip(model, c.suffix)
      outName = os.path.join(outDir, outModel)
      outFile = open(outName, 'w')

      with open(inFile) as f:
        for line in f:
          c.changeLine(line, outFile)
      outFile.close()
foo()  

   
