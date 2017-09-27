#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)
import my_functions as fs
import table_functions as tf
import math

# directory of models
modelDir = os.path.join("..", "models", "composition")
flowstar = os.path.join("..", "src", "flowstar")



regexp = re.compile(r'lin_[0-9]')
#regexp = re.compile(r'lin_dep_[0-9]')
#regexp = re.compile(r'sq_deg_[0-9]')
#regexp = re.compile(r'pair_dep_[0-9]')
comp = "_nocomp"


def getData(model):
  modelFile = os.path.join(modelDir, model)
  outputName = fs.getParam(modelFile, "output")
  infoFile = "infos/%s.txt" %outputName
  compTime = fs.getParam(infoFile, "computation time:")
  dim = tf.getDimension(modelFile)
  return compTime

def foo():
  prefix = "lin"
  #prefix = "lin_dep"
  prefix = "sq_deg"
  prefix = "sq_deg_long"
  #prefix = "pair_dep"
  fmaker = lambda prefix, dim, composition: "%s_%s_%s.model" %(prefix, dim, composition)
  print prefix
  l = []
  for i in range(1, 51, 1):
    cfile = fmaker(prefix, i, "comp")
    nfile = fmaker(prefix, i, "nocomp")
    if not os.path.isfile(os.path.join(modelDir, cfile)):
      continue
    #print cfile
    #print nfile
    print "%s,%s,%s" %(i, getData(cfile), getData(nfile))
    l += [(i, getData(cfile), getData(nfile))]
    
  
  cfile = fmaker(prefix, 100, "comp")
  nfile = fmaker(prefix, 100, "nocomp")
  print "%s,%s,%s" %(100, getData(cfile), getData(nfile))
  l += [(100, getData(cfile), getData(nfile))]
  rows = int(math.ceil(len(l) / 3.0))
  
  for i in range(rows):
    print "%s & %.2f & %.2f &" %(l[i][0],float(l[i][1]), float(l[i][2])),
    print "%s & %.2f & %.2f &" %(l[i + rows][0],float(l[i + rows][1]), float(l[i + rows][2])),
    if i+2*rows < len(l):
      print "%s & %.2f & %.2f" %(l[i+2*rows][0],float(l[i+2*rows][1]), float(l[i+2*rows][2])),
    else:
      print "& &",
    print "\\\\"
    print "\hline"
  
foo()
