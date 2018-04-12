#!/usr/bin/env python

import sys
import os
import re

# directory of models
modelDir = os.path.join(".", "base")

# directory of the outputs
outDir = os.path.join(".", "stage")

# filter out some of the files

def replaceVars(line, vs, dim, outFile):
  
  for i in range(dim):
    line2 = line
    for v in vs:
      line2 = line2.replace(v, v + str(i))
    #print line2
    outFile.write(line2)
  
def justPrint(line):
  if re.search('init|{|}', line):
    return True
    
  if len(line.strip()) == 0:
    return True
  return False  

def makeVars(m, dim, outFile):
  vs = re.split(',', m.group(2))
  outFile.write(m.group(1))
  strVars = ""
  for i in range(dim):
    for v in vs:
      strVars = strVars + "%s,"%(v + str(i))
  outFile.write(strVars[:-1])
  outFile.write('\n')
  return vs
  

def makeComposition(vs, dim, outFile):
  s = "["
  for i in range(dim):
    s = s + "["
    for v in vs:
      s += v + str(i) + ","
    s = s[:-1] + "],"
  s = s[:-1] + "]\n"
  outFile.write('  decomposition ')
  outFile.write(s)

def foo(dim):
  models = []
  for (dirpath, dirnames, filenames) in os.walk(modelDir):
    for (i, file) in enumerate(filenames):
      rs = []
      #rs.append(re.search('moore', file))
      rs.append(re.search('\.model', file))
      if reduce(lambda l, r: l and r, map(lambda r: r != None, rs)):
        models.append(file)
        
  
  #modify the remaining models with changer object
  for model in models:
    print model
    inFile = os.path.join(dirpath, model)
    outModelName = model.replace("nodim", str(dim))
    outModelName = outModelName.replace(".model", "_comp.model")
    print outModelName
    outName = os.path.join(outDir, outModelName)
    outFile = open(outName, 'w')

    place = 0
    with open(inFile) as f:
      for line in f:
        m1 = re.search('(\s*state var )(.*)', line)
        if m1:
          vs = makeVars(m1, dim, outFile)
          continue
          
        if re.search('gnuplot octagon', line):
          outFile.write('  gnuplot octagon %s0,%s0\n'%(vs[0],vs[0]))
          continue
          
        
        m2 = re.search('(\s*output .*)_nodim_(.*)', line)
        if m2:
          outFile.write("%s_%s_%s_comp\n"%(m2.group(1), str(dim), m2.group(2)))
          continue
          
        if re.search('no decomposition', line):
          makeComposition(vs, dim, outFile)
          continue
          
        if place == 0:
          outFile.write(line)
          
        
          
        if re.search('poly ode 1', line):
          place = 1
          continue
          
        if place == 1 and justPrint(line):
          outFile.write(line)
          continue
          
        if place == 1:
          replaceVars(line, vs, dim, outFile)
    outFile.close()
    
    outModelName2 = outModelName.replace("_comp.", "_nocomp.")
    outName2 = os.path.join(outDir, outModelName2)
    outFile2 = open(outName2, 'w')
    
    with open(outName) as f:
      for line in f:
        if re.search('(\s*output .*)', line):
          outFile2.write(line.replace("comp", "nocomp"))
          continue
        if re.search('(\s*decomposition)', line):
          outFile2.write("  no decomposition\n")
          continue
        outFile2.write(line)
    outFile2.close()
    
    
foo(1)
