#!/usr/bin/env python

import sys
import os
import re

# directory of models
modelDir = os.path.join(".")

# directory of the outputs
outDir = os.path.join(".", "upd")


def foo():
  models = []
  dirpath = ""
  # filter out some of the files
  for (dirpath, dirnames, filenames) in os.walk(modelDir):
    for (i, file) in enumerate(filenames):
      if dirpath != ".":
        continue
      rs = []
      #rs.append(re.search('moore', file))
      rs.append(re.search('\.model', file))
      #rs.append(re.search('and_or', file))
      if reduce(lambda l, r: l and r, map(lambda r: r != None, rs)):
        models.append(file)
  for model in models:
    print model
    inFile = os.path.join(modelDir, model)
    
    outModel = "a"
    
    temp = model.split(".")
    outModel = ".".join(temp[:-2] + [temp[-2]+"_id"] + [temp[-1]])

    outName = os.path.join(outDir, outModel)
    #print outName
    outFile = open(outName, 'w')

    with open(inFile) as f:
      matches = []
      for line in f:
        m1 = re.search("(.*?)fixed steps", line)
        if m1 != None:
          matches.append(1)
          outFile.write(m1.group(1) + "use cflow\n")
          outFile.write(m1.group(1) + "no components\n")
          outFile.write(line)
          continue
        m2 = re.search("(.*?output .*)", line)
        if m2 != None:
          matches.append(2)
          outFile.write(m2.group(1) + "_id\n")
          continue 
        m3 = re.search("alg_small_comp flow impl", line)
        if m3 != None:
          matches.append(3)
          #outFile.write(line)
          continue
        m4 = re.search("decomposition \[\[|no decomposition", line)
        if m4 != None:
          matches.append(4)
          #outFile.write(line)
          continue
        m5 = re.search("(.*?)QR precondition", line)
        if m5 != None:
          matches.append(5)
          outFile.write(m5.group(1) + "left model compositional\n")
          outFile.write(m5.group(1) + "identity precondition\n")
          continue
        outFile.write(line)
      if len(matches) != 5:
        raise Exception("model '%s' didn't have all matches (%s)" %(model, matches))
    outFile.close()
foo()  

   
