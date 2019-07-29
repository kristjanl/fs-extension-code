import sys
import math
import os
import re
import subprocess

import my_functions as fs
import table_functions as tf

infos =  ["computation time", "int_time", "pic_poly", "pic_decr", "pic_ref", "tr_remap1", "tr_precond", "tr_remap2"]

def get_processing_method(modelFile):
  if "_id_" in modelFile:
    return "ID"
  if "_pa_" in modelFile:
    return "PA"
  if "_nop_" in modelFile:
    return "NO"
  raise Exception("Couldn't determine processing method (%s)"%modelFile)

def get_comp(s):
  compLookup = {"fcomp": "FC", "lcomp":"LC", "nocomp": "NC"}
  if not s in compLookup:
    raise Exception("Couldn't determine composition type (%s)"%s)
  return compLookup[s]

def get_model_alias(s):
  if not s in tf.nameLookup:
    raise Exception("Couldn't determine alias (%s)"%s)
  return tf.nameLookup[s]

data = {}

def times(modelDir, groups, infoFields):
  for models in groups:
    for model in models:
      modelFile = os.path.join(modelDir, model)
      #print model

      m = re.search('(.*?)_(10_|20_){0,1}(nop|id|pa)_(.*)\.', model)
      modelPrefix = m.group(1)
      alias = get_model_alias(modelPrefix)
      compType = get_comp(m.group(4))

      
      outputName = fs.getParam(modelFile, "output")
      infoFile = "infos/%s.txt" %outputName
      #print outputName
      #print infoFile
      reason = fs.getParam(infoFile, "reason:")
      #print reason
      timeLookup = {}
      dim = tf.getDimension(modelFile)
      processing = get_processing_method(modelFile)
      #print processing
      for info in infos:
        value = fs.getParam(infoFile, "%s:"%info)
        if value == None:
          timeLookup[info] = None
          continue
        timeLookup[info] = float(value)
        #timeLookup[info] = value
      #print timeLookup
      #print dim

      if modelPrefix not in data:
        data[modelPrefix] = {}
      data[modelPrefix][compType] = [modelPrefix, alias, processing, compType, timeLookup]
  
  for l in [tf.cont, tf.arti, tf.hybr]:
    #print data
    compTypes = ['FC', 'LC', 'NC']
    for key in l:
      #print key 
      if key not in data:
        #print "--------------missing key--------------" + key
        continue
      for c in compTypes:
        if c not in data[key]:
          continue
        #printLine(data[key], c)
        printPercentage(data[key], c)
      print "\\hline"
      #sys.exit(0)

def formatTime(time):
  if time == None:
    return "--"
  return "{0:.3f}".format(float(time))
  #return "formatted"
  return time

def divFormatWrapper(t1, t2):
  if t1 == None:
    return "--"
  """
  if t1 == None and t2 == None:
    return "--(--)"
  elif t1 == None:
    return "--(%ss)"%(formatTime(t2))
  """
  return formatTime(t1 / t2) + "\\%"

def printLine(modelData, c):
  (modelPrefix, alias, processing, compType, timeLookup) = modelData[c]
  nameCol = ""
  processingCol = ""
  if (c == "FC" and "FC" in modelData) or (c == "LC" and "FC" not in modelData):
    nameCol = "\multirow{%s}{*}{%s}"%(len(modelData), alias)
    processingCol = "\multirow{%s}{*}{%s}"%(len(modelData), processing)
  line = "%s & %s & %s"%(nameCol, processingCol, compType)
  for t in infos:
    line += " & %s"%formatTime(timeLookup[t])
  line += " \\\\"
  print line


def printPercentage(modelData, c):
  (modelPrefix, alias, processing, compType, timeLookup) = modelData[c]
  nameCol = ""
  processingCol = ""
  if c == "NC":
    return
  if (c == "FC" and "FC" in modelData) or (c == "LC" and "FC" not in modelData):
    nameCol = "\multirow{%s}{*}{%s}"%(len(modelData) - 1, alias)
    processingCol = "\multirow{%s}{*}{%s}"%(len(modelData) - 1, processing)
  line = "%s & %s & %s"%(nameCol, processingCol, compType)
  for t in infos:
    #line += " & %s"%formatTime(timeLookup[t])
    line += " & %s"%(divFormatWrapper(timeLookup[t], modelData["NC"][4][t]))
  line += " \\\\"
  print line