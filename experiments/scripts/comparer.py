#!/usr/bin/env python 
import sys
import os

import table_functions as tf


#models constructed will have _id in the name
def compare(modelDir, scriptsDir, modelTypes, dims, type1, type2, infoFields):
  compareComb = "_%s_%s" %(type1,type2)
  nameF = lambda m, n, t: "%s_%s_id_%s.model"%(m, n, t)

  pairs = []
  for m in modelTypes:
    for i in dims:
      pairs.append( (nameF(m, i, type1), nameF(m, i, type2)) )
  
  tf.generateHtml(scriptsDir, modelDir, pairs, compareComb, infoFields)
#models constructed don't assume _id in the name
def compare2(modelDir, scriptsDir, modelTypes, dims, type1, type2, infoFields):
  compareComb = "_%s_%s" %(type1,type2)
  nameF = lambda m, n, t: "%s_%s_%s.model"%(m, n, t)

  pairs = []
  for m in modelTypes:
    for i in dims:
      pairs.append( (nameF(m, i, type1), nameF(m, i, type2)) )
  
  tf.generateHtml(scriptsDir, modelDir, pairs, compareComb, infoFields)
  
def compare3(modelDir, scriptsDir, modelTypes, dims, algos, infoFields):
  compareComb = reduce(lambda p,n: p + "_" + n, algos, "")
  nameF = lambda m, n, t: "%s_%s_%s.model"%(m, n, t)

  groups = []
  for m in modelTypes:
    for i in dims:
      groups.append(map(lambda s: nameF(m, i, s), algos))
  
  tf.generateGroupHtml(scriptsDir, modelDir, groups, compareComb, infoFields)
  

#no names in html file
def compare4(modelDir, scriptsDir, modelTypes, dims, algos, infoFields, suffix=""):
  compareComb = reduce(lambda p,n: p + "_" + n, algos, "")
  compareComb = ""
  nameF = lambda m, n, t: "%s_%s_%s.model"%(m, n, t)

  groups = []
  for m in modelTypes:
    if dims != None:
      for i in dims:
        groups.append(map(lambda s: nameF(m, i, s), algos))
    else:
      groups.append( map(lambda s: "%s_%s.model"%(m, s), algos) )
  
  tf.generateGroupHtml(scriptsDir, modelDir, groups, suffix, infoFields)
    
  










