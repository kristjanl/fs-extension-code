#!/usr/bin/env python 
import sys
import os

import table_functions as tf

def compare(modelDir, scriptsDir, modelTypes, dims, type1, type2, infoFields):
  compareComb = "_%s_%s" %(type1,type2)
  nameF = lambda m, n, t: "%s_%s_id_%s.model"%(m, n, t)

  pairs = []
  for m in modelTypes:
    for i in dims:
      pairs.append( (nameF(m, i, type1), nameF(m, i, type2)) )
  
  tf.generateHtml(scriptsDir, modelDir, pairs, compareComb, infoFields)

  
  
  










