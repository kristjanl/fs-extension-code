#!/usr/bin/env python 
import sys
import os

# directory of scripts
scriptsDir = os.path.join("..", "scripts")
sys.path.append(scriptsDir)

# directory of models
modelDir = os.path.join("..", ".." , "models", "compositional", "id")

import table_functions as tf

pairs = []


type1 = "flow"
type2 = "comp"

compareComb = "_%s_%s" %(type1,type2)

nameF = lambda m, n, t: "%s_%s_id_%s.model"%(m, n, t)

modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg", "sq_deg_long"]
dims = [2]# + range(10, 51, 10)
dims = [2, 20]

for m in modelTypes:
  for i in dims:
    pairs = pairs + [(nameF(m, i, type1), nameF(m, i, type2))]


tf.generateHtml(scriptsDir, modelDir, pairs, compareComb)

  
  
  










