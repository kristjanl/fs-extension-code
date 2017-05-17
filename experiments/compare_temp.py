#!/usr/bin/env python 
import sys
import os

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)

# directory of models
modelDir = os.path.join("..", "models", "bench")

import table_functions as tf

pairs = []

pairs = pairs + [('moore_rot_QR.model', 'moore_rot_QR2.model')]

tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_temp")

tf.write_table("temp.html", modelDir, pairs, "_temp")

  
  
  










