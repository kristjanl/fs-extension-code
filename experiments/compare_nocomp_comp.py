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

pairs = pairs + [('lin_3_comp.model', 'lin_3_nocomp.model')]
pairs = pairs + [('lin_20_comp.model', 'lin_20_nocomp.model')]


pairs = pairs + [('sq_deg_3_comp.model', 'sq_deg_3_nocomp.model')]
pairs = pairs + [('sq_deg_6_comp.model', 'sq_deg_6_nocomp.model')]
pairs = pairs + [('sq_deg_20_comp.model', 'sq_deg_20_nocomp.model')]

pairs = pairs + [('nl_4_comp.model', 'nl_4_nocomp.model')]
pairs = pairs + [('nl_20_comp.model', 'nl_20_nocomp.model')]

tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_nocomp_comp")

tf.write_table("experiments_nocomp_comp.html", modelDir, pairs, "_nocomp_comp")

  
  
  










