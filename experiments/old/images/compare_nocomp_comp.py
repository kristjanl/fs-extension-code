#!/usr/bin/env python 
import sys
import os

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)

# directory of models
modelDir = os.path.join("..", "models", "composition")

import table_functions as tf

#for i in range(1, 51):
#    print "pairs = pairs + [('lin_%s_comp.model', 'lin_%s_nocomp.model')]"%(i,i)


pairs = []

pairs = pairs + [('lin_1_comp.model', 'lin_1_nocomp.model')]
pairs = pairs + [('lin_10_comp.model', 'lin_10_nocomp.model')]
pairs = pairs + [('lin_20_comp.model', 'lin_20_nocomp.model')]
pairs = pairs + [('lin_30_comp.model', 'lin_30_nocomp.model')]
pairs = pairs + [('lin_40_comp.model', 'lin_40_nocomp.model')]
pairs = pairs + [('lin_50_comp.model', 'lin_50_nocomp.model')]
pairs = pairs + [('lin_100_comp.model', 'lin_100_nocomp.model')]


pairs = pairs + [('lin_dep_1_comp.model', 'lin_dep_1_nocomp.model')]
pairs = pairs + [('lin_dep_10_comp.model', 'lin_dep_10_nocomp.model')]
pairs = pairs + [('lin_dep_20_comp.model', 'lin_dep_20_nocomp.model')]
pairs = pairs + [('lin_dep_30_comp.model', 'lin_dep_30_nocomp.model')]
pairs = pairs + [('lin_dep_40_comp.model', 'lin_dep_40_nocomp.model')]
pairs = pairs + [('lin_dep_50_comp.model', 'lin_dep_50_nocomp.model')]
pairs = pairs + [('lin_dep_100_comp.model', 'lin_dep_100_nocomp.model')]


pairs = pairs + [('sq_deg_1_comp.model', 'sq_deg_1_nocomp.model')]
pairs = pairs + [('sq_deg_10_comp.model', 'sq_deg_10_nocomp.model')]
pairs = pairs + [('sq_deg_20_comp.model', 'sq_deg_20_nocomp.model')]
pairs = pairs + [('sq_deg_30_comp.model', 'sq_deg_30_nocomp.model')]
pairs = pairs + [('sq_deg_40_comp.model', 'sq_deg_40_nocomp.model')]
pairs = pairs + [('sq_deg_50_comp.model', 'sq_deg_50_nocomp.model')]
pairs = pairs + [('sq_deg_100_comp.model', 'sq_deg_100_nocomp.model')]

pairs = pairs + [('sq_deg_long_1_comp.model', 'sq_deg_long_1_nocomp.model')]
pairs = pairs + [('sq_deg_long_10_comp.model', 'sq_deg_long_10_nocomp.model')]
pairs = pairs + [('sq_deg_long_20_comp.model', 'sq_deg_long_20_nocomp.model')]
pairs = pairs + [('sq_deg_long_30_comp.model', 'sq_deg_long_30_nocomp.model')]
pairs = pairs + [('sq_deg_long_40_comp.model', 'sq_deg_long_40_nocomp.model')]
pairs = pairs + [('sq_deg_long_50_comp.model', 'sq_deg_long_50_nocomp.model')]
pairs = pairs + [('sq_deg_long_100_comp.model', 'sq_deg_long_100_nocomp.model')]

pairs = pairs + [('pair_dep_2_comp.model', 'pair_dep_2_nocomp.model')]
pairs = pairs + [('pair_dep_10_comp.model', 'pair_dep_10_nocomp.model')]
pairs = pairs + [('pair_dep_20_comp.model', 'pair_dep_20_nocomp.model')]
pairs = pairs + [('pair_dep_30_comp.model', 'pair_dep_30_nocomp.model')]
pairs = pairs + [('pair_dep_40_comp.model', 'pair_dep_40_nocomp.model')]
pairs = pairs + [('pair_dep_50_comp.model', 'pair_dep_50_nocomp.model')]
pairs = pairs + [('pair_dep_100_comp.model', 'pair_dep_100_nocomp.model')]

tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_nocomp_comp")

tf.write_table("experiments_nocomp_comp.html", modelDir, pairs, "_nocomp_comp")

  
  
  










