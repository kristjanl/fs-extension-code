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
pairs = pairs + [('lin_2_comp.model', 'lin_2_nocomp.model')]
pairs = pairs + [('lin_3_comp.model', 'lin_3_nocomp.model')]
pairs = pairs + [('lin_4_comp.model', 'lin_4_nocomp.model')]
pairs = pairs + [('lin_5_comp.model', 'lin_5_nocomp.model')]
pairs = pairs + [('lin_6_comp.model', 'lin_6_nocomp.model')]
pairs = pairs + [('lin_7_comp.model', 'lin_7_nocomp.model')]
pairs = pairs + [('lin_8_comp.model', 'lin_8_nocomp.model')]
pairs = pairs + [('lin_9_comp.model', 'lin_9_nocomp.model')]
pairs = pairs + [('lin_10_comp.model', 'lin_10_nocomp.model')]
pairs = pairs + [('lin_11_comp.model', 'lin_11_nocomp.model')]
pairs = pairs + [('lin_12_comp.model', 'lin_12_nocomp.model')]
pairs = pairs + [('lin_13_comp.model', 'lin_13_nocomp.model')]
pairs = pairs + [('lin_14_comp.model', 'lin_14_nocomp.model')]
pairs = pairs + [('lin_15_comp.model', 'lin_15_nocomp.model')]
pairs = pairs + [('lin_16_comp.model', 'lin_16_nocomp.model')]
pairs = pairs + [('lin_17_comp.model', 'lin_17_nocomp.model')]
pairs = pairs + [('lin_18_comp.model', 'lin_18_nocomp.model')]
pairs = pairs + [('lin_19_comp.model', 'lin_19_nocomp.model')]
pairs = pairs + [('lin_20_comp.model', 'lin_20_nocomp.model')]
pairs = pairs + [('lin_21_comp.model', 'lin_21_nocomp.model')]
pairs = pairs + [('lin_22_comp.model', 'lin_22_nocomp.model')]
pairs = pairs + [('lin_23_comp.model', 'lin_23_nocomp.model')]
pairs = pairs + [('lin_24_comp.model', 'lin_24_nocomp.model')]
pairs = pairs + [('lin_25_comp.model', 'lin_25_nocomp.model')]
pairs = pairs + [('lin_26_comp.model', 'lin_26_nocomp.model')]
pairs = pairs + [('lin_27_comp.model', 'lin_27_nocomp.model')]
pairs = pairs + [('lin_28_comp.model', 'lin_28_nocomp.model')]
pairs = pairs + [('lin_29_comp.model', 'lin_29_nocomp.model')]
pairs = pairs + [('lin_30_comp.model', 'lin_30_nocomp.model')]
pairs = pairs + [('lin_31_comp.model', 'lin_31_nocomp.model')]
pairs = pairs + [('lin_32_comp.model', 'lin_32_nocomp.model')]
pairs = pairs + [('lin_33_comp.model', 'lin_33_nocomp.model')]
pairs = pairs + [('lin_34_comp.model', 'lin_34_nocomp.model')]
pairs = pairs + [('lin_35_comp.model', 'lin_35_nocomp.model')]
pairs = pairs + [('lin_36_comp.model', 'lin_36_nocomp.model')]
pairs = pairs + [('lin_37_comp.model', 'lin_37_nocomp.model')]
pairs = pairs + [('lin_38_comp.model', 'lin_38_nocomp.model')]
pairs = pairs + [('lin_39_comp.model', 'lin_39_nocomp.model')]
pairs = pairs + [('lin_40_comp.model', 'lin_40_nocomp.model')]
pairs = pairs + [('lin_41_comp.model', 'lin_41_nocomp.model')]
pairs = pairs + [('lin_42_comp.model', 'lin_42_nocomp.model')]
pairs = pairs + [('lin_43_comp.model', 'lin_43_nocomp.model')]
pairs = pairs + [('lin_44_comp.model', 'lin_44_nocomp.model')]
pairs = pairs + [('lin_45_comp.model', 'lin_45_nocomp.model')]
pairs = pairs + [('lin_46_comp.model', 'lin_46_nocomp.model')]
pairs = pairs + [('lin_47_comp.model', 'lin_47_nocomp.model')]
pairs = pairs + [('lin_48_comp.model', 'lin_48_nocomp.model')]
pairs = pairs + [('lin_49_comp.model', 'lin_49_nocomp.model')]
pairs = pairs + [('lin_50_comp.model', 'lin_50_nocomp.model')]


"""
pairs = pairs + [('sq_deg_3_comp.model', 'sq_deg_3_nocomp.model')]
pairs = pairs + [('sq_deg_6_comp.model', 'sq_deg_6_nocomp.model')]
pairs = pairs + [('sq_deg_20_comp.model', 'sq_deg_20_nocomp.model')]

pairs = pairs + [('nl_4_comp.model', 'nl_4_nocomp.model')]
pairs = pairs + [('nl_20_comp.model', 'nl_20_nocomp.model')]
"""
tf.generate_comparision_plots(scriptsDir, modelDir, pairs, "_nocomp_comp")

tf.write_table("experiments_nocomp_comp.html", modelDir, pairs, "_nocomp_comp")

  
  
  










