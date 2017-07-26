#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os

# directory of scripts
scriptsDir = "scripts"
sys.path.append(scriptsDir)
import my_functions as fs
import table_functions as tf

# directory of models
modelDir = os.path.join("..", "models", "composition")
flowstar = os.path.join("..", "src", "flowstar")



models = []
models = models + ['lin_1_comp.model']
models = models + ['lin_1_nocomp.model']
models = models + ['lin_2_comp.model']
models = models + ['lin_2_nocomp.model']
models = models + ['lin_3_comp.model']
models = models + ['lin_3_nocomp.model']
models = models + ['lin_4_comp.model']
models = models + ['lin_4_nocomp.model']
models = models + ['lin_5_comp.model']
models = models + ['lin_5_nocomp.model']
models = models + ['lin_6_comp.model']
models = models + ['lin_6_nocomp.model']
models = models + ['lin_7_comp.model']
models = models + ['lin_7_nocomp.model']
models = models + ['lin_8_comp.model']
models = models + ['lin_8_nocomp.model']
models = models + ['lin_9_comp.model']
models = models + ['lin_9_nocomp.model']
models = models + ['lin_10_comp.model']
models = models + ['lin_10_nocomp.model']
models = models + ['lin_11_comp.model']
models = models + ['lin_11_nocomp.model']
models = models + ['lin_12_comp.model']
models = models + ['lin_12_nocomp.model']
models = models + ['lin_13_comp.model']
models = models + ['lin_13_nocomp.model']
models = models + ['lin_14_comp.model']
models = models + ['lin_14_nocomp.model']
models = models + ['lin_15_comp.model']
models = models + ['lin_15_nocomp.model']
models = models + ['lin_16_comp.model']
models = models + ['lin_16_nocomp.model']
models = models + ['lin_17_comp.model']
models = models + ['lin_17_nocomp.model']
models = models + ['lin_18_comp.model']
models = models + ['lin_18_nocomp.model']
models = models + ['lin_19_comp.model']
models = models + ['lin_19_nocomp.model']
models = models + ['lin_20_comp.model']
models = models + ['lin_20_nocomp.model']
models = models + ['lin_21_comp.model']
models = models + ['lin_21_nocomp.model']
models = models + ['lin_22_comp.model']
models = models + ['lin_22_nocomp.model']
models = models + ['lin_23_comp.model']
models = models + ['lin_23_nocomp.model']
models = models + ['lin_24_comp.model']
models = models + ['lin_24_nocomp.model']
models = models + ['lin_25_comp.model']
models = models + ['lin_25_nocomp.model']
models = models + ['lin_26_comp.model']
models = models + ['lin_26_nocomp.model']
models = models + ['lin_27_comp.model']
models = models + ['lin_27_nocomp.model']
models = models + ['lin_28_comp.model']
models = models + ['lin_28_nocomp.model']
models = models + ['lin_29_comp.model']
models = models + ['lin_29_nocomp.model']
models = models + ['lin_30_comp.model']
models = models + ['lin_30_nocomp.model']
models = models + ['lin_31_comp.model']
models = models + ['lin_31_nocomp.model']
models = models + ['lin_32_comp.model']
models = models + ['lin_32_nocomp.model']
models = models + ['lin_33_comp.model']
models = models + ['lin_33_nocomp.model']
models = models + ['lin_34_comp.model']
models = models + ['lin_34_nocomp.model']
models = models + ['lin_35_comp.model']
models = models + ['lin_35_nocomp.model']
models = models + ['lin_36_comp.model']
models = models + ['lin_36_nocomp.model']
models = models + ['lin_37_comp.model']
models = models + ['lin_37_nocomp.model']
models = models + ['lin_38_comp.model']
models = models + ['lin_38_nocomp.model']
models = models + ['lin_39_comp.model']
models = models + ['lin_39_nocomp.model']
models = models + ['lin_40_comp.model']
models = models + ['lin_40_nocomp.model']
models = models + ['lin_41_comp.model']
models = models + ['lin_41_nocomp.model']
models = models + ['lin_42_comp.model']
models = models + ['lin_42_nocomp.model']
models = models + ['lin_43_comp.model']
models = models + ['lin_43_nocomp.model']
models = models + ['lin_44_comp.model']
models = models + ['lin_44_nocomp.model']
models = models + ['lin_45_comp.model']
models = models + ['lin_45_nocomp.model']
models = models + ['lin_46_comp.model']
models = models + ['lin_46_nocomp.model']
models = models + ['lin_47_comp.model']
models = models + ['lin_47_nocomp.model']
models = models + ['lin_48_comp.model']
models = models + ['lin_48_nocomp.model']
models = models + ['lin_49_comp.model']
models = models + ['lin_49_nocomp.model']
models = models + ['lin_50_comp.model']
models = models + ['lin_50_nocomp.model']
models = models + ['lin_100_comp.model']
models = models + ['lin_100_nocomp.model']

models = models + ['lin_dep_1_comp.model']
models = models + ['lin_dep_1_nocomp.model']
models = models + ['lin_dep_2_comp.model']
models = models + ['lin_dep_2_nocomp.model']
models = models + ['lin_dep_3_comp.model']
models = models + ['lin_dep_3_nocomp.model']
models = models + ['lin_dep_4_comp.model']
models = models + ['lin_dep_4_nocomp.model']
models = models + ['lin_dep_5_comp.model']
models = models + ['lin_dep_5_nocomp.model']
models = models + ['lin_dep_6_comp.model']
models = models + ['lin_dep_6_nocomp.model']
models = models + ['lin_dep_7_comp.model']
models = models + ['lin_dep_7_nocomp.model']
models = models + ['lin_dep_8_comp.model']
models = models + ['lin_dep_8_nocomp.model']
models = models + ['lin_dep_9_comp.model']
models = models + ['lin_dep_9_nocomp.model']
models = models + ['lin_dep_10_comp.model']
models = models + ['lin_dep_10_nocomp.model']
models = models + ['lin_dep_11_comp.model']
models = models + ['lin_dep_11_nocomp.model']
models = models + ['lin_dep_12_comp.model']
models = models + ['lin_dep_12_nocomp.model']
models = models + ['lin_dep_13_comp.model']
models = models + ['lin_dep_13_nocomp.model']
models = models + ['lin_dep_14_comp.model']
models = models + ['lin_dep_14_nocomp.model']
models = models + ['lin_dep_15_comp.model']
models = models + ['lin_dep_15_nocomp.model']
models = models + ['lin_dep_16_comp.model']
models = models + ['lin_dep_16_nocomp.model']
models = models + ['lin_dep_17_comp.model']
models = models + ['lin_dep_17_nocomp.model']
models = models + ['lin_dep_18_comp.model']
models = models + ['lin_dep_18_nocomp.model']
models = models + ['lin_dep_19_comp.model']
models = models + ['lin_dep_19_nocomp.model']
models = models + ['lin_dep_20_comp.model']
models = models + ['lin_dep_20_nocomp.model']
models = models + ['lin_dep_21_comp.model']
models = models + ['lin_dep_21_nocomp.model']
models = models + ['lin_dep_22_comp.model']
models = models + ['lin_dep_22_nocomp.model']
models = models + ['lin_dep_23_comp.model']
models = models + ['lin_dep_23_nocomp.model']
models = models + ['lin_dep_24_comp.model']
models = models + ['lin_dep_24_nocomp.model']
models = models + ['lin_dep_25_comp.model']
models = models + ['lin_dep_25_nocomp.model']
models = models + ['lin_dep_26_comp.model']
models = models + ['lin_dep_26_nocomp.model']
models = models + ['lin_dep_27_comp.model']
models = models + ['lin_dep_27_nocomp.model']
models = models + ['lin_dep_28_comp.model']
models = models + ['lin_dep_28_nocomp.model']
models = models + ['lin_dep_29_comp.model']
models = models + ['lin_dep_29_nocomp.model']
models = models + ['lin_dep_30_comp.model']
models = models + ['lin_dep_30_nocomp.model']
models = models + ['lin_dep_31_comp.model']
models = models + ['lin_dep_31_nocomp.model']
models = models + ['lin_dep_32_comp.model']
models = models + ['lin_dep_32_nocomp.model']
models = models + ['lin_dep_33_comp.model']
models = models + ['lin_dep_33_nocomp.model']
models = models + ['lin_dep_34_comp.model']
models = models + ['lin_dep_34_nocomp.model']
models = models + ['lin_dep_35_comp.model']
models = models + ['lin_dep_35_nocomp.model']
models = models + ['lin_dep_36_comp.model']
models = models + ['lin_dep_36_nocomp.model']
models = models + ['lin_dep_37_comp.model']
models = models + ['lin_dep_37_nocomp.model']
models = models + ['lin_dep_38_comp.model']
models = models + ['lin_dep_38_nocomp.model']
models = models + ['lin_dep_39_comp.model']
models = models + ['lin_dep_39_nocomp.model']
models = models + ['lin_dep_40_comp.model']
models = models + ['lin_dep_40_nocomp.model']
models = models + ['lin_dep_41_comp.model']
models = models + ['lin_dep_41_nocomp.model']
models = models + ['lin_dep_42_comp.model']
models = models + ['lin_dep_42_nocomp.model']
models = models + ['lin_dep_43_comp.model']
models = models + ['lin_dep_43_nocomp.model']
models = models + ['lin_dep_44_comp.model']
models = models + ['lin_dep_44_nocomp.model']
models = models + ['lin_dep_45_comp.model']
models = models + ['lin_dep_45_nocomp.model']
models = models + ['lin_dep_46_comp.model']
models = models + ['lin_dep_46_nocomp.model']
models = models + ['lin_dep_47_comp.model']
models = models + ['lin_dep_47_nocomp.model']
models = models + ['lin_dep_48_comp.model']
models = models + ['lin_dep_48_nocomp.model']
models = models + ['lin_dep_49_comp.model']
models = models + ['lin_dep_49_nocomp.model']
models = models + ['lin_dep_50_comp.model']
models = models + ['lin_dep_50_nocomp.model']
models = models + ['lin_dep_100_comp.model']
models = models + ['lin_dep_100_nocomp.model']

models = models + ['sq_deg_1_comp.model']
models = models + ['sq_deg_1_nocomp.model']
models = models + ['sq_deg_2_comp.model']
models = models + ['sq_deg_2_nocomp.model']
models = models + ['sq_deg_3_comp.model']
models = models + ['sq_deg_3_nocomp.model']
models = models + ['sq_deg_4_comp.model']
models = models + ['sq_deg_4_nocomp.model']
models = models + ['sq_deg_5_comp.model']
models = models + ['sq_deg_5_nocomp.model']
models = models + ['sq_deg_6_comp.model']
models = models + ['sq_deg_6_nocomp.model']
models = models + ['sq_deg_7_comp.model']
models = models + ['sq_deg_7_nocomp.model']
models = models + ['sq_deg_8_comp.model']
models = models + ['sq_deg_8_nocomp.model']
models = models + ['sq_deg_9_comp.model']
models = models + ['sq_deg_9_nocomp.model']
models = models + ['sq_deg_10_comp.model']
models = models + ['sq_deg_10_nocomp.model']
models = models + ['sq_deg_11_comp.model']
models = models + ['sq_deg_11_nocomp.model']
models = models + ['sq_deg_12_comp.model']
models = models + ['sq_deg_12_nocomp.model']
models = models + ['sq_deg_13_comp.model']
models = models + ['sq_deg_13_nocomp.model']
models = models + ['sq_deg_14_comp.model']
models = models + ['sq_deg_14_nocomp.model']
models = models + ['sq_deg_15_comp.model']
models = models + ['sq_deg_15_nocomp.model']
models = models + ['sq_deg_16_comp.model']
models = models + ['sq_deg_16_nocomp.model']
models = models + ['sq_deg_17_comp.model']
models = models + ['sq_deg_17_nocomp.model']
models = models + ['sq_deg_18_comp.model']
models = models + ['sq_deg_18_nocomp.model']
models = models + ['sq_deg_19_comp.model']
models = models + ['sq_deg_19_nocomp.model']
models = models + ['sq_deg_20_comp.model']
models = models + ['sq_deg_20_nocomp.model']
models = models + ['sq_deg_21_comp.model']
models = models + ['sq_deg_21_nocomp.model']
models = models + ['sq_deg_22_comp.model']
models = models + ['sq_deg_22_nocomp.model']
models = models + ['sq_deg_23_comp.model']
models = models + ['sq_deg_23_nocomp.model']
models = models + ['sq_deg_24_comp.model']
models = models + ['sq_deg_24_nocomp.model']
models = models + ['sq_deg_25_comp.model']
models = models + ['sq_deg_25_nocomp.model']
models = models + ['sq_deg_26_comp.model']
models = models + ['sq_deg_26_nocomp.model']
models = models + ['sq_deg_27_comp.model']
models = models + ['sq_deg_27_nocomp.model']
models = models + ['sq_deg_28_comp.model']
models = models + ['sq_deg_28_nocomp.model']
models = models + ['sq_deg_29_comp.model']
models = models + ['sq_deg_29_nocomp.model']
models = models + ['sq_deg_30_comp.model']
models = models + ['sq_deg_30_nocomp.model']
models = models + ['sq_deg_31_comp.model']
models = models + ['sq_deg_31_nocomp.model']
models = models + ['sq_deg_32_comp.model']
models = models + ['sq_deg_32_nocomp.model']
models = models + ['sq_deg_33_comp.model']
models = models + ['sq_deg_33_nocomp.model']
models = models + ['sq_deg_34_comp.model']
models = models + ['sq_deg_34_nocomp.model']
models = models + ['sq_deg_35_comp.model']
models = models + ['sq_deg_35_nocomp.model']
models = models + ['sq_deg_36_comp.model']
models = models + ['sq_deg_36_nocomp.model']
models = models + ['sq_deg_37_comp.model']
models = models + ['sq_deg_37_nocomp.model']
models = models + ['sq_deg_38_comp.model']
models = models + ['sq_deg_38_nocomp.model']
models = models + ['sq_deg_39_comp.model']
models = models + ['sq_deg_39_nocomp.model']
models = models + ['sq_deg_40_comp.model']
models = models + ['sq_deg_40_nocomp.model']
models = models + ['sq_deg_41_comp.model']
models = models + ['sq_deg_41_nocomp.model']
models = models + ['sq_deg_42_comp.model']
models = models + ['sq_deg_42_nocomp.model']
models = models + ['sq_deg_43_comp.model']
models = models + ['sq_deg_43_nocomp.model']
models = models + ['sq_deg_44_comp.model']
models = models + ['sq_deg_44_nocomp.model']
models = models + ['sq_deg_45_comp.model']
models = models + ['sq_deg_45_nocomp.model']
models = models + ['sq_deg_46_comp.model']
models = models + ['sq_deg_46_nocomp.model']
models = models + ['sq_deg_47_comp.model']
models = models + ['sq_deg_47_nocomp.model']
models = models + ['sq_deg_48_comp.model']
models = models + ['sq_deg_48_nocomp.model']
models = models + ['sq_deg_49_comp.model']
models = models + ['sq_deg_49_nocomp.model']
models = models + ['sq_deg_50_comp.model']
models = models + ['sq_deg_50_nocomp.model']
models = models + ['sq_deg_100_comp.model']
models = models + ['sq_deg_100_nocomp.model']


models = models + ['sq_deg_long_1_comp.model']
models = models + ['sq_deg_long_1_nocomp.model']
models = models + ['sq_deg_long_2_comp.model']
models = models + ['sq_deg_long_2_nocomp.model']
models = models + ['sq_deg_long_3_comp.model']
models = models + ['sq_deg_long_3_nocomp.model']
models = models + ['sq_deg_long_4_comp.model']
models = models + ['sq_deg_long_4_nocomp.model']
models = models + ['sq_deg_long_5_comp.model']
models = models + ['sq_deg_long_5_nocomp.model']
models = models + ['sq_deg_long_6_comp.model']
models = models + ['sq_deg_long_6_nocomp.model']
models = models + ['sq_deg_long_7_comp.model']
models = models + ['sq_deg_long_7_nocomp.model']
models = models + ['sq_deg_long_8_comp.model']
models = models + ['sq_deg_long_8_nocomp.model']
models = models + ['sq_deg_long_9_comp.model']
models = models + ['sq_deg_long_9_nocomp.model']
models = models + ['sq_deg_long_10_comp.model']
models = models + ['sq_deg_long_10_nocomp.model']
models = models + ['sq_deg_long_11_comp.model']
models = models + ['sq_deg_long_11_nocomp.model']
models = models + ['sq_deg_long_12_comp.model']
models = models + ['sq_deg_long_12_nocomp.model']
models = models + ['sq_deg_long_13_comp.model']
models = models + ['sq_deg_long_13_nocomp.model']
models = models + ['sq_deg_long_14_comp.model']
models = models + ['sq_deg_long_14_nocomp.model']
models = models + ['sq_deg_long_15_comp.model']
models = models + ['sq_deg_long_15_nocomp.model']
models = models + ['sq_deg_long_16_comp.model']
models = models + ['sq_deg_long_16_nocomp.model']
models = models + ['sq_deg_long_17_comp.model']
models = models + ['sq_deg_long_17_nocomp.model']
models = models + ['sq_deg_long_18_comp.model']
models = models + ['sq_deg_long_18_nocomp.model']
models = models + ['sq_deg_long_19_comp.model']
models = models + ['sq_deg_long_19_nocomp.model']
models = models + ['sq_deg_long_20_comp.model']
models = models + ['sq_deg_long_20_nocomp.model']
models = models + ['sq_deg_long_21_comp.model']
models = models + ['sq_deg_long_21_nocomp.model']
models = models + ['sq_deg_long_22_comp.model']
models = models + ['sq_deg_long_22_nocomp.model']
models = models + ['sq_deg_long_23_comp.model']
models = models + ['sq_deg_long_23_nocomp.model']
models = models + ['sq_deg_long_24_comp.model']
models = models + ['sq_deg_long_24_nocomp.model']
models = models + ['sq_deg_long_25_comp.model']
models = models + ['sq_deg_long_25_nocomp.model']
models = models + ['sq_deg_long_26_comp.model']
models = models + ['sq_deg_long_26_nocomp.model']
models = models + ['sq_deg_long_27_comp.model']
models = models + ['sq_deg_long_27_nocomp.model']
models = models + ['sq_deg_long_28_comp.model']
models = models + ['sq_deg_long_28_nocomp.model']
models = models + ['sq_deg_long_29_comp.model']
models = models + ['sq_deg_long_29_nocomp.model']
models = models + ['sq_deg_long_30_comp.model']
models = models + ['sq_deg_long_30_nocomp.model']
models = models + ['sq_deg_long_31_comp.model']
models = models + ['sq_deg_long_31_nocomp.model']
models = models + ['sq_deg_long_32_comp.model']
models = models + ['sq_deg_long_32_nocomp.model']
models = models + ['sq_deg_long_33_comp.model']
models = models + ['sq_deg_long_33_nocomp.model']
models = models + ['sq_deg_long_34_comp.model']
models = models + ['sq_deg_long_34_nocomp.model']
models = models + ['sq_deg_long_35_comp.model']
models = models + ['sq_deg_long_35_nocomp.model']
models = models + ['sq_deg_long_36_comp.model']
models = models + ['sq_deg_long_36_nocomp.model']
models = models + ['sq_deg_long_37_comp.model']
models = models + ['sq_deg_long_37_nocomp.model']
models = models + ['sq_deg_long_38_comp.model']
models = models + ['sq_deg_long_38_nocomp.model']
models = models + ['sq_deg_long_39_comp.model']
models = models + ['sq_deg_long_39_nocomp.model']
models = models + ['sq_deg_long_40_comp.model']
models = models + ['sq_deg_long_40_nocomp.model']
models = models + ['sq_deg_long_41_comp.model']
models = models + ['sq_deg_long_41_nocomp.model']
models = models + ['sq_deg_long_42_comp.model']
models = models + ['sq_deg_long_42_nocomp.model']
models = models + ['sq_deg_long_43_comp.model']
models = models + ['sq_deg_long_43_nocomp.model']
models = models + ['sq_deg_long_44_comp.model']
models = models + ['sq_deg_long_44_nocomp.model']
models = models + ['sq_deg_long_45_comp.model']
models = models + ['sq_deg_long_45_nocomp.model']
models = models + ['sq_deg_long_46_comp.model']
models = models + ['sq_deg_long_46_nocomp.model']
models = models + ['sq_deg_long_47_comp.model']
models = models + ['sq_deg_long_47_nocomp.model']
models = models + ['sq_deg_long_48_comp.model']
models = models + ['sq_deg_long_48_nocomp.model']
models = models + ['sq_deg_long_49_comp.model']
models = models + ['sq_deg_long_49_nocomp.model']
models = models + ['sq_deg_long_50_comp.model']
models = models + ['sq_deg_long_50_nocomp.model']
models = models + ['sq_deg_long_100_comp.model']
models = models + ['sq_deg_long_100_nocomp.model']

models = models + ['pair_dep_2_comp.model']
models = models + ['pair_dep_2_nocomp.model']
models = models + ['pair_dep_4_comp.model']
models = models + ['pair_dep_4_nocomp.model']
models = models + ['pair_dep_6_comp.model']
models = models + ['pair_dep_6_nocomp.model']
models = models + ['pair_dep_8_comp.model']
models = models + ['pair_dep_8_nocomp.model']
models = models + ['pair_dep_10_comp.model']
models = models + ['pair_dep_10_nocomp.model']
models = models + ['pair_dep_12_comp.model']
models = models + ['pair_dep_12_nocomp.model']
models = models + ['pair_dep_14_comp.model']
models = models + ['pair_dep_14_nocomp.model']
models = models + ['pair_dep_16_comp.model']
models = models + ['pair_dep_16_nocomp.model']
models = models + ['pair_dep_18_comp.model']
models = models + ['pair_dep_18_nocomp.model']
models = models + ['pair_dep_20_comp.model']
models = models + ['pair_dep_20_nocomp.model']
models = models + ['pair_dep_22_comp.model']
models = models + ['pair_dep_22_nocomp.model']
models = models + ['pair_dep_24_comp.model']
models = models + ['pair_dep_24_nocomp.model']
models = models + ['pair_dep_26_comp.model']
models = models + ['pair_dep_26_nocomp.model']
models = models + ['pair_dep_28_comp.model']
models = models + ['pair_dep_28_nocomp.model']
models = models + ['pair_dep_30_comp.model']
models = models + ['pair_dep_30_nocomp.model']
models = models + ['pair_dep_32_comp.model']
models = models + ['pair_dep_32_nocomp.model']
models = models + ['pair_dep_34_comp.model']
models = models + ['pair_dep_34_nocomp.model']
models = models + ['pair_dep_36_comp.model']
models = models + ['pair_dep_36_nocomp.model']
models = models + ['pair_dep_38_comp.model']
models = models + ['pair_dep_38_nocomp.model']
models = models + ['pair_dep_40_comp.model']
models = models + ['pair_dep_40_nocomp.model']
models = models + ['pair_dep_42_comp.model']
models = models + ['pair_dep_42_nocomp.model']
models = models + ['pair_dep_44_comp.model']
models = models + ['pair_dep_44_nocomp.model']
models = models + ['pair_dep_46_comp.model']
models = models + ['pair_dep_46_nocomp.model']
models = models + ['pair_dep_48_comp.model']
models = models + ['pair_dep_48_nocomp.model']
models = models + ['pair_dep_50_comp.model']
models = models + ['pair_dep_50_nocomp.model']


regexp = re.compile(r'lin_[0-9]')
#regexp = re.compile(r'lin_dep_[0-9]')
#regexp = re.compile(r'sq_deg_[0-9]')
#regexp = re.compile(r'pair_dep_[0-9]')
comp = "_nocomp"


def getData(model):
  modelFile = os.path.join(modelDir, model)
  outputName = fs.getParam(modelFile, "output")
  infoFile = "infos/%s.txt" %outputName
  compTime = fs.getParam(infoFile, "computation time:")
  dim = tf.getDimension(modelFile)
  return compTime

def foo():
  prefix = "lin"
  #prefix = "lin_dep"
  #prefix = "sq_deg"
  #prefix = "sq_deg_long"
  #prefix = "pair_dep"
  fmaker = lambda prefix, dim, composition: "%s_%s_%s.model" %(prefix, dim, composition)
  print prefix
  for i in range(1, 51, 1):
    cfile = fmaker(prefix, i, "comp")
    nfile = fmaker(prefix, i, "nocomp")
    if not os.path.isfile(os.path.join(modelDir, cfile)):
      continue
    #print cfile
    #print nfile
    print "%s,%s,%s" %(i, getData(cfile), getData(nfile))
    
  
  cfile = fmaker(prefix, 100, "comp")
  nfile = fmaker(prefix, 100, "nocomp")
  print "%s,%s,%s" %(100, getData(cfile), getData(nfile))
  
foo()

"""
for model in models:
  if bool(regexp.search(model)) == False:
    continue
  if comp not in model:
    continue
  modelFile = os.path.join(modelDir, model)
  outputName = fs.getParam(modelFile, "output")
  infoFile = "infos/%s.txt" %outputName
  compTime = fs.getParam(infoFile, "computation time:")
  dim = tf.getDimension(modelFile)
  print "%s,%s" %(dim, compTime)
  #print model
"""
