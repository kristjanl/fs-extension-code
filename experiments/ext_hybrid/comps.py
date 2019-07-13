#!/usr/bin/env python 
import sys
import subprocess
import math
import re
import os
import argparse


# directory of scripts
scriptsDir = os.path.join("..", "scripts")
sys.path.append(scriptsDir)

import flowstar_runner
import comparer

import print_times

# directory of models
modelDir = os.path.join("..", "..", "models", "compositional", "extended", "comp2")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = []

modelTypes += ["bouncing_ball_id"]
modelTypes += ["cruise_control_id"]
modelTypes += ["diabetic_1_id"]
modelTypes += ["diabetic_2_id"]
modelTypes += ["filtered_oscillator_4_id"]
modelTypes += ["filtered_oscillator_8_id"]
modelTypes += ["filtered_oscillator_16_id"]
modelTypes += ["filtered_oscillator_32_id"]
modelTypes += ["nonholonomic_id"]
modelTypes += ["rod_reactor_id"]
modelTypes += ["two_tanks_id"]


#modelTypes = ["vanderpol_nop"]
#modelTypes += ["switching_5_id"]
#modelTypes += ["jet_engine_pa"]
#modelTypes += ["buckling_column_nop"]
#modelTypes = ["vehicle_platoon_3_nop"]
#modelTypes = ["filtered_oscillator_32_id"]

comps = ["fcomp", "lcomp", "nocomp"]
#comps = ["nocomp"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
#discard gr2 and fcomp combo
groups = [filter(lambda s: not "_pa_" in s or not "_fcomp" in s , g) for g in groups]
groups = [filter(lambda s: not "_nop_" in s or not "_lcomp" in s , g) for g in groups]
groups = [filter(lambda s: not "_tqr_" in s or not "_fcomp" in s , g) for g in groups]

models = [m for g in groups for m in g]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
    ("int_time", "Int t"),
    ("pic_poly", "p poly"),
    ("pic_decr", "p decr"),
    ("pic_ref", "p ref"),
    ("tr_precond", "Precond t"), 
    ("tr_remap1", "RM 1"), 
    ("tr_remap2", "RM 2"), 
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare', 'both', 'times'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models, doLog=True)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps", vars=[0])
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps", vars=[0])
elif args.action == 'times':
  print_times.times(modelDir, groups, infoFields)




















