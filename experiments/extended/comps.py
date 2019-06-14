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

# directory of models
modelDir = os.path.join("..", "..", "models", "compositional", "extended", "comp")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = []
modelTypes += ["Lotka_Volterra_10_pa"]
modelTypes += ["sq_deg_long_10_id"]
modelTYpes += ["lin_dep_20_id"]
modelTypes += ["and_or_v2_id"]
modelTypes += ["and_v3_id"]
modelTypes += ["pair_dep_id"]
modelTypes += ["jet_engine_pa"]
modelTypes += ["Brusselator_id"]
modelTypes += ["buckling_column_id"]
modelTypes += ["moore_rot_pa"]
modelTypes += ["moore_rot_point_tqr"]
modelTypes += ["cruise_control_id"]
modelTypes += ["diabetic_1_id"]
modelTypes += ["diabetic_2_id"]
modelTypes += ["filtered_oscillator_4_id"]
modelTypes += ["filtered_oscillator_8_id"]
modelTypes += ["filtered_oscillator_16_id"]
modelTypes += ["Lorentz_id"]
modelTypes += ["neuron_I_id"]
modelTypes += ["neuron_I_id"]
modelTypes += ["nonholonomic_id"]
modelTypes += ["rod_reactor_id"]
modelTypes += ["Roessler_id"]
modelTypes += ["switching_5_id"]
modelTypes += ["two_tanks_id"]
modelTypes += ["vanderpol_nop"]
modelTypes += ["vehicle_platoon_3_id"]

modelTypes += ["lin_dep_id"]
modelTypes += ["lin_id"]
modelTypes += ["bouncing_ball_id"]



#modelTypes = ["sq_deg_long_10_id"]
#modelTypes = ["Lotka_Volterra_10_pa"]
comps = ["fcomp", "lcomp", "nocomp"]
#comps = ["fcomp"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
    
#discard gr2 and fcomp combo
groups = [filter(lambda s: not "_pa_" in s or not "_fcomp" in s , g) for g in groups]    

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
parser.add_argument('action', choices=['run', 'compare', 'both'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models, doLog=True)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps")
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps")




















