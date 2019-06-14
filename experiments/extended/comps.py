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
modelTypes += ["lin_dep_20_id"]
modelTypes += ["and_or_v2_id"]
modelTypes += ["and_v3_id"]
modelTypes = []
modelTypes += ["pair_dep_10_id"]
modelTypes += ["jet_engine_10_pa"]
modelTypes += ["Brusselator_10_id"]
modelTypes += ["buckling_column_10_id"]
modelTypes += ["moore_rot_10_pa"]
modelTypes += ["cruise_control_10_id"]
modelTypes += ["diabetic_1_10_id"]
modelTypes += ["diabetic_2_10_id"]
modelTypes += ["filtered_oscillator_4_10_id"]
modelTypes += ["filtered_oscillator_8_10_id"]
modelTypes += ["filtered_oscillator_16_10_id"]
modelTypes += ["Lorentz_10_id"]
modelTypes += ["neuron_I_10_id"]
modelTypes += ["neuron_I_10_id"]
modelTypes += ["nonholonomic_10_id"]
modelTypes += ["rod_reactor_10_id"]
modelTypes += ["Roessler_10_id"]
modelTypes += ["switching_5_10_id"]
modelTypes += ["two_tanks_10_id"]
modelTypes += ["vehicle_platoon_3_10_id"]
modelTypes += ["vanderpol_10_nop"]

modelTypes += ["lin_dep_20_id"]
modelTypes += ["lin_10_id"]
modelTypes += ["bouncing_ball_10_id"]

#modelTypes = []
#modelTypes = ["lin_dep_20_id"]
#modelTypes += ["moore_rot_point_2_tqr"]

modelTypes = []
#modelTypes += ["jet_engine_10_pa"]
modelTypes += ["moore_rot_10_pa"]


#modelTypes = ["sq_deg_long_10_id"]
#modelTypes = ["Lotka_Volterra_10_pa"]
comps = ["fcomp", "lcomp", "nocomp"]
#comps = ["lcomp"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
#discard gr2 and fcomp combo
groups = [filter(lambda s: not "_pa_" in s or not "_fcomp" in s , g) for g in groups]
groups = [filter(lambda s: not "_nop_" in s or not "_lcomp" in s , g) for g in groups]
groups = [filter(lambda s: not "_tqr_" in s or not "_fcomp" in s , g) for g in groups]

groups[0].append("moore_rot_pa.model")
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
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps", vars=[0])
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps", vars=[0])




















