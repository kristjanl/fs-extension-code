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
modelDir = os.path.join("..", "..", "models", "compositional", "paper", "stage")
flowstar = os.path.join("..", "..", "src", "flowstar")


modelTypes = ["Lotka_Volterra_10_qr2"]
modelTypes += ["sq_deg_long_10_id"]
modelTypes += ["lin_dep_20_id"]
modelTypes += ["and_v3"]
modelTypes += ["and_or_v2"]

#modelTypes = ["sq_deg_long_10_id"]
#modelTypes += ["lin_dep_20_id"]
#modelTypes = ["Lotka_Volterra_10_qr2"]
#modelTypes = ["sq_deg_long_10_id"]


comps = ["fcomp", "comp", "nocomp"]
#comps = ["fcomp"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
    
#discard gr2 and fcomp combo
groups = [filter(lambda s: not "qr2" in s or not "fcomp" in s , g) for g in groups]    
    
    
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
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps")




















