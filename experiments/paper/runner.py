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

modelTypes = []
#modelTypes = ["and_fast_out_high_10_id", "Lotka_Volterra_10_qr2"]
#modelTypes = ["Lotka_Volterra_10_qr2"]
#modelTypes = ["and_fast_out_high_10_id"]
modelTypes += ["lin_10_id", "pair_dep_10_id", "sq_deg_long_10_id"]
modelTypes += ["lin_10_qr2", "pair_dep_10_qr2", "sq_deg_long_10_qr2"]
modelTypes += ["lin_dep_20_id"]

#modelTypes = ["jet_engine_10_qr2"]
#modelTypes = ["and_fast_out_high_10_id"]

comps = ["comp", "nocomp", "flow"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
    
#discard gr2 and flow combo
groups = [filter(lambda s: not "qr2" in s or not "flow" in s , g) for g in groups]    
    
    
models = [m for g in groups for m in g]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
    ("int time", "Int t"), 
    ("picard poly", "p poly"),
    ("picard decreasing", "p decr"),
    ("picard refining", "p ref"),
    ("precond time", "Precond t"), 
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_comps")




















