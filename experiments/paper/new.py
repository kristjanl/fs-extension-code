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
modelDir = os.path.join("..", "..", "models", "compositional", "paper", "stage2")
flowstar = os.path.join("..", "..", "src", "flowstar")



modelTypes = ["and_v3", "and_or_v2"]
#modelTypes = ["and_v3"]
#modelTypes = ["and_or_v2"]

comps = ["comp", "nocomp"]

groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
    
#discard gr2 and flow combo
#groups = [filter(lambda s: not "qr2" in s or not "flow" in s , g) for g in groups]    
    
    
models = [m for g in groups for m in g]


infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
    ("int time", "Int t"), 
    ("remap 1", "Remap1"), 
    ("evaluate t", "Eval@t"),
    ("precond time", "Precond t"), 
    ("remap 2", "Remap2"),
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_new")




















