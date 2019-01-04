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
#modelDir = os.path.join("..", "..", "models", "full_point_int")
modelDir = os.path.join("..", "..", "models", "temp5")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["temp", "and_v3", "and_or"]
modelTypes = ["and_or_full"]
#modelTypes = ["and_v3"]
#modelTypes = ["temp"]
modelTypes = ["sq_deg_long"]
modelTypes = ["mult_comp"]

param = ["pret", "pdis", "flow", "notfull"]
param = ["my", "flow"]
#param = ["my"]

groups = [["%s_%s.model"%(t,p) for p in param] for t in modelTypes]
groups = [["%s_%s_id.model"%(t,p) for p in param] for t in modelTypes]
    
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
parser.add_argument('action', choices=['run', 'compare', 'both'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)
elif args.action == 'both':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)
























