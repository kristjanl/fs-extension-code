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
modelDir = os.path.join("..", "..", "models", "compositional", "paper", "method_stage")
flowstar = os.path.join("..", "..", "src", "flowstar")


modelTypes = ["Brusselator", "Lotka_Volterra", "sq_deg_long", "lin_dep"]
modelTypes += ["and_or_v2"]
modelTypes += ["and_v3"]

algos = ["qrflow", "qrplain", "qr1", "qr2", "qr3", "id"]
#algos = ["qrplain"]
#algos = ["qr1", "qr2", "qr3"]
dims = ["nodim"]


groups = [["%s_nodim_%s.model"%(t,a) for a in algos] for t in modelTypes]
    
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
  comparer.compare5(modelDir, scriptsDir, groups, infoFields, suffix="_method")




















