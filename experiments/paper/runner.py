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


modelTypes = ["Lotka_Volterra_1_qr2", "Lotka_Volterra_5_qr2", "Lotka_Volterra_10_qr2"]
#modelTypes = ["Lotka_Volterra_1_qr2", "Lotka_Volterra_10_qr2"]

comps = ["comp", "nocomp"]

models = ["%s_%s.model"%(modelType, comp) \
    for modelType in modelTypes \
    for comp in comps]
    
    
groups = [["%s_%s.model"%(modelType, comp) for comp in comps] \
    for modelType in modelTypes]
    
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
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)




















