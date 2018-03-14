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
modelDir = os.path.join("..", "..", "models", "compositional", "dqr", "stage")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["Brusselator"]


algos = ["qrflow"]#, "qr1", "qr2", "qr3"]
#algos = ["qr1", "qr2", "qr3"]
dims = ["nodim"]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
#    ("int time", "Int t"), 
#    ("remap 1", "Remap1"), 
#    ("evaluate t", "Eval@t"),
#    ("precond time", "Precond t"), 
#    ("remap 2", "Remap2"),
]

models = ["%s_%s_%s.model"%(modelType, dim, algo) \
    for modelType in modelTypes \
    for dim in dims \
    for algo in algos]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare2(modelDir, scriptsDir, modelTypes, dims, "qrflow", \
      "qrflow", infoFields)




















