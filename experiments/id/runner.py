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
modelDir = os.path.join("..", "..", "models", "compositional", "id")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg", "sq_deg_long"]
modelTypes = ["pair_dep"]
dims = [2] + range(10, 51, 10)
dims = [2]
algos = ["flow", "comp", "nocomp"]
algos = ["flow", "comp"]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Computation time"), 
    ("int time", "Integration time"), 
    ("remap 1", "Remap 1"), 
    ("evaluate t", "Evaluate at t"),
    ("precond time", "Precondition time"), 
    ("remap 2", "Remap 2"),
]

models = ["%s_%s_id_%s.model"%(modelType, dim, algo) \
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
  if len(set(["flow", "comp","nocomp"]).intersection(set(args.rest))) != 2 or\
      len(args.rest) != 2:
    raise ValueError("compare needs flow, comp or nocomp as arguments")
  comparer.compare(modelDir, scriptsDir, modelTypes, dims, args.rest[0], \
      args.rest[1], infoFields)




















