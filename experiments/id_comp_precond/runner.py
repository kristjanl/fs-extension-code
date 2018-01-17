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
modelDir = os.path.join("..", "..", "models", "compositional", "id_comp_precond")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg", "sq_deg_long"]
modelTypes = ["lin", "lin_dep", "pair_dep", "sq_deg_long"]
#modelTypes = ["lin"]
#modelTypes = ["pair_dep"]
#modelTypes = ["lin", "pair_dep"]
#modelTypes = ["sq_deg_long"]
#modelTypes = ["lin_dep"]

dims = [2] + range(10, 51, 10)
#dims = [40]
algos = ["id_flow", "id_comp", "cid_comp"]
algos = ["id_comp", "cid_comp"]

algos = ["id_flow"]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"), 
    ("int time", "Int t"), 
    ("remap 1", "Remap1"), 
    ("evaluate t", "Eval@t"),
    ("precond time", "Precond t"), 
    ("remap 2", "Remap2"),
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
  if len(set(["id_flow", "id_comp","cid_comp"]).intersection(set(args.rest))) != 2\
      or len(args.rest) != 2:
    raise ValueError("compare needs flow, comp or nocomp as arguments")
  comparer.compare2(modelDir, scriptsDir, modelTypes, dims, args.rest[0], \
      args.rest[1], infoFields)




















