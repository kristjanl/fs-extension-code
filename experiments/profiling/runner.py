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
modelDir = os.path.join("..", "..", "models", "profiling")
flowstar = os.path.join("..", "..", "src", "flowstar")

modelTypes = ["lower_path"]
modelTypes = ["full"]

param = ["nocomp", "comp"]

groups = [["%s_%s.model"%(m,p) for p in param] for m in modelTypes]
    
models = [m for g in groups for m in g]

infoFields = [
    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
    ("int_time", "Int t"), #("int time", "Int t"), 
    ("eval_t", "eval t"), 
    ("pic_poly", "p poly"),
    ("pic_decr", "p decr"),
    ("pic_ref", "p ref"),

    ("pic_ref_start", "refstart"), 
    ("pic_ref_first_picard", "reffirst"), 
    ("pic_ref_rem", "refrem"), 
    ("pic_ref_subset", "refsub"), 
]
infoFields = [
    ("pic_poly", "p poly"),
    ("pic_decr", "p decr"),
    ("pic_ref", "p ref"),

    ("pic_ref_start", "refstart"), 
    ("pic_ref_first", "reffirst"), 
    ("pic_ref_rem", "refrem"), 
    ("pic_ref_subset", "refsub"),
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare','both'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  print models
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)
elif args.action == 'both':
  print models
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)

























