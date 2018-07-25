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
modelDir = os.path.join("..", "..", "models", "point_int")
flowstar = os.path.join("..", "..", "src", "flowstar2")

modelTypes = ["temp", "and_v3", "and_or"]
#modelTypes = ["and_or"]
#modelTypes = ["and_v3"]
modelTypes = ["temp"]

param = ["pret", "pdis", "flow"]

groups = [["%s_%s.model"%(t,p) for p in param] for t in modelTypes]
    
models = [m for g in groups for m in g]

infoFields = [

    ("int progress", "Int Progress"), 
    ("computation time", "Comp t"),  
    ("reason", "Stop Reason"), 
    ("int time", "Int t"), 
    ("eval t", "eval t"), 
    ("picard poly", "p poly"),
    ("picard decreasing", "p decr"),
    ("picard refining", "p ref"),
    ("precond time", "Precond t"), 

    ("refstart", "refstart"), 
    ("reffirst", "reffirst"), 
    ("refrem", "refrem"), 
    ("refsub", "refsub"), 

    ("pstart", "pstart"), 
    ("pm", "pm"), 
    ("pltr", "pltr"), 
    ("prrange", "prrange"), 
    ("pinsert", "pinsert"), 
    ("pscaling", "pscaling"), 
    ("plintra", "plintra"), 
    ("pleft", "pleft"), 
    ("pend", "pend"), 


    
]

parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['run', 'compare'])
parser.add_argument('rest', nargs=argparse.REMAINDER)
args = parser.parse_args()

if args.action == 'run':
  flowstar_runner.runFlowstar(modelDir, flowstar, models)
elif args.action == 'compare':
  comparer.compare5(modelDir, scriptsDir, groups, infoFields)
























