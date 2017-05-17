#!/usr/bin/env python 
import sys
import subprocess
import math
import os

for m in sorted(os.listdir('../models/bench')):
  if not 'QR.' in m:
    continue
  print m
  outName = m.replace("QR.", "QR_flow.")
  outFile = open(os.path.join('../models/bench', outName), 'w')
  with open(os.path.join('../models/bench', m)) as f:
    for line in f:
      if 'alg_small_comp' in line:
        print line
        outFile.write(line.replace("alg_small_comp", "alg_small_comp flow impl"))
        continue
      if 'out' in line:
        outFile.write(line.replace("QR", "QR_flow"))
        continue
      outFile.write(line)
    outFile.close()
