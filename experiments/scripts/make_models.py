#!/usr/bin/env python 
import sys
import subprocess
import math
import os

for m in sorted(os.listdir('../models/bench')):
  if not 'QR.' in m:
    continue
  print m
  outName = m.replace("QR.", "QR_fullflow.")
  outFile = open(os.path.join('../models/bench', outName), 'w')
  with open(os.path.join('../models/bench', m)) as f:
    for line in f:
      if 'alg_small_comp' in line:
        continue
      if 'no decomposition' in line:
        continue
      if 'out' in line:
        outFile.write(line.replace("QR", "QR_fullflow"))
        continue
      outFile.write(line)
    outFile.close()