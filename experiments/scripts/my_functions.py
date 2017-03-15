import sys
import subprocess
import math
import os
  
def find_file_max_time(fileName):
  max = -1
  if not os.path.isfile(fileName):
    return -1
    
  with open(fileName) as inFile:
    for line in inFile:
      data = line.split(",") 
      time = float(data[0])
      max = time
  return max

def find_max_time(f1, f2):
  time1 = find_file_max_time(f1)
  time2 = find_file_max_time(f2)
  max_common = "%s" %min(time1, time2)
  return max_common
  
