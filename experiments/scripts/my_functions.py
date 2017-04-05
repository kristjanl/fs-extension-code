import sys
import subprocess
import math
import os
import re
  
def find_file_max_time(fileName):
  max = -1
  if not os.path.isfile(fileName):
    return 100000000000000
    
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


def get_range_at_time(time, csv):
  if not os.path.isfile(csv):
    return []
  with open(csv) as f:
    for line in f:
      #-1 since last one is newline
      data = line.split(',')[:-1]
      #line with needed time
      if float(data[0]) == float(time):
        #get variable ranges
        return map(lambda x: float(x), data)
        


def get_range_up_to(time, csv):
  if not os.path.isfile(csv):
    return []
  ret = []
  with open(csv) as f:
    for line in f:
      #-1 since last one is newline
      data = line.split(',')[:-1]
      #line with needed time
      ret.append(map(lambda x: float(x), data))
      if float(data[0]) == float(time):
        return ret
        
def get_bounds(values):
  #filter out when no data
  minVar = maxVar = None
  
  minBounds = []
  maxBounds = []
  for exp in values:
    for data in exp:
      if len(minBounds) == 0:
        minBounds = list(data)
        maxBounds = data
        continue
      
      for (i,d) in enumerate(data):
        if minBounds[i] > d:
          minBounds[i] = d
        if maxBounds[i] < d:
          maxBounds[i] = d
  return [minBounds, maxBounds]
  

def getParam(filename, param):
  with open(filename) as f:
    for line in f:
      #print line,
      m = re.search('%s (.*)' %param, line)
      if m != None:
        return m.group(1)
  
  
  
  
