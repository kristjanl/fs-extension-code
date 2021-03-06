import sys
import subprocess
import math
import os
import re
  
def find_file_max_time(fileName):
  max = -1
  if not os.path.isfile(fileName):
    return None
    
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
        
def get_range_bounds(time, csvs, modelFiles):
  all = map(lambda s: get_single_range_bounds(time, s), csvs)
  
  #filter out empty data
  all = filter(lambda x:x!=[], all)
  
  #for d in all:
  #  print d[1][:9]
  #  print d[1]
  
  if len(all) == 0:
    raise Exception("empty ranges for models %s"%modelFiles)
  minBounds = all[0][0]
  maxBounds = all[0][1]
  
  for (singleMin, singleMax) in all[1:]:
    #print len(singleMin)
    #print len(minBounds)
    for (i, d) in enumerate(singleMin):
      if minBounds[i] > d:
        minBounds[i] = d
    for (i, d) in enumerate(singleMax):
      if maxBounds[i] < d:
        maxBounds[i] = d
  for (i, _) in enumerate(minBounds):
    if(i == 0):
      continue
    size = maxBounds[i] - minBounds[i]
    minBounds[i] = minBounds[i] - size/10
    maxBounds[i] = maxBounds[i] + size/10
  
  #print maxBounds[:6]
  return [minBounds,maxBounds]
  
def get_single_range_bounds(time, csv):
  
  if not os.path.isfile(csv):
    return []
  minBounds = []
  maxBounds = []

  lastMax = 1
  with open(csv) as f:
    #print csv
    first = None
    for line in f:
      #-1 since last one is newline
      data = line.split(',')[:-1]
      dData = map(lambda x: float(x), data)
      if minBounds == []:
        minBounds = list(dData)
        maxBounds = list(dData)
        continue
      
      #print dData
      toInspect =  (len(dData) - 2)/2
      maxWidth = max(dData[-toInspect:])
      if first == None and maxWidth > 0.1:
        first = maxWidth
      #print maxWidth / lastMax
      #print maxWidth / first
      if maxWidth / lastMax > 1000 or (first != None and maxWidth / first > 100):
        #print "breaking"
        break
        pass
      lastMax = maxWidth
      for (i, d) in enumerate(dData):
        #print (i, d)
        #print "(%s, %s, %s, %s)" %(i, d, minBounds[i], minBounds[i] > d)
        if minBounds[i] > d:
          minBounds[i] = d
        if maxBounds[i] < d:
          maxBounds[i] = d
      
      if float(data[0]) >= float(time):
        break
  return [minBounds, maxBounds]

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
    if exp == None:
      continue
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
        return m.group(1).strip()
  
  
  
  
