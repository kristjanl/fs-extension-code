#!/usr/bin/env python 

import sympy

def statevars(vars, f):
  f.write("  state var ")
  for v in vars[0:-1]:
    f.write("%s,"%v)
  f.write("%s"%vars[-1])
  f.write('\n')
  
def no_comp(vars, f):
  f.write("    no decomposition\n")
def fullcomp(vars, f):
  f.write("    decomposition [")
  for v in vars[0:-1]:
    f.write("[%s],"%v)
  f.write("[%s]"%vars[-1])
  f.write("]\n")
  #f.write('\n')
def paircomp(vars, f):
  f.write("    decomposition [")
  pairs = [(vars[i],vars[i+1]) for i in range(0,len(vars),2)]
  for v1,v2 in pairs[0:-1]:
    f.write("[%s,%s],"%(v1,v2))
  f.write("[%s,%s]"%(pairs[-1][0],pairs[-1][1]))
  f.write("]\n")
  #f.write('\n')
  
def odes(vars, f):
  for v in vars:
    f.write("    %s' = -0.1*%s\n"%(v,v))
  #f.write('\n')
  
def lin(vars, f):
  for v in vars:
    f.write("    %s' = -0.1*%s\n"%(v,v))
    
def sqdeg(vars, f):
  for v in vars:
    f.write("    %s' = -0.1*%s*%s\n"%(v,v,v))
    
def lindep(vars, f):
  prev = None
  for v in vars:
    if prev == None:
      f.write("    %s' = -0.1*%s\n"%(v,v))
    else:
      f.write("    %s' = 0.1*%s - 0.1*%s\n"%(v,prev,v))
    prev = v
def pairdep(vars, f):
  pairs = [(vars[i],vars[i+1]) for i in range(0,len(vars),2)]
  for v1,v2 in pairs:
    f.write("    %s' = -0.1*%s*%s\n"%(v1,v1,v2))
    f.write("    %s' = -0.1*%s*%s\n"%(v2,v1,v2))

def initcond(vars, f):
  for v in vars:
    f.write("    %s in [0.5,1.0]\n"%v)

def name(vars, f, fileprefix, comp):
  #print ("    output %s_%s_%s\n"%(fileprefix, len(vars), comp))
  f.write("    output %s_%s_%s\n"%(fileprefix, len(vars), comp))

def foo(dim):
  vars = map(lambda i:"x%d"%(i+1), range(dim))
  
  
  fileprefix,odef,compf = "lin", lin, fullcomp
  #fileprefix,odef,compf = "lin_dep", lindep, fullcomp
  #fileprefix,odef,compf = "sq_deg", sqdeg, fullcomp
  #fileprefix,odef,compf = "sq_deg_long", sqdeg, fullcomp
  #fileprefix,odef,compf = "pair_dep", pairdep, paircomp
  
  print "%s, %s" %(dim, fileprefix);
  filename1 = 'models/%s_%d_nocomp.model'%(fileprefix, dim)
  f1 = open(filename1, 'w')
  lookup1 = {'$VARS$': statevars, '$COMP$':no_comp, '$ODES$': odef, '$INIT$':initcond, '$NAME$':(lambda v,f:name(v,f,fileprefix, "nocomp"))}
  
  filename2 = 'models/%s_%d_comp.model'%(fileprefix, dim)
  f2 = open(filename2, 'w')
  lookup2 = {'$VARS$': statevars, '$COMP$':compf, '$ODES$': odef, '$INIT$':initcond, '$NAME$':(lambda v,f:name(v,f,fileprefix, "comp"))}
  
  with open('template.model') as f3:
    for line in f3:
      if line.strip() in lookup1.keys():
        lookup1[line.strip()](vars, f1)
        lookup2[line.strip()](vars, f2)
      else:
        f1.write(line)
        f2.write(line)
  
  f1.close()
  f2.close()
  
  #print "============"
  #with open(filename1) as f:
  # for line in f:
  #     print line,
  
def bar():
  for i in range(1,51,1):
  #for i in range(2,51,2):
    foo(i)
  foo(100)