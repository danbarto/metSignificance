import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
from helpers import *
from math import *

import json

outfile = 'MC_tune.txt'

# Define working points etc
presel = 'Sum$(jet_pt>0)>1'
sigCut = 9.

#samples = [WW,WZ,ZZ]#,ST_top,ST_antitop]
samples = allMCSamples
#samples = [data]
isData = False

# load chain to list
el = eventlist( samples, presel )

el.doJetSmearing(False)
if isData: #protection for data
  el.doJetSmearing(False)

# Do the minimization
print 'Total events:',len(el.evlist)

gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
gmin.SetTolerance(10.0)
gmin.SetStrategy(0)
gmin.SetPrintLevel(3)

LL = minLL(el)

gmin.SetFunction(LL)

variable  = ['a1','a2','a3','a4','a5','N1','S1']
step      = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
start     = [1.0,  1.0,  1.0,  1.0, 1.0, 0., 0.5]

print 'Minimizing parameters',variable
print 'With stepsize of', step
print 'And starting values', start


for i in range(0,7):
  gmin.SetVariable(i,variable[i],start[i], step[i])

gmin.Minimize()

#save full event list
full_el = el

#filter events with high significance
el.evlist = [x  for x in el.evlist if x.sig < sigCut]

print
print 'Now fitting after applying significance cut'
print 'Total events:',len(el.evlist)

gmin.SetStrategy(1)
gmin.Minimize()
gmin.Hesse()

pars = [gmin.X()[i] for i in range(0,7)]

with open(outfile, 'w') as of:
  json.dump(pars, of)
