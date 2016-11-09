import ROOT

from loadAllSamples import *
from eventlist import *
from likelihood import *
from helpers import *
from math import *

import json

def main():

  outfile = 'data/data_tune_old_2jet15.txt'
  
  # Define working points etc
  presel = 'Sum$(jet_pt>15&&abs(jet_eta)<2.5&&jet_passid)>=2'
  sigCut = 9.
  
  #samples = [WW,WZ,ZZ,ST_top,ST_antitop]
  samplesMC   = allMCSamples
  #samplesData = [ICHEP]
  samplesData = [data]
  isData = False
  
  # load chain to list
  el_data = eventlist( samplesData, presel )
  #el_MC   = eventlist( samplesMC, presel )
  
  del samplesData, samplesMC
  
  el_data.getPileUpDist()
  #el_MC.getPileUpDist()
  
  #el_MC.doPileUpReweight(el_data.PUhist)
  
  #el_MC.doJetSmearing(False)
  el_data.doJetSmearing(False)
  
  el = el_data
  
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
