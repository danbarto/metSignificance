import ROOT

from loadAllSamples_puppi import *
#from loadAllSamples_August import *
from eventlist import *
from likelihood import *
from helpers import *
from math import *

import json

def main():

  tuneName = 'august_2016BH_njet0p_puppi_v2'
  
  # Define working points etc
  presel = 'Sum$(jet_pt>30&&abs(jet_eta)<2.5&&jet_passid)>=0'
  sigCut = 9.
  
  #samplesMC = [WW_M17,WZ_M17,ZZ_M17,ST_top_M17,ST_antitop_M17]
  #samplesMC = [TTJets_M17]
  samplesMC   = allMCSamples
  #samplesMC   = [DY]
  #samplesData = [ICHEP]
  samplesData = [data]
  tightZwindow = False
  
  # load chain to list
  el_data = eventlist( samplesData, presel, isData=True, tiny=True )
  el_MC   = eventlist( samplesMC, presel, isData=False, tiny=True )
  
  del samplesData, samplesMC
  
  el_data.getPileUpDist()
  el_MC.getPileUpDist()
  
  el_MC.doPileUpReweight(el_data.PUhist)
  #el_MC.doSmearing(useRand=False)
  #el_MC.doSmearing()
  
  samples = {'MC':el_MC, 'data':el_data}
  #samples = {'data':el_data}
  #samples = {'MC':el_MC}
  for sample in samples:
    print
    print 'Tuning',sample
    el = samples[sample]
    outfile = 'data/'+sample+'_'+tuneName
    outfileUnc = 'data/'+sample+'_Uncertainty_'+tuneName
    
    # Protection
    el.evlist = [x for x in el.evlist if not math.isnan(x.met_pt)]
    
    # Use additional cut
    if tightZwindow:
      for x in el.evlist: x.getMuonInvMass()
      el.evlist = [x for x in el.evlist if 80 < x.muonInvMass < 100]
    
    # Do the minimization
    print 'Total events:',len(el.evlist)
    
    gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
    gmin.SetTolerance(10.0)
    gmin.SetStrategy(0)
    gmin.SetPrintLevel(3)
    
    LL = minLL(el)
    
    gmin.SetFunction(LL)
    
    variable  = ['a1','a2','a3','a4','a5','u1','u2']
    #variable  = ['a1','a2','a3','a4','a5','N1','S1']
    step      = [0.05]*7
    start     = [1.0,  1.0,  1.0,  1.0, 1.0, 0., .5]
    #start     = [1.0]*7
    
    print 'Minimizing parameters',variable
    print 'With stepsize of', step
    print 'And starting values', start
    
    
    for i in range(0,7):
      gmin.SetVariable(i,variable[i],start[i], step[i])
    
    gmin.Minimize()
    
    ##save full event list
    #full_el = el
    
    #filter events with high significance
    el.evlist = [x  for x in el.evlist if x.sig < sigCut]
    
    print
    print 'Now fitting after applying significance cut'
    print 'Total events:',len(el.evlist)
    
    gmin.SetStrategy(1)
    gmin.Minimize()
    gmin.Hesse()
    
    pars = [gmin.X()[i] for i in range(0,7)]
    uncs = [gmin.Errors()[i] for i in range(0,7)]
    
    with open(outfile, 'w') as of:
      json.dump(pars, of)
    with open(outfileUnc, 'w') as of:
      json.dump(uncs, of)
    
    del el
