import ROOT

from loadAllSamples import *
from helpers import *
from math import *

import json

rand = ROOT.TRandom3(10**6+1)

# Define working points etc
jet_pt_split = 15
presel = 'Sum$(jet_pt>0)>1'
sigCut = 9.

etabins = [(0,0.8),(0.8,1.3),(1.3,1.9),(1.9,2.5),(2.5,100)]

#samples = [WW,WZ,ZZ]#,ST_top,ST_antitop]
samples = allMCSamples
#samples = [data]

isData = False
smear = False
if isData:
  smear = False


class event:
  def __init__(self, weight=1,jet_pt=[],jet_phi=[],jet_eta=[],jet_sigmapt=[],jet_sigmaphi=[],jet_sf=[],met_pt=0,met_phi=0,met_sumpt=0,significance=0,cov_xx=0,cov_xy=0,cov_yy=0,det=0):
    self.weight       = weight
    self.jet_pt       = jet_pt
    self.jet_phi      = jet_phi
    self.jet_eta      = jet_eta
    self.jet_sigmapt  = jet_sigmapt
    self.jet_sigmaphi = jet_sigmaphi
    self.jet_sf       = jet_sf
    self.met_pt       = met_pt
    self.met_phi      = met_phi
    self.met_sumpt    = met_sumpt
    self.sig          = significance
    self.cov_xx       = cov_xx
    self.cov_xy       = cov_xy
    self.cov_yy       = cov_yy
    self.det          = det
  def setSignif(self, significance=0,cov_xx=0,cov_xy=0,cov_yy=0,det=0):
    self.sig          = significance
    self.cov_xx       = cov_xx
    self.cov_xy       = cov_xy
    self.cov_yy       = cov_yy
    self.det          = det


#load chain to list, hope this increases speed
event_list = []
for s in samples:

  weight = s.weight
  s.chain.Draw('>>eList',presel)
  elist = ROOT.gDirectory.Get("eList")
  number_events = elist.GetN()
  print "Sample",s.name,", looping over " + str(number_events) + " events"

  #Event Loop starts here
  for i in range(number_events):
    s.chain.GetEntry(elist.GetEntry(i))
    if i>0 and (i%100000)==0:
      print "Filled ",i

    jet_pts   = [x for x in s.chain.jet_pt]
    jet_sigmapts  = [x for x in s.chain.jet_sigmapt]
    jet_phis  = [x for x in s.chain.jet_phi]
    jet_sigmaphis = [x for x in s.chain.jet_sigmaphi]
    jet_etas  = [abs(x) for x in s.chain.jet_eta]
    jet_sfs   = [x for x in s.chain.jet_sf]
    met_pt    = s.chain.met_pt
    met_phi   = s.chain.met_phi
    met_sumpt = s.chain.met_sumpt
    event_list.append(event(weight,jet_pts,jet_phis,jet_etas,jet_sigmapts,jet_sigmaphis,jet_sfs,met_pt,met_phi,met_sumpt))
    

def getSig(args, ev):
  cov_xx = 0
  cov_xy = 0
  cov_yy = 0

  dmet_x = 0
  dmet_y = 0

  for i,j in enumerate(ev.jet_pt):
    jpt = j
    if j > jet_pt_split:
      index = 0
      found = False
      for ia, a in enumerate(etabins):
        if a[0] <= abs(ev.jet_eta[i]) < a[1]:
          index = ia
          found = True

      if not found: print 'jet eta outside (0,100), sth went wrong'
      # jet smearing

      cj = cos(ev.jet_phi[i])
      sj = sin(ev.jet_phi[i])

      if smear:
        jetsf = ev.jet_sf[i]
        if( jetsf < 1 ): jetsf = 1
        sm = rand.Gaus(0, sqrt((jetsf**2)-1) * ev.jet_sigmapt[i]*jpt)
        #print sm
        dmet_x -= cj*sm
        dmet_y -= sj*sm

        jpt += sm;

      dpt = args[index] * jpt * ev.jet_sigmapt[i]
      dph =               jpt * ev.jet_sigmaphi[i]

      dpt2 = dpt**2
      dph2 = dph**2

      cov_xx += dpt2*cj*cj + dph2*sj*sj
      cov_xy += (dpt2-dph2)*cj*sj
      cov_yy += dph2*cj*cj + dpt2*sj*sj

  # pseudo-jet stuff
  cov_tt = args[5]**2 + args[6]**2*ev.met_sumpt
  cov_xx += cov_tt
  cov_yy += cov_tt

  det = cov_xx*cov_yy - cov_xy*cov_xy

  ncov_xx = cov_yy / det
  ncov_xy = -cov_xy / det
  ncov_yy = cov_xx / det

  met_x = ev.met_pt * cos(ev.met_phi)
  met_y = ev.met_pt * sin(ev.met_phi)

  if smear:
    met_x += dmet_x
    met_y += dmet_y
    met_pt = sqrt(met_x**2 + met_y**2)

  sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy
  ev.setSignif(sig,cov_xx,cov_xy,cov_yy,det)

def getLL(args):
  LL = 0.
  for ev in event_list:
    getSig(args, ev)
    LL += ev.weight * (ev.sig + log(ev.det))
  return LL


class minLL( ROOT.TPyMultiGenFunction ):
  def __init__( self ):
    ROOT.TPyMultiGenFunction.__init__( self, self )

  def NDim( self ):
    return 7

  def DoEval( self, args ):
    LL = getLL(args)
    return LL


print 'Total events:',len(event_list)

gmin = ROOT.Math.Factory.CreateMinimizer("Minuit2")
gmin.SetTolerance(10.0)
gmin.SetStrategy(0)
gmin.SetPrintLevel(3)

LL = minLL()

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
full_elist = event_list

#filter events with high significance
event_list = [x  for x in event_list if x.sig < sigCut]

print
print 'Now fitting after applying significance cut'
print 'Total events:',len(event_list)

gmin.SetStrategy(1)
gmin.Minimize()
gmin.Hesse()

#min.X()
