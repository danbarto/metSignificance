import ROOT
import math
import gc

rand = ROOT.TRandom3(10**6+1)

etabins = [(0,0.8),(0.8,1.3),(1.3,1.9),(1.9,2.5),(2.5,100)]
etabins = [0.8,1.3,1.9,2.5,100]

def getBin(abseta):
  for i, a in enumerate(etabins):
    if abseta < a:
      return i
      break

class event:
  def __init__(self, weight=1,jet_pt=[],jet_phi=[],jet_eta=[],jet_etabin=[],jet_sigmapt=[],jet_sigmaphi=[],jet_sf=[],met_pt=0,met_phi=0,met_sumpt=0,nvert=0,significance=0,cov_xx=0,cov_xy=0,cov_yy=0,det=0,group=None,isData=False):
    self.weight       = weight
    self.jet_pt       = jet_pt
    self.jet_phi      = jet_phi
    self.jet_eta      = jet_eta
    self.jet_etabin   = jet_etabin
    self.jet_sigmapt  = jet_sigmapt
    self.jet_sigmaphi = jet_sigmaphi
    self.jet_sf       = jet_sf
    self.met_pt       = met_pt
    self.met_phi      = met_phi
    self.met_sumpt    = met_sumpt
    self.nvert        = nvert
    self.sig          = significance
    self.cov_xx       = cov_xx
    self.cov_xy       = cov_xy
    self.cov_yy       = cov_yy
    self.det          = det
    #self.dmet_x       = 0
    #self.dmet_y       = 0
    self.group        = group
    self.isData       = isData

  def setSignif(self, significance=0,cov_xx=0,cov_xy=0,cov_yy=0,det=0):
    self.sig          = significance
    self.cov_xx       = cov_xx
    self.cov_xy       = cov_xy
    self.cov_yy       = cov_yy
    self.det          = det
  
  def getSig(self, args, smear):
    cov_xx       = 0
    cov_xy       = 0
    cov_yy       = 0
    dmet_x       = 0
    dmet_y       = 0
    jet_pt = self.jet_pt
    
    i = 0
    for j in jet_pt:
      j_pt = j
      j_phi = self.jet_phi[i]
      j_sigmapt = self.jet_sigmapt[i]
      j_sigmaphi = self.jet_sigmaphi[i]
      j_sf = self.jet_sf[i]
      #j_eta = abs(self.jet_eta[i])
      index = self.jet_etabin[i]
      #index = 0
      #ia = 0
      #found = False
      #for a in etabins:
      #  if j_eta < a:
      #    index = ia
      #    found = True
      #    break
      #  ia += 1
      #  #if a[0] <= j_eta < a[1]:
      #  #  index = ia
      #  #  found = True

      #if not found: print 'jet eta outside (0,100), sth went wrong'
      ## jet smearing
      
      cj = math.cos(j_phi)
      sj = math.sin(j_phi)

      if smear:
        if( j_sf < 1 ): j_sf = 1
        rd = rand.Gaus(0, math.sqrt((j_sf*j_sf)-1))
        #print rd
        sm = rd * j_sigmapt*j_pt
        dmet_x -= cj*sm
        dmet_y -= sj*sm

        j_pt += sm;
      
      if not self.isData:
        dpt = args[index] * j_pt * j_sigmapt# * j_sf
      else:
        dpt = args[index] * j_pt * j_sigmapt
      dph =               j_pt * j_sigmaphi

      dpt *= dpt
      dph *= dph

      cov_xx += dpt*cj*cj + dph*sj*sj
      cov_xy += (dpt-dph)*cj*sj
      cov_yy += dph*cj*cj + dpt*sj*sj
      
      i += 1

    # pseudo-jet stuff
    cov_tt = args[5]*args[5] + args[6]*args[6]*self.met_sumpt
    cov_xx += cov_tt
    cov_yy += cov_tt

    det = cov_xx*cov_yy - cov_xy*cov_xy

    ncov_xx =  cov_yy / det
    ncov_xy = -cov_xy / det
    ncov_yy =  cov_xx / det

    met_x = self.met_pt * math.cos(self.met_phi)
    met_y = self.met_pt * math.sin(self.met_phi)

    if smear:
      met_x += dmet_x
      met_y += dmet_y
      self.met_pt = math.sqrt(met_x*met_x + met_y*met_y)

    self.cov_xx = cov_xx
    self.cov_xy = cov_xy
    self.cov_yy = cov_yy
    self.det = det
    self.sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy



class eventlist:
  def __init__(self, samples, cut):
    self.evlist   = []
    #self.samples  = samples
    self.cut      = cut
    self.smear    = False
    self.args     = []

    for s in samples:
    
      weight = s.weight
      s.chain.Draw('>>eList',self.cut)
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
        jet_etabins = [getBin(abs(x)) for x in s.chain.jet_eta]
        jet_sfs   = [x for x in s.chain.jet_sf]
        met_pt    = s.chain.met_pt
        met_phi   = s.chain.met_phi
        met_sumpt = s.chain.met_sumpt
        gc.disable()
        self.evlist.append(event(weight,jet_pts,jet_phis,jet_etas,jet_etabins,jet_sigmapts,jet_sigmaphis,jet_sfs,met_pt,met_phi,met_sumpt,s.chain.nvertices,group=s.subGroup,isData=s.isData))
        gc.enable()
  
  def doJetSmearing(self,smear):
    self.smear = smear
  
  def getPileUpDist(self,nvert_max=100):
    self.PUhist = []
    self.nvert_max = nvert_max
    for i in range(0,nvert_max):
      self.PUhist.append(0)
    for ev in self.evlist:
      self.PUhist[ev.nvert] += ev.weight
    norm = sum(self.PUhist)
    self.PUhist = [i/norm for i in self.PUhist]
  
  def doPileUpReweight(self, PUhist_other):
    if len(PUhist_other) != len(self.PUhist):
      print 'PU histograms do have different number of bins, aborting'
      return
    ratio = [PUhist_other[i]/self.PUhist[i] if self.PUhist[i]>0 else 0 for i in range(0,self.nvert_max)]
    minx = min([r for r in ratio if r>0])
    for ev in self.evlist:
      PUweight = ratio[ev.nvert]
      if PUweight <= 0: PUweight = minx
      ev.weight *= PUweight
  
  def getLL(self, args):
    self.LL   = 0.
    self.args = args
    for ev in self.evlist:
      ev.getSig(args, self.smear)
      self.LL += ev.weight * (ev.sig + math.log(ev.det))
    return self.LL
  
