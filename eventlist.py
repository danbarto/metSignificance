import ROOT
from math import *
rand = ROOT.TRandom3(10**6+1)

etabins = [(0,0.8),(0.8,1.3),(1.3,1.9),(1.9,2.5),(2.5,100)]

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
    self.dmet_x       = 0
    self.dmet_y       = 0

  def setSignif(self, significance=0,cov_xx=0,cov_xy=0,cov_yy=0,det=0):
    self.sig          = significance
    self.cov_xx       = cov_xx
    self.cov_xy       = cov_xy
    self.cov_yy       = cov_yy
    self.det          = det
  
  def getSig(self, args, smear):
    self.sig          = 0
    self.cov_xx       = 0
    self.cov_xy       = 0
    self.cov_yy       = 0
    self.det          = 0
    self.dmet_x       = 0
    self.dmet_y       = 0
    jet_pt_split = 15.
    for i,j in enumerate(self.jet_pt):
      jpt = j
      if j > jet_pt_split:
        index = 0
        found = False
        for ia, a in enumerate(etabins):
          if a[0] <= abs(self.jet_eta[i]) < a[1]:
            index = ia
            found = True

        if not found: print 'jet eta outside (0,100), sth went wrong'
        # jet smearing

        cj = cos(self.jet_phi[i])
        sj = sin(self.jet_phi[i])

        if smear:
          jetsf = self.jet_sf[i]
          if( jetsf < 1 ): jetsf = 1
          sm = rand.Gaus(0, sqrt((jetsf**2)-1) * self.jet_sigmapt[i]*jpt)
          #print sm
          self.dmet_x -= cj*sm
          self.dmet_y -= sj*sm

          jpt += sm;

        dpt = args[index] * jpt * self.jet_sigmapt[i]
        dph =               jpt * self.jet_sigmaphi[i]

        dpt2 = dpt**2
        dph2 = dph**2

        self.cov_xx += dpt2*cj*cj + dph2*sj*sj
        self.cov_xy += (dpt2-dph2)*cj*sj
        self.cov_yy += dph2*cj*cj + dpt2*sj*sj
      else:
        print 'not doing anything'
    # pseudo-jet stuff
    cov_tt = args[5]**2 + args[6]**2*self.met_sumpt
    self.cov_xx += cov_tt
    self.cov_yy += cov_tt

    self.det = self.cov_xx*self.cov_yy - self.cov_xy*self.cov_xy

    ncov_xx = self.cov_yy / self.det
    ncov_xy = -self.cov_xy / self.det
    ncov_yy = self.cov_xx / self.det

    met_x = self.met_pt * cos(self.met_phi)
    met_y = self.met_pt * sin(self.met_phi)

    if smear:
      met_x += self.dmet_x
      met_y += self.dmet_y
      self.met_pt = sqrt(met_x**2 + met_y**2)

    self.sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy



class eventlist:
  def __init__(self, samples, cut):
    self.evlist   = []
    self.samples  = samples
    self.cut      = cut
    self.smear    = False

    for s in self.samples:
    
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
        jet_sfs   = [x for x in s.chain.jet_sf]
        met_pt    = s.chain.met_pt
        met_phi   = s.chain.met_phi
        met_sumpt = s.chain.met_sumpt
        self.evlist.append(event(weight,jet_pts,jet_phis,jet_etas,jet_sigmapts,jet_sigmaphis,jet_sfs,met_pt,met_phi,met_sumpt))
  
  def doJetSmearing(self,smear):
    self.smear = smear
  
  def getLL(self, args):
    self.LL = 0.
    for ev in self.evlist:
      ev.getSig(args, self.smear)
      self.LL += ev.weight * (ev.sig + log(ev.det))
    return self.LL
  
