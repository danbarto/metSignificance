import ROOT
import math
from metSignificance.tools.jackknife import *
import random

rand = ROOT.TRandom3(10**6+1)

etabins = [(0,0.8),(0.8,1.3),(1.3,1.9),(1.9,2.5),(2.5,100)]
etabins = [0.8,1.3,1.9,2.5,100]

def getBin(abseta):
  for i, a in enumerate(etabins):
    if abseta < a:
      return int(i)
      break

def cartesian(pt, phi):
    return (pt*math.cos(phi), pt*math.sin(phi))

class event:
  #def __init__(self, weight, sample, i, elist, isData, tiny=False, unc_cov=True):
  #  sample.chain.GetEntry(elist.GetEntry(i))
  #  self.weight       = weight
  #  self.jet_pt       = [x for x in sample.chain.jet_pt]
  #  self.jet_sigmapt  = [x for x in sample.chain.jet_sigmapt]
  #  self.jet_phi      = [x for x in sample.chain.jet_phi]
  #  self.jet_sigmaphi = [x for x in sample.chain.jet_sigmaphi]
  #  self.jet_etabin   = [getBin(abs(x)) for x in sample.chain.jet_eta]

  #  if not isData:
  #    self.jet_sf     = [x for x in sample.chain.jet_sf]

  #  self.met_pt       = sample.chain.met_pt
  #  self.met_phi      = sample.chain.met_phi
  #  self.met_sumpt    = sample.chain.met_sumpt
  #  self.nvert        = sample.chain.nvertices
  #  self.sig          = 0.
  #  self.det          = 0.

  #  if unc_cov:
  #      cands_vec = []
  #      for j,c in enumerate(sample.chain.cand_pt):
  #          if c>0:
  #              cands_vec.append(cartesian(sample.chain.cand_pt[j], sample.chain.cand_phi[j]))
  #      cov = jackknifeMultiDim(cands_vec,1)
  #      self.unc_cov_xx     = cov[0][0]
  #      self.unc_cov_xy     = cov[0][1]
  #      self.unc_cov_yy     = cov[1][1]
  #      print len(cands_vec)
  #      del cands_vec, cov
  #      print self.unc_cov_xx,self.unc_cov_xy,self.unc_cov_yy

  #  if not tiny:
  #    self.group        = sample.subGroup
  #    self.muon_pt      = [x for x in sample.chain.muon_pt]
  #    self.muon_eta     = [x for x in sample.chain.muon_eta]
  #    self.muon_phi     = [x for x in sample.chain.muon_phi]
  
  def __init__(self, weight, jet_pt,jet_sigmapt,jet_phi,jet_sigmaphi,jet_etabin,jet_sf,met_pt,met_phi,met_sumpt,nvertices,c_xx,c_xy,c_yy):
    self.weight       = weight
    self.jet_pt       = jet_pt
    self.jet_sigmapt  = jet_sigmapt
    self.jet_phi      = jet_phi
    self.jet_sigmaphi = jet_sigmaphi
    self.jet_etabin   = jet_etabin
    self.jet_sf       = jet_sf 

    self.met_pt       = met_pt
    self.met_phi      = met_phi
    self.met_sumpt    = met_sumpt
    self.nvert        = nvertices
    self.sig          = 0.
    self.det          = 0.

    self.unc_cov_xx     = c_xx#cov[0][0]
    self.unc_cov_xy     = c_xy#cov[0][1]
    self.unc_cov_yy     = c_yy#cov[1][1]

  
  def getMuonInvMass(self):
    self.muonInvMass = math.sqrt(2*self.muon_pt[0]*self.muon_pt[1]*(math.cosh(self.muon_eta[0]-self.muon_eta[1])-math.cos(self.muon_phi[0]-self.muon_phi[1])))
  
  def smearJets(self, useRand=True):
    dmet_x       = 0
    dmet_y       = 0
    jet_pt = self.jet_pt
    for i,j in enumerate(jet_pt):
      j_pt = j
      j_sf = self.jet_sf[i]
      j_phi = self.jet_phi[i]
      j_sigmapt = self.jet_sigmapt[i]
      if( j_sf < 1 ): j_sf = 1
      if useRand: rd = rand.Gaus(0, math.sqrt((j_sf*j_sf)-1))
      else: rd = 0.
      sm = rd * j_sigmapt*j_pt
      cj = math.cos(j_phi)
      sj = math.sin(j_phi)
      dmet_x -= cj*sm
      dmet_y -= sj*sm
      self.jet_pt[i] += sm
      self.jet_sigmapt[i] = j_sigmapt * j_sf

    del self.jet_sf
    met_x = self.met_pt * math.cos(self.met_phi)
    met_y = self.met_pt * math.sin(self.met_phi)
    met_x += dmet_x
    met_y += dmet_y

    self.met_pt = math.sqrt(met_x*met_x + met_y*met_y)
    self.met_phi = (met_y/abs(met_y)) * math.acos(met_x/self.met_pt)
    
  def getSig(self, args):
    cov_xx       = 0
    cov_xy       = 0
    cov_yy       = 0
    jet_pt = self.jet_pt
    i = 0
    for j in jet_pt:
      j_pt = j
      j_phi = self.jet_phi[i]
      j_sigmapt = self.jet_sigmapt[i]
      j_sigmaphi = self.jet_sigmaphi[i]
      index = self.jet_etabin[i]

      cj = math.cos(j_phi)
      sj = math.sin(j_phi)
      dpt = args[index] * j_pt * j_sigmapt
      dph =               j_pt * j_sigmaphi

      dpt *= dpt
      dph *= dph

      cov_xx += dpt*cj*cj + dph*sj*sj
      cov_xy += (dpt-dph)*cj*sj
      cov_yy += dph*cj*cj + dpt*sj*sj
      
      i += 1

    # unclustered energy
    cov_tt = args[5]*args[5] + args[6]*args[6]*self.met_sumpt
    cov_xx += cov_tt
    cov_yy += cov_tt

    det = cov_xx*cov_yy - cov_xy*cov_xy

    if det>0:
        ncov_xx =  cov_yy / det
        ncov_xy = -cov_xy / det
        ncov_yy =  cov_xx / det
    else:
        print cov_xx, cov_yy, cov_xy
        ncov_xx = cov_xx if cov_xx > 0 else 1
        ncov_yy = cov_yy if cov_yy > 0 else 1
        ncov_xy = cov_xy if cov_xy > 0 else 1

    met_x = self.met_pt * math.cos(self.met_phi)
    met_y = self.met_pt * math.sin(self.met_phi)

    self.det = det
    self.sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy

  def getSigAlt(self, args):
    cov_xx       = 0
    cov_xy       = 0
    cov_yy       = 0
    jet_pt = self.jet_pt
    i = 0
    for j in jet_pt:
      j_pt = j
      j_phi = self.jet_phi[i]
      j_sigmapt = self.jet_sigmapt[i]
      j_sigmaphi = self.jet_sigmaphi[i]
      index = self.jet_etabin[i]

      cj = math.cos(j_phi)
      sj = math.sin(j_phi)
      dpt = args[index] * j_pt * j_sigmapt
      dph =               j_pt * j_sigmaphi

      dpt *= dpt
      dph *= dph

      cov_xx += dpt*cj*cj + dph*sj*sj
      cov_xy += (dpt-dph)*cj*sj
      cov_yy += dph*cj*cj + dpt*sj*sj
      #print cov_xx, cov_xy, cov_yy
      
      i += 1

    # unclustered energy
    cov_xx += args[5]**2 * self.unc_cov_xx
    cov_xy += args[5]*args[6] * self.unc_cov_xy
    cov_yy += args[6]**2 * self.unc_cov_yy
    #cov_xy += args[6]**2 * self.unc_cov_xy
    #cov_yy += args[7]**2 * self.unc_cov_yy

    det = cov_xx*cov_yy - cov_xy*cov_xy
    #print cov_xx, cov_xy, cov_yy
    
    ncov_xx =  cov_yy / det
    ncov_xy = -cov_xy / det
    ncov_yy =  cov_xx / det

    met_x = self.met_pt * math.cos(self.met_phi)
    met_y = self.met_pt * math.sin(self.met_phi)

    self.det = det
    self.sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy



class eventlist:
  def __init__(self, samples, cut, isData, tiny=False, reduction=True, reductionFactor=1e6):
    self.evlist   = []
    self.cut      = cut
    self.args     = []

    for s in samples:
      total_events = s.chain.GetEntries()
      weight = s.weight
      s.chain.Draw('>>eList',self.cut)
      elist = ROOT.gDirectory.Get("eList")
      number_events = elist.GetN()
      print "Sample",s.name,", looping over " + str(number_events) + " events"
      if reduction:
        redFactor = int(number_events/reductionFactor)
        if redFactor ==0: redFactor = 1 
      else:
        redFactor = 1
      if isData: print "will reduce sample size by", redFactor
      #Event Loop starts here
      #temp = (event(weight,s,i,elist) for i in range(number_events))
      #temp = [event(weight,s,i,elist, isData, tiny=tiny) for i in range(number_events) if (i%int(redFactor)==0 or not isData) ]
      for i in range(number_events):
        if (i%int(redFactor)==0 or not isData):
            if i%10000==0:
                print "Done with %i, corresponds to entry %i"%(i,elist.GetEntry(i))
            s.chain.GetEntry(elist.GetEntry(i))
            #s.chain.GetEntry(i)
            #rint = random.randint(1,total_events)
            #print rint
            #s.chain.GetEntry(rint)
            #print elist.GetEntry(i)
            jet_pt       = [x for x in s.chain.jet_pt]
            jet_sigmapt  = [x for x in s.chain.jet_sigmapt]
            jet_phi      = [x for x in s.chain.jet_phi]
            jet_sigmaphi = [x for x in s.chain.jet_sigmaphi]
            jet_etabin   = [getBin(abs(x)) for x in s.chain.jet_eta]

            if not isData:
              jet_sf     = [x for x in s.chain.jet_sf]
            else: jet_sf = 1

            met_pt       = s.chain.met_pt
            met_phi      = s.chain.met_phi
            met_sumpt    = s.chain.met_sumpt
            nvert        = s.chain.nvertices

            #cands_vec = []
            #for j,c in enumerate(s.chain.cand_pt):
            #    if c>0:
            #        cands_vec.append(cartesian(s.chain.cand_pt[j], s.chain.cand_phi[j]))
            #cov = jackknifeMultiDim(cands_vec,1)
            #print cov[0][0]
            #del cands_vec
            #if not tiny:
            #  group        = s.subGroup
            #  muon_pt      = [x for x in s.chain.muon_pt]
            #  muon_eta     = [x for x in s.chain.muon_eta]
            #  muon_phi     = [x for x in s.chain.muon_phi]
            c_xx = s.chain.cov_xx
            c_xy = s.chain.cov_xy
            c_yy = s.chain.cov_yy
            self.evlist.append(event(weight, jet_pt,jet_sigmapt,jet_phi,jet_sigmaphi,jet_etabin,jet_sf,met_pt,met_phi,met_sumpt,nvert,c_xx, c_xy, c_yy))

      #self.evlist = self.evlist + temp
      #del temp
  
  def getPileUpDist(self,nvert_max=100):
    self.PUhist = []
    self.nvert_max = nvert_max
    for i in range(0,nvert_max):
      self.PUhist.append(0)
    for ev in self.evlist:
      if ev.nvert < nvert_max:
        self.PUhist[ev.nvert] += ev.weight
      else: self.PUhist[-1] += ev.weight
    norm = sum(self.PUhist)
    self.PUhist = [i/norm for i in self.PUhist]
  
  def doPileUpReweight(self, PUhist_other):
    if len(PUhist_other) != len(self.PUhist):
      print 'PU histograms do have different number of bins, aborting'
      return
    ratio = [PUhist_other[i]/self.PUhist[i] if self.PUhist[i]>0 else 0 for i in range(0,self.nvert_max)]
    minx = min([r for r in ratio if r>0])
    for ev in self.evlist:
      if ev.nvert < self.nvert_max: PUweight = ratio[ev.nvert]
      else: PUweight = ratio[-1]
      if PUweight <= 0: PUweight = minx
      ev.weight *= PUweight
  
  def doSmearing(self, useRand=True):
    for x in self.evlist: x.smearJets( useRand = useRand )
    if useRand: print 'Smeared jets'
    else: print 'Applied JER scale factors, no smearing done!'
  
  def getLL(self, args):
    self.LL   = 0.
    self.args = args
    for ev in self.evlist:
      #ev.getSigAlt(args)
      ev.getSig(args)
      try:
        math.log(ev.det)
      except:
        #print args[5],args[6],args[7]
        #print ev.det
        #print ev.unc_cov_xx
        #print ev.unc_cov_yy
        #print ev.unc_cov_xy
        #print ev.met_pt
        ev.det = 0.001
        
      self.LL += ev.weight * (ev.sig + math.log(ev.det))
    return self.LL
  
