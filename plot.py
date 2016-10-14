import ROOT

from loadAllSamples import *
from helpers import *
from math import *

import json



# Load parameters for significance calculation

print 'Loaded Parameters for MC from file'
with open('data/paraMC.txt', 'r') as paraMC:
  parMC = json.load(paraMC, encoding="utf-8")
  #print(json.dumps(parMC,sort_keys=True,indent=2))

#print
print 'Loaded Parameters for MC from file'
with open('data/paraData.txt', 'r') as paraData:
  parData = json.load(paraData, encoding="utf-8")
  #print(json.dumps(parData,sort_keys=True,indent=2))

print

types = ['Zmumu','top','EWK','data']

# Define working points etc
jet_pt_split = 15

# Define histograms
sig_hist = {}
for t in types:
  sig_hist[t] = ROOT.TH1F('sig_'+t,t,50,0,100)

# Calculate Significance, fill histograms

presel = 'Sum$(jet_pt>0)>1'

#allMCSamples = [allMCSamples[2]]

rand = ROOT.TRandom3(101)

# create histogram with parameters out of json
nBins = 0
bins = []
for a in parMC:
  if 'a' in a:
    nBins += 1
    bins.append(parMC[a]['eta'][1])
bins.append(0.)
bins = sorted(bins)

allSamples = allMCSamples + [data]

for s in allMCSamples:
  
  weight = s.weight
  
  s.chain.Draw('>>eList',presel)
  elist = ROOT.gDirectory.Get("eList")
  number_events = elist.GetN()
  print "Sample",s.name,", looping over " + str(number_events) + " events"

  #Event Loop starts here
  first = True

  for i in range(number_events):

    s.chain.GetEntry(elist.GetEntry(i))
    if i>0 and (i%100000)==0:
      print "Filled ",i
    
    cov_xx = 0
    cov_xy = 0
    cov_yy = 0

    dmet_x = 0
    dmet_y = 0

    jet_pts   = [x for x in s.chain.jet_pt]
    jet_dpts  = [x for x in s.chain.jet_sigmapt]
    jet_phis  = [x for x in s.chain.jet_phi]
    jet_dphis = [x for x in s.chain.jet_sigmaphi]
    jet_etas  = [abs(x) for x in s.chain.jet_eta]
    jet_sfs   = [x for x in s.chain.jet_sf]
    met_pt    = s.chain.met_pt
    met_phi   = s.chain.met_phi

    for i,j in enumerate(jet_pts):
      jpt = j
      if j > jet_pt_split:
        for a in parMC:
          if parMC[a]['eta'][0] <= jet_etas[i] < parMC[a]['eta'][1]:
            parA = parMC[a]['v']
        # jet smearing
        
        cj = cos(jet_phis[i])
        sj = sin(jet_phis[i])
        
        if not s.isData:
          jetsf = jet_sfs[i]
          if( jetsf < 1 ): jetsf = 1
          sm = rand.Gaus(0, sqrt(jetsf**2-1) * jet_dpts[i]*j);
          dmet_x -= cj*sm;
          dmet_y -= sj*sm;

          jpt += sm;
        
        dpt = parA * jpt * jet_dpts[i]
        dph =        jpt * jet_dphis[i]
        
        dpt2 = dpt**2
        dph2 = dph**2
        
        cov_xx += dpt2*cj*cj + dph2*sj*sj
        cov_xy += (dpt2-dph2)*cj*sj
        cov_yy += dph2*cj*cj + dpt2*sj*sj

    # pseudo-jet stuff
    cov_tt = parMC['N']['v']**2 + parMC['S']['v']*s.chain.met_sumpt
    cov_xx += cov_tt
    cov_yy += cov_tt
    
    det = cov_xx*cov_yy - cov_xy*cov_xy

    ncov_xx = cov_yy / det
    ncov_xy = -cov_xy / det
    ncov_yy = cov_xx / det
    
    met_x = met_pt * cos(met_phi)
    met_y = met_pt * sin(met_phi)
    
    if not s.isData:
      met_x += dmet_x
      met_y += dmet_y
      met_pt = sqrt(met_x**2 + met_y**2)
    
    sig = met_x*met_x*ncov_xx + 2*met_x*met_y*ncov_xy + met_y*met_y*ncov_yy

    sig_hist[s.subGroup].Fill(sig, weight)
            

# Save results

can = ROOT.TCanvas('can','',700,700)
can.SetLogy()
sig_stack = ROOT.THStack()

sig_hist['EWK'].SetFillColor(ROOT.kRed-9)
sig_hist['top'].SetFillColor(ROOT.kYellow-9)
sig_hist['Zmumu'].SetFillColor(ROOT.kBlue-9)

sig_stack.Add(sig_hist['top'])
sig_stack.Add(sig_hist['EWK'])
sig_stack.Add(sig_hist['Zmumu'])

sig_stack.Draw('hist')




