#include <TH1F.h>
#include <TFile.h>
#include <TDirectory.h>

TFile *f = new TFile("data/nvert_2016BCDEFG_0jet.root");
TH1F * h = (TH1F*)f->Get("data")->Clone("h"); 

double nvertReweight(double val = 0)
{
   if (!h) 
       return 0.;
   return h->GetBinContent(h->FindBin(val));     
}
