import ROOT

from loadAllSamples import *

import json



print 'Loaded Parameters for MC from file'
with open('paraMC.txt', 'r') as paraMC:
  par = json.load(paraMC)
  print(json.dumps(par,sort_keys=True,indent=2))

