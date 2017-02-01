import ROOT

from eventlist import *
from likelihood import *
from helpers import *
from math import *

import json

import pickle
import scipy
import scipy.optimize
import numpy

samples = pickle.load(file('test_dump_list_pkl'))
el = samples['MC']
LL = minLL(el)

start = numpy.array([1.,1.,1.,1.,1.,0.,0.5])
#el.getLL(start)

res = scipy.optimize.minimize(LL,start,method='Powell',tol=0.001)
