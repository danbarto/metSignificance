import math
import itertools
import random
import numpy

def binomial(n, k):
  assert n is not type(1), "n must be integer"
  assert k is not type(1), "k must be integer"
  return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

#itertools.combinations([1,2,3,4,5],3)
def jackknife(s, d):

    N_s = len(s)
    N_jk = binomial(N_s,d)
    #print N_jk
    y_bar = 0

    for it in itertools.combinations(s,d):
        y_bar += sum(it)
    y_bar = y_bar/N_jk

    C_ii = 0
    for it in itertools.combinations(s,d):
        C_ii += (sum(it) - y_bar)**2

    C_ii = (N_s-d)*C_ii/(d*N_jk)
    return C_ii

def jackknifeMultiDim(s, d):
    '''Jackknife delete-d method. Takes list with n-dim tuples as elements, deletes subsets of size d. Returns estimated covariance matrix as numpy array
    '''
    dim = len(s[0])
    N_s = len(s)
    N_jk = binomial(N_s,d)
    y_bar = [0]*dim

    for it in itertools.combinations(s,d):
        # add the components of the vectors within 'it', then sum up the total
        y = map(sum,itertools.izip(*it))
        y_bar = map(sum,itertools.izip(y,y_bar))

    # average the sums
    y_bar = [ y_b/N_jk for y_b in y_bar ]

    cov = numpy.zeros((dim,dim))
    n = range(0,dim)
    for ind in itertools.combinations_with_replacement(n, 2):
        #print ind[0], ind[1]
        for it in itertools.combinations(s,d):
            y  = map(sum,itertools.izip(*it)) # y_i, y_j,...
            cov[ind[0]][ind[1]] += (y[ind[0]] - y_bar[ind[0]])*(y[ind[1]] - y_bar[ind[1]])
            if ind[0] != ind[1]: cov[ind[1]][ind[0]] += cov[ind[0]][ind[1]]

    cov = cov*(N_s-d)/(d*N_jk)
    return cov


def nth(iterable, n, default=None):
    "Returns the nth item or a default value"
    return next(itertools.islice(iterable, n, None), default)

def randomJackknife(s,d,Nmax):

    N_s = len(s)
    N_jk = binomial(N_s,d)
    y_bar = 0
    rseeds = [random.randrange(0,N_jk) for i in range(Nmax)]
    rseeds = sorted(rseeds)

    for r in rseeds:
        comb = itertools.combinations(s,d)
        y_bar += sum(nth(comb,r))
    y_bar = y_bar/Nmax

    C_ii = 0
    for r in rseeds:
        comb = itertools.combinations(s,d)
        C_ii += (sum(nth(comb,r)) - y_bar)**2

    C_ii = (N_s-d)*C_ii/(d*Nmax)
    return C_ii

