from random import random
import matplotlib.pyplot as plt
from basic import *

""" After looking at your figure 1 answer these questions:
 Q: Which (s,v) combinations are the most cost effective?
(16,24) is most effective with low probability p (0 <= p <= 0.01).
(8,48) is most effective with moderate probability p (0.01 < p <= 0.04).
(4,96) is most effective with high probability p (0.04 < p <= 1).

 Q: Where would the Basic and P-Best costs be on this plot?
The Basic cost would be a straight line at 384 because 384 PCR tests are conducted
no matter what the results of the tests are. The P-Best plot would be lower than
all dorfman combinations for low p values and would rise exceedingly for high p values.

After looking at your figure 2 answer these questions:
 Q: How do the mean=average number of tests needed for
    the (16,24) and the (4,96) versions compare?
(16,24) has more variability in its distribution, with an average cost of about 125 tests.
(4,24) is more densely clustered and also has an average cost of around 125 tests.

 Q: Which design is has more variability in the
    number of tests required? Guesses as to why?
The (16,24) design has more variability. There is greater variability because
there are s=16 samples per vial, so any positive PCR test on each vial increases
the cost by 4x as much as the (4,96) design.
"""

def mean( L ):
    if len(L)>0:
        return sum(L) / len(L)
    else:
        return 0

def nDorfman(p, s, v):
    """Conducts a simulation of a batch with s samples in v vials at prob p
    and returns the number of PCR tests that must be conducted
    >>> nDorfman(0, 16, 24)
    24
    >>> nDorfman(1, 16, 24)
    408
    """
    # *Simulation* of a vial is the process of counting the number of
    # positive tests are in each vial.  Any time that is > 0, the vial must
    # be re-tested.

    cost = 0
    for vial in range(v):
        cost += 1 # initial PCR test conducted
        samples = [random() < p for sample in range(s)]
        # iterate through every sample in the vial
        vialinfected = any(samples) # search for positive tests
        if vialinfected:
            cost += s # all samples in vial must be PCR tested
    return cost

def test():
    """Run document tests."""
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    test()
    Nbatch = 1000

    counts = dict()
    pvals = [0.02*(2**i) for i in range(-5,3)]
    svVals = [(16,24), (8,48), (4,96)]
    for p in pvals:
        for (s,v) in svVals:
            counts[p,s,v] = [ nDorfman(p, s, v) for _ in range(Nbatch)]
    plt.figure(1) #explicitly numbers the figures
    for (s,v) in svVals:
        plt.plot( pvals, [mean(counts[p, s, v]) for p in pvals], label="{},{}".format(s,v))
    plt.xlabel('probability p')
    plt.ylabel('Dorfman cost')
    plt.legend()
    plt.savefig('Dorfman-cost.pdf')
    plt.show()

    plt.figure(2)    #switch to second figure
    p = 0.02
    newsvVals = [(16,24), (4,96)] # comparing two methods
    for (s,v) in newsvVals:
        plt.hist( counts[p,s,v] , bins=range(0,408+1),
                      align='left', alpha=0.5, label=str( (s,v) ))
        # this histogram zooms in on Figure 1 at p = 0.02 and examines the
        # variability of the number of tests required per simulation of Nbatch
    plt.xlabel('number of tests, p={:.4f}'.format(p))
    plt.ylabel('Count')
    plt.legend()
    plt.savefig('Dorfman-histogram.pdf')
    plt.show()
