from random import random
import matplotlib.pyplot as plt

"""
(a) how many pairs of individuals have zero vials in common,
    only one vial that overlaps, two vials that overlap,
    three vials that overlap?
[33792, 25344, 14400, 0, 0, 0, 0]
33792 pairs of individuals have zero vials in common.
25344 pairs of individuals have one vial in common.
14400 pairs of individuals have two vials in common.

(b) At about what p are the average number of True Positives and
    False Positives equal for the Shental design?
At p=0.5, ATP = AFP
"""
def mean( L ):
    if len(L)>0:
        return sum(L) / len(L)
    else:
        return 0

def overlaps(pools):
    """returns a list with a count of the number of vials that
    overlaps between sample i and sample j, with i < j
    >>> p = [ {1,2,3,4,5,6}, {1,7,8,9,10,11}, {11,12,13,14,15,16} ]
    >>> overlaps(p)
    [1, 2, 0, 0, 0, 0, 0]
    >>> pD = [ {1,2,3}, {4,5,6}, {7,8,9}, {3,6,9} ]
    >>> overlaps(pD)
    [3, 3, 0, 0]
    >>> pC = [ {5,6,7,8}, {5,6,7,8}, {5,6,7,8}, {5,6,7,8} ]
    >>> overlaps(pC)
    [0, 0, 0, 0, 6]
    """

    # count initializes list of len {vials} + 1 indexes (all with zero value).
    # if patient has sample put into 4 vials, count returns list of 5 zeroes.
    count = [0]*(len(pools[0])+1)
    for i in range(len(pools)): # iterate through each sample
        for j in range(i+1, len(pools)): # iterate through adj sample
            same = len(pools[i] & pools[j])
            # & finds all equ integers in adj sets... same = len of set returned
            count[same] += 1 # len of set returned determines index incremented
    return count

def pBest(p, pools):
    """conducts one simulation and returns the number of
    (TP, FP, FN, TN)=True and False Positives and Negatives.
    >>> pools = getPools('pooling384-48-by-sample.txt')
    >>> pBest(0, pools)
    (0, 0, 0, 384)
    >>> pools = getPools('pooling384-48-by-sample.txt')
    >>> pBest(1, pools)
    (384, 0, 0, 0)
    """

    # initialize count at 0 for TP, FP, FN, and TN
    TP = 0
    FP = 0
    FN = 0
    TN = 0
    vials = [False]*48 # vials initially negative
    samples = [random() < p for _ in range(len(pools))] # true distribution
    for i in range(len(samples)):
        if samples[i]: # patient actually sick
            for vial in pools[i]:
                vials[vial] = True # change index in vials from False to True
    # vials now contains correct labeling (T or F) based on samples in vial
    broadResults = [] # sample distribution
    for i in range(len(samples)):
        vialList = []
        for vial in pools[i]:
            vialList.append(vials[vial])
        result = all(vialList) # if any vials are False, result = False
        broadResults.append(result)

    for i in range(len(samples)):
        if samples[i] and broadResults[i]:
        # true and sample distrib are true
            TP += 1
        elif not samples[i] and broadResults[i]:
        # true distrib == False and sample distrib == True
            FP += 1
        elif samples[i] and not broadResults[i]:
        # true distrib == True and sample distrib == False
            FN += 1 # can never happen
        else:
        # true and sample distrib == False
            TN += 1
    return (TP, FP, FN, TN)

def simulate(p, pools, numTrials=100):
    """conducts numTrials simulations to accumulate the total
    counts of (TP, FP, FN, TN)
    >>> pools = getPools('pooling384-48-by-sample.txt')
    >>> simulate(0, pools, 100)
    (0, 0, 0, 38400)
    >>> pools = getPools('pooling384-48-by-sample.txt')
    >>> simulate(1, pools, 100)
    (38400, 0, 0, 0)
    """

    TP = 0
    FP = 0
    FN = 0
    TN = 0
    simulations = [pBest(p, pools) for _ in range(numTrials)]
    # returns numTrials tuples (TP, FP, FN, TN)
    for i in simulations: # iterate across numTrials tuples in simulations
        TP += i[0]
        FP += i[1]
        FN += i[2]
        TN += i[3]
    return (TP, FP, FN, TN) # each index of numTrials tuples added together

def printStats(TP, FP, FN, TN):
    print('True Positive ={}'.format(TP))
    print('False Positive ={}'.format(FP))
    print('False Negative ={}'.format(FN))
    print('True Negative ={}'.format(TN))

def getPools(filename='pooling384-48-by-sample.txt'):
    """returns a list containing the sets of vials"""
    pools = [ ]
    with open(filename) as fin:
        for line in fin:
            sL = line.strip().split()
            pools.append(set([int(v) for v in sL[1:]]))
            # change sets in pools from strings to ints
    return pools

def test():
    """Run document tests."""
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    test()

    # build pools dictionary
    pools = getPools('pooling384-48-by-sample.txt')
    print('the overlap array:',overlaps(pools))
    #pvals for x axis
    pvals = [0.001*i for i in range(0,16+1)]
    #FPR values matching p in pvals
    FPR = []
    for p in pvals:
        indFPR = simulate(p, pools, numTrials=5000)
        FPR.append(indFPR[1]/sum(indFPR))
        # indFPR: FP/(TP+FP+FN+TN)
    plt.figure(1)
    plt.plot( pvals, FPR, label='PBest: 384 samples, 48 wells, 6 wells/sample')
    plt.xlabel('True Positive Rate')
    plt.ylabel('False Positive Rate')
    plt.legend()
    plt.savefig('PBest.pdf')
    plt.show()
