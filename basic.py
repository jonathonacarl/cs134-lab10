from random import random
import matplotlib.pyplot as plt

"""
(a) Question: at p=0.02=2%, would the testing center likely observe:
  0/384=0% positives? No. Not very likely
  12/384=3% positives? Yes. Somewhat likely.
  15/384=4% positives? Not likely.
(b) Question: at p=0.002=0.2%, would the testing center likely observe:
  4/384=1% positives? No. Not as likely.
  0/384=0% positives? Yes. Very likely.
"""

def basic(p, samples=384):
    """simulates 1 batch with sick probability p
    and returns the number of sick individuals in that batch"""

    posTest = 0
    for c in range(samples):
        r = random()
        if r < p:
            posTest += 1
    return posTest

    # also can be expressed as the following:
    # tests = [0,1] # negTest = 0 & posTest = 1
    # distrib = [1-p, p] # negTest have 1-p probability
    # results = (random.choices(tests, distrib, k=384))
    # random.choices takes in tests, probability mapped to each number in tests, and size of output
    # return sum(results) # since posTest = 1, sum returns total posTests

def simBasic(p, samples, numTrials=1000):
    """conducts numTrials simulations and
    returns a list of the number of infected individuals in each batch"""
    return [basic(p) for _ in range(numTrials)] # len(list) == numTrials

if __name__ == '__main__':
    samples = 384
    Nbatch = 1000
    pvals = [0.0002, 0.002, 0.02]

    plt.figure(1)
    for p in pvals:
        plt.hist( simBasic(p, samples, Nbatch), bins=range(0,25+1),
                      align='left', alpha=0.5, label=str(p))
    #This histogram's data goes in bins of width 1 from 0 to 25
    #alpha=0.5 makes data 50% transparent
    
    #align='left' centers the bins over the integer values
    plt.xlabel('number of positives')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig('Basic-histogram.pdf')
    plt.show()
