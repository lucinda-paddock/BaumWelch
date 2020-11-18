from forwardBackward import *

WITHINTIMESLARGER = 50
RANDOMBETWEEN = 100
RANDOMWITHIN  = RANDOMBETWEEN * WITHINTIMESLARGER
RANDOMTOT = RANDOMWITHIN + RANDOMBETWEEN

def initTransD():
    ''' Initializes a transition probability matrix. '''
    transD = {}
    for nuc1 in ["A","C","G","T","a","c","g","t"]:
        # multinomial gives a dist of numbers that add up to the first parameter
        outside = list(np.random.multinomial(RANDOMBETWEEN, np.ones(4)/4, 1)[0])
        inside = list(np.random.multinomial(RANDOMWITHIN, np.ones(4)/4, 1)[0])
        # converts numbers from multinom to probabilities (add to 1 for each starting nuc)
        outsideProb = list(map(lambda x: x/RANDOMTOT, outside))
        insideProb = list(map(lambda x: x/RANDOMTOT, inside))
        outsideI = 0
        insideI = 0
        # assign probabilites to each transition pair
        for nuc2 in ["A","C","G","T","a","c","g","t"]:
            if nuc1.isupper():
                if nuc2.isupper():
                    transD[nuc1 + nuc2] = math.log(insideProb[insideI])
                    insideI += 1
                else:
                    transD[nuc1 + nuc2] = math.log(outsideProb[outsideI])
                    outsideI += 1
            else:
                if nuc2.isupper():
                    transD[nuc1 + nuc2] = math.log(outsideProb[outsideI])
                    outsideI += 1
                else:
                    transD[nuc1 + nuc2] = math.log(insideProb[insideI])
                    insideI += 1
    return transD



def makeTransIndexD(seq):
    '''Given a sequnce, makes a dictionary where the keys are all
    possible trnasitions and the values are all indexes where that 
    transition can occur. This allows us to only loop through the whole
    sequnce one, improving efficiency.'''
    transIndexD = {}
    for i in range(len(seq)-1): #all possible nucleotide pairs
        nuc1 = seq[i]
        nuc2 = seq[i+1]
        nucPairList = [nuc1+nuc2, nuc1+nuc2.lower(), nuc1.lower()+nuc2, nuc1.lower()+nuc2.lower()]
        for nucPair in nucPairList: #all possible transitions
            if not nucPair in transIndexD:
                transIndexD[nucPair] = (i,)
            else: 
                transIndexD[nucPair] = transIndexD[nucPair] + (i,)
    return transIndexD

def makeCountD(seq,transD,FTable):
    '''Given a matrix of transition frequencies and a sequence, uses
    the forward and backward algorithms to construct a matrix of expected
    counts according to the Baum-Welch algorithm'''
    countD = {}
    BTable = backward(seq, transD)
    transIndexD = makeTransIndexD(seq)
    PofX = np.logaddexp(FTable[0][-1],FTable[1][-1])
    # call helper to get counts for each transition pair
    for key,value in transD.items():
        countD[key] = calcCount(key, value, FTable, BTable, transIndexD, PofX)
    return countD

def calcCount(nucPair, transProb, FTable, BTable, transIndexD, PofX):
    '''For a given transtion (like Ag for example) calculate the expected
    count'''
    if not nucPair in transIndexD:
        count = 0
    else:
        sumTrans = -float('inf') #log of 0
        transIndexTup = transIndexD[nucPair]
        for i in transIndexTup:
            nuc1 = nucPair[0]
            nuc2 = nucPair[1]
            # get forward entry
            if nuc1.islower():
                forward = FTable[1][i]
            else: 
                forward = FTable[0][i]
            # get backward entry
            if nuc2.islower():
                backward = BTable[1][i+1]
            else:
                backward = BTable[0][i+1]
            # use equation to get counts
            oneTrans = forward + backward + transProb
            sumTrans = np.logaddexp(sumTrans,oneTrans)
        logCount = sumTrans-PofX
        count = np.exp(logCount) #get out of log space
    return count

def convertCountDToTransD(countD):
    '''Converts a matrix of counts of transitions to a matrix of 
    transition probabilities'''
    transD = {}
    for nuc1 in ['A','G','C','T','a','g','c','t']:
        # get sum for transitions from nuc1
        countSum = 0
        for nuc2 in ['A','G','C','T','a','g','c','t']:
            countSum += countD[nuc1+nuc2]
        for nuc3 in ['A','G','C','T','a','g','c','t']:
            # get prob for each transition pair
            transD[nuc1+nuc3] = math.log(countD[nuc1+nuc3]/countSum)
    return transD
