from forwardBackward import *

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


def makeCountD(seq,transD):
    '''Given a matrix of transition frequencies and a sequence, uses
    the forward and backward algorithms to construct a matrix of expected
    counts according to the Baum-Welch algorithm'''
    countD = {}
    FTable = forward(seq,transD)
    BTable = backward(seq, transD)
    transIndexD = makeTransIndexD(seq)
    PofX = np.logaddexp(forward[0][-1],forward[1][-1])
    for key,value in transD.items():
        countD[key] = calcCount(key, value, FTable, BTable, transIndexD, PofX)
    return countD

def calcCount(nucPair, transProb, FTable, BTable, transIndexD, PofX):
    '''For a given transtion (like Ag for example) calculate the expected
    count'''
    if not nucPair in transIndexD:
        count = 1 #add 1 as a pseudocount
    else:
        sumTrans = -float('inf') #log of 0
        transIndexTup = transIndexD[nucPair]
        for i in transIndexTup:
            nuc1 = nucPair[0]
            nuc2 = nucPair[1]
            if nuc1.islower():
                forward = FTable[1][i]
            else: 
                forward = FTable[0][i]
            if nuc2.islower():
                backward = BTable[1][i+1]
            else:
                backward = BTable[0][i+1]
            oneTrans = forward + backward + log(transProb)
            sumTrans = np.logaddexp(sumTrans,oneTrans)
        logCount = sumTrans-PofX
        count = np.exp(logCount) #get out of log space
    return count

def convertCountDToTransD(countD):
    '''Converts a matrix of counts of transitions to a matrix of 
    transition probabilities'''
    transD = {}
    for nuc1 in ['A','G','C','T','a','g','c','t']:
        countSum = 0
        for nuc2 in ['A','G','C','T','a','g','c','t']:
            countSum += countD[nuc1+nuc2]
        for nuc3 in ['A','G','C','T','a','g','c','t']:
            transD[nuc1+nuc3] = countD[nuc1+nuc3]/countSum
    return transD






