import math

def makeFrags(fragLen, alphabet):
    """Return all possible DNA fragments of len fragLen."""
    if fragLen==0:
        return [""]
    else:
        outL=[]
        for base in alphabet:
            for frag in makeFrags(fragLen-1,alphabet):
                outL.append(base+frag)
        return outL

def context(seqL,k,alphabet,pseudocount):
    """Take list of seqs as input, make probability matrix for k order
    markov chain (as dictionary) and return it. E.g. if k is 2 and
    alphabet is the nucleotides, we'll have keys like "AC" in our
    dictionary. This will tell us the probability given an A that the
    next base is a C. We can add pseudocount to all cases to make sure
    there are no transitions with zero probability."""
    fragLen=k+1
    fragL=makeFrags(fragLen,alphabet)
    countD={}
    for frag in fragL: countD[frag]=pseudocount
    originFragD={}
    for frag in makeFrags(fragLen-1,alphabet):
        originFragD[frag]=len(alphabet)*pseudocount
    for seq in seqL:
        for i in range(len(seq)-1):
            if seq[i:i+fragLen] in fragL:
                countD[seq[i:i+fragLen]]+=1
                originFragD[seq[i:(i+fragLen-1)]]+=1
    probD={}
    for key in countD:
        probD[key] = countD[key]/float(originFragD[key[:fragLen-1]])
    return probD,countD

def logDictValues(D):
    '''Return a new dict whose values are the log of the input dict's values.'''
    newD={}
    for key in D:
        newD[key] = math.log(D[key])
    return newD

def exponentiateDictValues(D):
    '''Return a new dict whose values have been exponentiated.'''
    newD={}
    for key in D:
        newD[key] = math.exp(D[key])
    return newD

def printProbs(probD,nucs,formatString):
    """Print probD, organized according to the prefixes of its
    keys. Print to precision given in formatString (e.g. '.3f' etc.)."""
    S=set()
    for key in probD.keys():
        S.add(key[:-1])
    prefixes=sorted(list(S))
    for prefix in prefixes:
        for base in nucs:
            print(prefix+base+":"+format(probD[prefix+base],formatString),end=' ')
        print()
