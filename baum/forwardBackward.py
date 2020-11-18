import math
import numpy as np

def forward(seq,transD):
    '''Takes a sequence seq and a matrix of transition probabilities transD. 
    Uses the forward algorithm to fill out a score table.'''
    tsScM = [[0 for i in range(len(seq))] for j in range(2)]

    for i in range(len(seq)):
        if i == 0:
            tsScM[0][i] = math.log(0.5)
            tsScM[1][i] = math.log(0.5)
        else:
            # Nonisland entry
            nonToNon = tsScM[0][i-1] + transD[seq[i-1].upper() + seq[i].upper()]
            isToNon = tsScM[1][i-1] + transD[seq[i-1].lower() + seq[i].upper()]
            toNon = np.logaddexp(nonToNon,isToNon)
            tsScM[0][i] = toNon

            # Island entry
            nonToIs = tsScM[0][i-1] + transD[seq[i-1].upper() + seq[i].lower()]
            isToIs = tsScM[1][i-1] + transD[seq[i-1].lower() + seq[i].lower()]
            toIs = np.logaddexp(nonToIs,isToIs)
            tsScM[1][i] = toIs

    return tsScM

def backward(seq,transD):
    '''Takes a sequence seq and a matrix of transition probabilities transD. 
    Uses the forward algorithm to fill out a score table.'''
    tsScM = [[0 for i in range(len(seq))] for j in range(2)]

    for i in range(len(seq)-1,-1,-1):
        if i == len(seq) - 1:
            tsScM[0][i] = 0
            tsScM[1][i] = 0
        else:
            # Nonisland entry
            nonToNon = tsScM[0][i+1] + transD[seq[i].upper() + seq[i+1].upper()]
            isToNon = tsScM[1][i+1] + transD[seq[i].upper() + seq[i+1].lower()]
            toNon = np.logaddexp(nonToNon,isToNon)
            tsScM[0][i] = toNon

            # Island entry
            nonToIs = tsScM[0][i+1] + transD[seq[i].lower() + seq[i+1].upper()]
            isToIs = tsScM[1][i+1] + transD[seq[i].lower() + seq[i+1].lower()]
            toIs = np.logaddexp(nonToIs,isToIs)
            tsScM[1][i] = toIs

    return tsScM
