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

'''
transD = {'AA': 0.4396583462189824, 'AC': 0.22516751343789118, 'AG': 0.16478904351667772, 'AT': 0.14910536779324055, 'Aa': 0.010529416096016493, 'Ac': 0.01045578381562477, 'Ag': 0.0002208968411751712, 'At': 7.363228039172373e-05, 'CA': 0.2725360998177485, 'CC': 0.34550679938314877, 'CG': 0.08495724099256975, 'CT': 0.2576054955839058, 'Ca': 0.010234123089864012, 'Cc': 0.010304219823356232, 'Cg': 0.009322865554465163, 'Ct': 0.009533155754941819, 'GA': 0.20253164556962025, 'GC': 0.25724237532919886, 'GG': 0.34661456120975276, 'GT': 0.15164387052926684, 'Ga': 0.009939682269985558, 'Gc': 0.010619318664514484, 'Gg': 0.010789227763146716, 'Gt': 0.010619318664514484, 'TA': 0.10183440753594447, 'TC': 0.2887456618740704, 'TG': 0.36331184928111054, 'TT': 0.20555280118988598, 'Ta': 0.00932077342588002, 'Tc': 0.011303916707982151, 'Tg': 0.00932077342588002, 'Tt': 0.010609816559246405, 'aA': 0.008203725858827551, 'aC': 0.01076739018971116, 'aG': 0.010596479234318919, 'aT': 0.011792855922064605, 'aa': 0.13570329858143906, 'ac': 0.2948213980516151, 'ag': 0.4385575115364895, 'at': 0.0895573406255341, 'cA': 0.010843967358562901, 'cC': 0.008488964346349746, 'cG': 0.010624897310915165, 'cT': 0.010953502382386768, 'ca': 0.17257243003450354, 'cc': 0.3759242017635139, 'cg': 0.2547236979024043, 'ct': 0.1558683389013637, 'gA': 4.99325910021471e-05, 'gC': 0.0031956858241374143, 'gG': 0.013132271433564687, 'gT': 0.00948719229040795, 'ga': 0.05202975982423728, 'gc': 0.3416887202276926, 'gg': 0.4641733659559594, 'gt': 0.11624307185299845, 'tA': 0.01012210796915167, 'tC': 0.009158097686375322, 'tG': 0.009318766066838046, 'tT': 0.010925449871465296, 'ta': 0.0586439588688946, 'tc': 0.3693766066838046, 'tg': 0.5072300771208226, 'tt': 0.025224935732647814}
newTransD = {}
for key,value in transD.items():
    newTransD[key] = math.log(value)
full = 1.8497407732746777 * 10**(-28)
seq = "GATGGCCCCGTCACTGGCCTCTCTACCTTGGCCCGCACTGGCGCCGCGCA"
'''