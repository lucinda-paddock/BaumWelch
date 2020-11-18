import math,random,sys
import fasta
import forwardBackward as fb
import baumWelch as bw
import numpy as np
import cpgViterbi as vit

DELTA = 0.01

## Main

if __name__ == "__main__":

    if len(sys.argv) !=3:
        sys.stderr.write("""
        Usage: python testSeq.fa outputPredictSeq.fa

""")
        sys.exit(-1)

    testSeqFN = sys.argv[1]
    outFN = sys.argv[2]

    ## Load sequences
    testSeq = fasta.load(testSeqFN)[0][1]
    solSeq = fasta.load("solutionSeq.fa")[0][1]
    # start loop
    transD = bw.initTransD()
    FTable = fb.forward(testSeq,transD)
    crit = np.logaddexp(FTable[0][-1],FTable[1][-1])
    while True:
        countD = bw.makeCountD(testSeq,transD,FTable)
        newTransD = bw.convertCountDToTransD(countD)
        newFTable = fb.forward(testSeq, newTransD)
        newCrit = np.logaddexp(newFTable[0][-1],newFTable[1][-1])
        # check criteria for substantial improvement
        if newCrit - crit < DELTA:
            break
        else:
            transD = newTransD
            FTable = newFTable
            crit = newCrit
    # predict island states with viterbi
    scoreTab, btTab = vit.viterbi(testSeq, transD)
    predictSeq = vit.bt(testSeq, scoreTab, btTab)
    vit.islandPredictionSummaryStats(predictSeq,solSeq)
        
    # write prediction (in form of sequence predictSeq, which has
    # the states as upper and lower case) to file.
    f=open(outFN,'w')
    f.write('>predictSeq\n'+predictSeq+'\n')
    f.close()
    