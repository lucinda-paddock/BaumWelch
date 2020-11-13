import math
import fasta
from context import *


def islandPedictionSummaryStats(predict,correct):
    """Print some statistics summarizing how good the prediction is
    relative to the correct solution. Our definitions of sensitivity
    and specificity assume the 'lower case' model represents the
    feature of interest."""
    trueLower=0
    trueUpper=0
    predictLowerAndTrueLower=0
    predictUpperAndTrueUpper=0
    for i in range(len(predict)):
        if correct[i].islower(): trueLower+=1
        else: trueUpper+=1

        if predict[i].islower() and correct[i].islower(): predictLowerAndTrueLower+=1
        if predict[i].isupper() and correct[i].isupper(): predictUpperAndTrueUpper+=1

    sensitivity = 1.0 * predictLowerAndTrueLower / trueLower
    specificity = 1.0 * predictUpperAndTrueUpper / trueUpper
    
    print("Sensitivity (TPR):", format(sensitivity,".3f"))
    print("Specificity (TNR):", format(specificity,".3f"))
    print("Youden's J       :", format(sensitivity + specificity - 1,".3f"))
    

def runExample():
    """Wrapper to run the provided cpg example."""
    
    states = ["A","C","G","T","a","c","g","t"]
    
    # calculate transition matrix from using A solution for training
    Asolution=fasta.load("A-cpgSolution.fa")[0][1]
    cpgTransD,cpgCountD=context([Asolution],1,states,1)
    cpgTransD = logDictValues(cpgTransD)
    
    # now use viterbi to get sol for B
    B=fasta.load("B-region.fa")[0][1]
    scM,btM=viterbi(B,cpgTransD)
    solB=bt(B,scM,btM)

    f=open("Bpredict.fa","w")
    print("> prediction",file=f)
    print(solB,file=f)
    f.close()

    # test out our results on B vs. the provided solution
    Bsolution=fasta.load("B-cpgSolution.fa")[0][1]
    islandPedictionSummaryStats(solB,Bsolution)

def viterbi(seq,transD):
    '''takes a sequence seq and a matrix of transition probabilities transD (in the form of a 
    dictionary) as input. It uses the viterbi algorithm to fill out a score table and a backtracking 
    table with these, and then returns these tables.'''
    scoreTab = [[],[]]
    scoreTab[0].append(math.log(0.5)) #placeholder
    scoreTab[1].append(math.log(0.5)) #placeholder
    btTab = [[],[]]
    btTab[0].append("NA")
    btTab[1].append("NA")

    for i in range(len(seq)-1):

        #some values to use in our calculations
        oldMinus = scoreTab[0][i]
        oldPlus = scoreTab[1][i]
        oldNuc = seq[i]
        newNuc = seq[i+1]

        #assuming -
        staysMinus = oldMinus + transD[oldNuc+newNuc]
        becomesMinus = oldPlus + transD[oldNuc.lower()+newNuc]
        if staysMinus > becomesMinus:
            scoreTab[0].append(staysMinus)
            btTab[0].append("-")
        else:
            scoreTab[0].append(becomesMinus)
            btTab[0].append("+")

        #assuming +
        staysPlus = oldPlus + transD[oldNuc.lower()+newNuc.lower()]
        becomesPlus = oldMinus + transD[oldNuc+newNuc.lower()]
        if staysPlus > becomesPlus:
            scoreTab[1].append(staysPlus)
            btTab[1].append("+")
        else:
            scoreTab[1].append(becomesPlus)
            btTab[1].append("-")
    
    return scoreTab, btTab

def bt(seq,scM,btM):
    '''Takes the score and backtrack tables (scM and btM) as input along with a sequence. 
    It uses these to calculate the most probable path (series of states) through seq. It 
    returns this as a string, with CpG islands indicated by lower case letters'''
    
    path = ""

    #start backtracking w the highest of the final scores
    if max(scM[0][-1],scM[1][-1]) == scM[0][-1]:
        lastVal = "-"
    else: 
        lastVal = "+"

    for i in range(len(seq)-1,-1,-1):
        if lastVal == "-":
            path = seq[i] + path
            if btM[0][i] == "+":
                lastVal = "+"
        else:
            path = seq[i].lower() + path
            if btM[1][i] == "-":
                lastVal = "-"
    return path




'''

ts="AAGCAGTA"
states = ["A","C","G","T","a","c","g","t"]
Asolution=fasta.load("A-cpgSolution.fa")[0][1]
cpgTransD,cpgCountD=context([Asolution],1,states,1)
cpgTransD = logDictValues(cpgTransD)
tsScM,tsbtM=viterbi(ts.upper(),cpgTransD)

tsbtM = [['NA', '-', '-', '-', '-', '-', '-', '+'], ['NA', '+', '+', '+', '-', '+', '+', '+']]
bt(ts,tsScM,tsbtM)
'''
