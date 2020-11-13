import math,numpy,random,sys
import fasta


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
    testSeq=fasta.load(testSeqFN)[0][1]













    
    # write prediction (in form of sequence predictSeq, which has
    # the states as upper and lower case) to file. You will want to
    # comment out these lines until your code actually produces
    # predictSeq.
    f=open(outFN,'w')
    f.write('>predictSeq\n'+predictSeq+'\n')
    f.close()
    
