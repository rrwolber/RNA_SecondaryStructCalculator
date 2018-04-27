


import sys, random
sys.setrecursionlimit(100000)
from rna import *
from hivSeqs import *

def rnaWin(RNA,wsize,step):
    """Runs a sliding window over RNA and uses mfold5 to calculate the optimal 
    number of matches in each window. WIll return a tuple comprising the score of 
    the optimal window, its start position, and its sequence."""
    bestSoFar=0
    bestTuple=(0,0,'')
    for i in range(0,len(RNA),step):
        memo={}
        score = mfold5(RNA[i:i+wsize],memo)
        if score>bestSoFar:
            bestSoFar=score
            bestTuple=(bestSoFar,i,RNA[i:i+wsize])
    return bestTuple

def randSeq(RNA):
    '''Takes an RNA string as input and returns a new string formed by
    randomly shuffling the symbols in the given string.'''
    L=list(RNA)
    random.shuffle(L)
    return "".join(L)

def randomWins(RNA, wsize, step, trials):
    """Takes an RNA input, shuffles it, calculates the best scoring window using
    rnaWin(wsize, step), repeats trials times."""
    outputList = []
    counter=0
    while counter < trials:
        newRNA=randSeq(RNA)
        aScore=rnaWin(newRNA,wsize,step)
        counter += 1
        outputList.append(aScore[0])
    return outputList

def pval(trueVal,randomWinsList):
    """Tells what proportion of the scores in the output from randomWins that 
    are greater than or equal to trueVal."""
    someList=[]
    for i in range(len(randomWinsList)):
        if randomWinsList[i] >= trueVal:
            someList.append(randomWinsList[i])
    output = (len(someList)*1.0)/(len(randomWinsList)*1.0)
    return output

def findSecStrucWrapper():
    '''Wrapper function to search HIV Pol gene for secondary structure.'''
    bestScoreHIV, bestPosHIV, bestSeqHIV = rnaWin(hivPol, 60, 30)
    randomWinsList = randomWins(hivPol,60,30,100)
    p = pval(bestScoreHIV,randomWinsList)
    print "Best window score:", bestScoreHIV
    print "Best window position:", bestPosHIV
    print "Best window seq:", bestSeqHIV
    print "Best window pvalue:", p