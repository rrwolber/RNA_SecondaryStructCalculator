

def isComplement(base1,base2):
    """Returns boolean indicating if 2 RNA bases are complementary."""
    if base1=="A" and base2=="U":
        return True
    elif base1=="U" and base2=="A":
        return True
    elif base1=="C" and base2=="G":
        return True
    elif base1=="G" and base2=="C":
        return True
    if base1=="G" and base2=="U":
        return True
    elif base1=="U" and base2=="G":
        return True
    else:
        return False

def getStruct(RNA,memo):
    """Returns a tuple in which the first element is the maximum number of matches
    and the second element is (x,y) where x and y are the indices of two nucleotides
    that are matched in an optimal folding."""
    if RNA in memo:
        return memo[RNA]
    elif len(RNA)<6:
        return (0,[])
    else:
        bestSoFar=getStruct(RNA[1:],memo) # lose it case
        bestSoFar=(bestSoFar[0], adjust(bestSoFar[1],1))
        for i in range(5,len(RNA)):       # use it cases
            if isComplement(RNA[0],RNA[i]):
                score1=getStruct(RNA[1:i],memo)
                score1=(score1[0], adjust(score1[1],1))
                score2=getStruct(RNA[(i+1):],memo)
                score2=(score2[0], adjust(score2[1], i+1))
                aThing=(score1[0]+score2[0]+1, [(0,i)] + score1[1]+score2[1])
                if aThing[0]>=bestSoFar[0]:
                    bestSoFar=aThing
        memo[RNA]=bestSoFar
        return bestSoFar

def adjust(pairs, k):
    """Takes as input a list of pairs of numbers and an integer k and returns a 
    new tuple that increases every number in the pairs by k."""
    output = []  
    for i in range(len(pairs)):
        newPairs=(pairs[i][0]+k,pairs[i][1]+k)
        output.append(newPairs)
    return output