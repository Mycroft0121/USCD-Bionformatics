def frequentWords(text, k):
    FrequentPatterns = []
    count = []
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        count.append(patterncount(text, pattern))
    maxcount = max(count)
    for i in range(len(text)-k+1):
        if count[i] == maxcount:
            FrequentPatterns.append(text[i:i+k])    
    return set(FrequentPatterns)
    
def patterncount(text, pattern):
    count = 0
    i = 0
    end = len(text) - len(pattern)
    for i in range(end+1):
        if text[i:i+len(pattern)] == pattern:
            count +=1
        i +=1
    return count


def PatternToNumber(pattern):
    DNAdict = {'A':0, 'C':1,'G':2,'T':3}
    pattern = pattern[::-1]
    digit = len(pattern)
    num = 0
    for i in range(digit):
        x = pattern[i]
        num += 4**i*DNAdict[x]
    return (num, digit)
   
def NumberToPattern(num, digit):
    DNAdict = {0:'A', 1:'C', 2:'G', 3:'T'}
    if digit == 1:
        return DNAdict[num]
    pref, r = num/4, DNAdict[num%4]
    prefPattern = NumberToPattern(pref, digit-1)
    return prefPattern+r
 
def ComputingFrequencies(text, k):
    n = 4**k
    frequencyArray = [0]*n
    location = {m:[] for m in range(n)}
    for i in range(len(text)-k+1):
        pattern = text[i:i+k]
        j = PatternToNumber(pattern)[0]
        frequencyArray[j] +=1
        location[j].append(i)
    return frequencyArray, location
    
def reversecomplement(pattern):
    complement = {"A":"T", "C":"G","G":"C","T":"A"}
    revcom = ""
    for char in pattern:
        revcom +=complement[char]
    return revcom[::-1]
        
def findpattern(pattern, text):
    location =""
    length = len(pattern)
    for i in range(len(text)-length):
        if pattern == text[i:i+length]:
            location += str(i) +" "
    return location
    
def ClumpFinding(genome, k, L, t):
    freq, loc = ComputingFrequencies(genome, k)
    rough = []
    result =[]
    
    for i in range(len(freq)):
        if freq[i]>=t:
            rough.append(i)
            
    for char in rough:
        loci =loc[char]
        for j in range(len(loci)-t+1):
            if loci[j+t-1]-loci[j]<=L-k:
                result.append(NumberToPattern(char,k))
                break
    return result, len(result)
                                                                   
def minskew(text):
    skew = [0]
    skewness={"A":0,"C":-1,"G":1,"T":0,"\n":0}
    for i in range(len(text)):
        skew.append(skew[i]+skewness[text[i]])
    return skew.index(min(skew))
    
def hammingdist(text1, text2):
    score = 0
    for i in range(len(text1)):
        if text1[i]!=text2[i]:
            score+=1
    return score
    
def approPatternLoci(pattern, text, score):
    loci = []
    for i in range(len(text)-len(pattern)+1):
        if hammingdist(pattern,text[i:i+len(pattern)])<=score:
            loci.append(i)
    return loci
    
def countd(text, pattern, score):
    return len(approPatternLoci(pattern, text, score))
    
def findingFrequencyBySortying(text, k):
    frequentPatterns = []
    tlen = len(text)
    index= []
    count = []
    for i in xrange(tlen-k+1):
        pattern = text[i:i+k]
        index.append(PatternToNumber(pattern)[0])
        count.append(1)
    sortedIndex = sorted(index)
    for i in xrange(1,tlen-k+1):
        if sortedIndex[i] == sortedIndex[i-1]:
            count[i] = count[i-1]+1
    maxCount = max[count]
    for i in xrange(tlen-k+1):
        if count[i]==maxCount:
            pattern = NumberToPattern(sortedIndex(i),k)
            frequentPatterns.append(pattern)
    return set(frequentPatterns)
        
def immediateNeighbors(pattern,d):
    neighborhood = set([])
    Code = ["A","C","G","T"]
    for i in xrange(1,len(pattern)):
        symbol = pattern[i]
        for x in Code and x!= symbol:
            neibor = list(pattern[:])
            neibor[i]=x
            neighborhood.update([str(neibor)])
    return neighborhood

def neighbors(pattern, d):
    if d==0:
        return {pattern}
    if len(pattern)==1:
        return {"A","C","T","G"}
    neighborhood = set([])
    suffixNeighbors = neighbors(pattern[1:],d)
    for x in suffixNeighbors:
        if hammingdist(pattern[1:],x)<d:
            for y in ["A","C","T","G"]:
                neighborhood.update([y+x])
        else:
            neighborhood.update([pattern[0]+x])
    return list(neighborhood)

def FrequentWordsMismatches(text,k,d):
    freqpattern = []
    close, frequencyArray = [0]*(4**k), [0]*(4**k)
    for i in xrange(len(text)-k+1):
        neighborhood = neighbors(text[i:i+k],d)
        for Pattern in neighborhood:
            index = PatternToNumber(Pattern)[0]
            close[index] = 1
    for i in xrange(4**k):
        if close[i] == 1:
            Pattern = NumberToPattern(i,k)
            frequencyArray[i]= countd(text, Pattern, d)
    maxCount = max(frequencyArray)
    for i in xrange(4**k):
        if frequencyArray[i] == maxCount:
            freqpattern.append(NumberToPattern(i, k))
    return freqpattern

        
def FrequentWordsMismatchesSorting(text,k,d):
    freqpattern = []
    Neighborhood = []
    for i in xrange(len(text)-k+1):
        Neighborhood.append(neighbors(text[i:i+k],d))
    Neighborhoods = []
    for y in Neighborhood:
        for x in y:
            Neighborhoods.append(x)
    
    index, count = [0]*len(Neighborhoods),[0]*len(Neighborhoods)
    for i in xrange(len(Neighborhoods)):
        pattern = str(Neighborhoods[i])      
        index[i] = PatternToNumber(pattern)[0]
        count[i] = 1
    sortedindex = sorted(index)
    for i in xrange(len(Neighborhoods)-1):
        if sortedindex[i] == sortedindex[i+1]:
            count[i+1]=count[i]+1
    maxcount = max(count)
    for i in xrange(len(Neighborhoods)):
        if count[i] == maxcount:
            freqpattern.append(NumberToPattern(sortedindex[i],k))
    return freqpattern
    
def FrequentWordsMismatchesReverseComplementsSorting(text,k,d): 
    freqpattern = []
    Neighborhood = []
    textrc = reversecomplement(text)
    for i in xrange(len(text)-k+1):
        Neighborhood.append(neighbors(text[i:i+k],d))
        Neighborhood.append(neighbors(textrc[i:i+k],d))
    Neighborhoods = []
    for y in Neighborhood:
        for x in y:
            Neighborhoods.append(x)
    
    index, count = [0]*len(Neighborhoods),[0]*len(Neighborhoods)
    for i in xrange(len(Neighborhoods)):
        pattern = str(Neighborhoods[i])      
        index[i] = PatternToNumber(pattern)[0]
        count[i] = 1
    sortedindex = sorted(index)
    for i in xrange(len(Neighborhoods)-1):
        if sortedindex[i] == sortedindex[i+1]:
            count[i+1]=count[i]+1
    maxcount = max(count)
    for i in xrange(len(Neighborhoods)):
        if count[i] == maxcount:
            freqpattern.append(NumberToPattern(sortedindex[i],k))
    return freqpattern            
       
text = "CACTTACACTATTCTTCTTCTATTATTACACCACCACCACCACTTCCACTATTATTACACTTATATACACCACTATACACTTCTTATTCCACCACCACTTCTTCTATTACACTTACACCACTTCTTACACTACACCACTATATATTCCACTATTACACTTATTACACTATTCTTCCACTTCCACTTCCACCACTACACTATTCCACCACTATTCTTCTTACACTTCTTCTTCTA"
k = 10
d = 3
print FrequentWordsMismatchesReverseComplementsSorting(text,k,d)
    