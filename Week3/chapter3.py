# -*- coding: utf-8 -*-

#chapter 3

count = 0
import random
from scipy import stats
#Helper functions
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

def hammingdist(text1, text2):
    score = 0
    for i in range(len(text1)):
        if text1[i]!=text2[i]:
            score+=1
    return score

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

def NumberToPattern(num, digit):
    DNAdict = {0:'A', 1:'C', 2:'G', 3:'T'}
    if digit == 1:
        return DNAdict[num]
    pref, r = num/4, DNAdict[num%4]
    prefPattern = NumberToPattern(pref, digit-1)
    return prefPattern+r

def PatternToNumber(pattern):
    DNAdict = {'A':0, 'C':1,'G':2,'T':3}
    pattern = pattern[::-1]
    digit = len(pattern)
    num = 0
    for i in range(digit):
        x = pattern[i]
        num += 4**i*DNAdict[x]
    return (num, digit)

#Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
#Input: A collection of strings Dna in list, and integers k and d.
#Output: All (k, d)-motifs in Dna.
def MotifEnumeration(DNA, k,d):
    patterns = [[] for x in xrange(len(DNA))]
    numKmer = len(DNA[0])-k+1
    for i in xrange(len(DNA)):
        dna = DNA[i]
        for j in xrange(numKmer):
            neighbor = neighbors(dna[j:j+k], d)
            patterns[i] += neighbor
    result = set(patterns[0])
    for p in patterns[1:]:
        result.intersection_update(p)
    result = list(result)  
    result = " ".join(str(x)for x in result)
    return result       
    

#Implement DistanceBetweenPatternAndStrings.
#Input: A string Pattern followed by a collection of strings Dna.
#Output: d(Pattern, Dna)    
def DistanceBetweenPatternAndStrings(Pattern, DNA):
    plen = len(Pattern)
    dist = 0
    numPattern = len(DNA[0])-len(Pattern)+1
    for dna in DNA:            
        mindist = float("inf")
        for i in xrange(numPattern):
            score = hammingdist(Pattern, dna[i:i+plen])
            if score <mindist:
                mindist = score
        dist +=mindist
    return dist
    
#Median String Problem: Find a median string.
#Input: A collection of strings Dna and an integer k.
#Output: A k-mer Pattern that minimizes d(Pattern, Dna) among all k-mers Pattern.                
def medianString(DNA,k):
    distancemin = float("inf")
    for i in xrange(4**k):
        pattern = NumberToPattern(i, k)
        dist = DistanceBetweenPatternAndStrings(pattern, DNA)
        if dist < distancemin:
            distancemin = dist
            median = pattern
    return median
    
    
#Find a Profile-most probable k-mer in a string.
# Input: A string Text, an integer k, and a 4 × k matrix Profile.
# Output: A Profile-most probable k-mer in Text.    
def mostProbableKmer(Text, k, probMatrix):
    likelykmer = ""
    numKmer = len(Text)-k+1
    DNAdict = {'A':0, 'C':1,'G':2,'T':3}
    prob = -1    
    for i in xrange(numKmer):
        kmer = Text[i:i+k]
        probkmer = 1
        for j in xrange(k):
            probkmer *= float(probMatrix[DNAdict[kmer[j]]][j])
        if probkmer > prob:
            prob = probkmer
            likelykmer = kmer
    return likelykmer
'''
f=open('dataset_159_3.txt','r')
Text=f.readline()
Text=Text[:len(Text)-1]
k=int(f.readline())
probMatrix =[[] for i in range(4)]
for j in range(4): 
    p = f.readline().split()
    for num in p:
        probMatrix[j].append(float(num))
f.close()'''

#return the score
def score(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]
    max_count = sum([max([c.count(nucleotide) for nucleotide in 'ACGT']) for c in columns])
    return len(motifs[0])*len(motifs) - max_count

#Profile input: a matrix of motifs, output: a matrix of probability
def profile(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]
    p = [[ float(col.count(nuc)) / float(len(col)) for nuc in 'ACGT'] for col in columns]
    return [list(x) for x in zip(*p)]
#Profile input: a matrix of motifs, output: a matrix of probability using pseudo probability
def profileWithSuccession(motifs):
    columns = [''.join(seq) for seq in zip(*motifs)]

    p = [[ (float(col.count(nuc))+1) / (float(len(col))*2) for nuc in 'ACGT'] for col in columns]
    return [list(x) for x in zip(*p)]


#GREEDYMOTIFSEARCH
#input: a list of DNA strings, int k of len of Kmer, int t
#Output: A collection of strings BestMotifs resulting from applying GREEDYMOTIFSEARCH(Dna,k,t).
#If at any step you find more than one Profile-most probable k-mer in a given string, use the
#one occurring first
def GreedyMotifSearch(DNA, k, t):
    BestMotifs = [DNA[i][:k] for i in xrange(t)]
    Bestscore = score(BestMotifs)
    numkmer= len(DNA[0])-k+1
    for i in xrange(numkmer):
        Motif = [DNA[0][i:i+k]]
        
        for j in xrange(1,t):
            Profile= profile(Motif)
            Motif.append(mostProbableKmer(DNA[j], k,Profile ))
        newscore = score(Motif)
        if newscore < Bestscore:
            BestMotifs = Motif
            Bestscore = newscore
    BestMotifs = " ".join(str(x)for x in BestMotifs)
    return BestMotifs
        
def GreedyMotifSearchwithSuccession(DNA, k, t):
    BestMotifs = [DNA[i][:k] for i in xrange(t)]
    Bestscore = score(BestMotifs)
    numkmer= len(DNA[0])-k+1
    for i in xrange(numkmer):
        Motif = [DNA[0][i:i+k]]
        
        for j in xrange(1,t):
            Profile= profileWithSuccession(Motif)
            Motif.append(mostProbableKmer(DNA[j], k,Profile ))
        newscore = score(Motif)
        if newscore < Bestscore:
            BestMotifs = Motif
            Bestscore = newscore
    BestMotifs = " ".join(str(x)for x in BestMotifs)
    return BestMotifs        

        
                        
'''f=open('dataset_160_9.txt','r')
text = f.readline()
k, t = int(text[0:2]), int(text[3:5])
DNA = []
for line in f:
    DNA.append(line[:-1])
f.close()
GreedyMotifSearchwithSuccession(DNA, k, t)'''

#input: a probability profile matrix, DNA list of strings
#output: motifs in list of string
def motifs(Profile, DNA,k):
    Motifs = []
    for dna in DNA:
        Motifs.append(mostProbableKmer(dna, k, Profile))
    return Motifs
    
#randomizedMotifSearch
#input: a list of DNA strings, int k of len of Kmer, int t
#output: best motifs as string after iterations

 

def randomizedMotifSearch(DNA, k, t):
    numkmer = len(DNA[0])-k+1
    Motifs = []
    for dna in DNA:
        start = random.randint(0,numkmer-1)
        Motifs.append(dna[start:start+k])
    BestMotifs = Motifs    
    while True:
        #Profile = ProfileOfMotifs(Motifs,1)
        Profile = profileWithSuccession(Motifs)
        #Motifs = bestMotifsForProfile(DNA, Profile)
        Motifs = motifs(Profile, DNA, k)
        Score = score(Motifs)
        print Score
        if Score < score(BestMotifs):
            BestMotifs = Motifs
        else:
            return BestMotifs

#randomizedMotifSearch
#input: a list of DNA strings, int k of len of Kmer, int t, int N of number of runs
#output: best motifs as string after N runs        
def randSearchRepeat(DNA, k, t,n):   
    result = []
    numkmer = len(DNA[0])-k+1
    bestMotifs = []
    for dna in DNA:
        start = random.randint(0,numkmer-1)
        bestMotifs.append(dna[start:start+k])
    minscore = score(bestMotifs)
    print minscore
    for i in xrange(n+1):
        newmotif = randomizedMotifSearch(DNA, k, t)
        newscore = score(newmotif)
        if newscore< minscore:
            minscore = newscore
            result = newmotif
    result = " ".join(str(x)for x in result)

    return minscore, result

'''n =100
f=open('dataset_161_5.txt','r')
text = f.readline()
k, t = int(text[0:2]), int(text[3:5])
DNA = []
for line in f:
    DNA.append(line[:-1])
f.close()
print randSearchRepeat(DNA, k, t,n)'''

#gibbsRandom input: int t, output: a int generated from the generated CDF
def gibbsRandom(t):
    base = [random.randint(1,t) for x in range(t)] 
    csum = sum(base)
    pfile = [float(x)/csum for x in base]
    xk = range(t)
    cdf = stats.rv_discrete(values=(xk, pfile))
    return cdf.rvs()


#random select 1 k-mer from 1 DNA
def randKmer(dna, k):
    nkmer = len(dna)-k+1
    index = random.randint(0,nkmer-1)
    return dna[index:index+k]
#random select t k-mer from t DNA, creating random Motifs
def randMotifs(DNA, k):
    return [randKmer(dna, k) for dna in DNA]

#profRandKmer input: a dna string, int k, and a prob profile matrix
def profRandKmer(dna, k,t, Profile): 
     DNAdict = {'A':0, 'C':1,'G':2,'T':3}  
     numkmer= len(dna)-k+1
     kmerProb = []
     kmers = [dna[i:i+k] for i in xrange(numkmer)]
     i = 0
     for kmer in kmers:
         prob =1.0
         for idx, char in enumerate(kmer):
             prob *= Profile[DNAdict[char]][idx]
         kmerProb.append((prob, kmers[i]))
         i+=1
     total = sum(x[0] for x in kmerProb)
     randval = random.uniform(0, total)
     probability_acc = 0
     chosen_kmer = ''
     while probability_acc < randval:
         prob, chosen_kmer = kmerProb.pop()
         probability_acc += prob

     return chosen_kmer
     '''for i in xrange(numkmer):
         kmer = dna[i:i+k]
         prob = 1.0
         for j in xrange(k):
             kmernum = DNAdict[kmer[j]]
             #print kmernum,len(Profile[kmernum]), j
             prob *= float(Profile[kmernum][j])
         kmerProb.append(prob)
     csum = sum(kmerProb)
     xk = range(t)
     CDF = [float(prob)/csum for prob in kmerProb]
     gibbsRandom =  stats.rv_discrete(values=(xk,CDF))
     randnum = gibbsRandom.rvs()
     return dna[randnum:randnum+k]'''
     
 
                                       
#Input: Integers k, t, and N, followed by a collection of strings Dna.
#Output: The strings BestMotifs resulting from running GIBBSSAMPLER(Dna, k, t, N) with
#20 random starts. Remember to use pseudocounts!
def gibbsSampler(DNA, k, t, N):
    start = [randMotifs(DNA, k) for x in xrange(20)]
    startscore = [score(motif) for motif in start]
    bestMotifs = sorted(zip(startscore, start))[0][1]
    bestscore = score(bestMotifs)
    Motifs = [x for x in bestMotifs]
    for j in xrange(N):
        i = random.randint(0,t-1)
        Motifs.pop(i)

        Profile = profileWithSuccession(Motifs)
        a = profRandKmer(DNA[i],k,t,Profile)
        Motifs.insert(i,a)
        if score(Motifs) < bestscore:
            bestscore = score(Motifs)
            bestMotifs = [x for x in Motifs]
                 
    return bestMotifs, score(bestMotifs)

def repeatGibbs(DNA, k, t, N, repeats):
    best = float('inf')
    bestMotif = []
    for i in xrange(repeats):
        m, s = gibbsSampler(DNA, k, t, N)
        if s < best:
            best = s
            bestMotif = m
    return best, bestMotif
    
f=open('dataset_163_4.txt','r')
text = f.readline()
k, t= int(text[0:2]), int(text[3:5])
N=2000
DNA = []
for line in f:
    DNA.append(line[:-1])
f.close()

'''DNA = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
     'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
k = 8
t=5
N = 100'''
repeats = 10
#print gibbsSampler(DNA, k, t, N)
a = repeatGibbs(DNA, k, t, N, repeats)
print a[0], " ".join(str(x)for x in a[1])
