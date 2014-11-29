# You do not need to edit this for the lab.

# common.py
# This module defines some common functions for dealing with motifs and
# profiles, as in Chapter 3 of our textbook.

# Convert back and forth between acids (A,C,G,T) and indices (01,2,3).
def acidOfIndex(i):
    "Return acid symbol (A,C,G,T) of the given index (0,1,2,3)."
    return "ACGT"[i]

def indexOfAcid(a):
    "Return index (0,1,2,3) of the given acid symbol (A,C,G,T)."
    return "ACGT".index(a)

# Given a Motifs array, compute its Profile.
def ProfileOfMotifs(motifs, initCount):
    """
    Given Motifs, a list of DNA k-mers, return its Profile array.
    Profile will be a list of k lists, each a list of 4 floats.
    Each length 4 list is a probability distribuition over the
    4 nucleic acids (A, C, G, T).
    The initCount parameter is passed on to CountsOfMotifs.
    """
    counts = CountsOfMotifs(motifs, initCount)
    profile = []
    for count in counts:
        total = float(sum(count))
        probs = [c/total for c in count]
        profile.append(probs)
    return profile

# Given a Motifs array, compute its Counts.
def CountsOfMotifs(motifs, initCount):
    """
    Given Motifs, a list of DNA k-mers, return its Counts array.
    Counts will be a list of k lists, each a list of 4 ints.
    Each int counts the number of occurrences of the corresponding
    nucleic acids (A, C, G, T).
    We start each counter with the initial value initCount, so set
    this to 0 if you want true counts, or 1 for pseudocounts.
    """
    k = len(motifs[0])
    for kmer in motifs:
        assert k == len(kmer)
    counts = []
    for i in range(k):
        count_i = [initCount] * 4
        for kmer in motifs:
            count_i[indexOfAcid(kmer[i])] += 1
        counts.append(count_i)
    return counts

def Score(motifs):
    """
    Return the score of Motifs, a list of k-mers.
    That is, the number of symbols in the Motifs array
    that would not be in the consensus string.
    """
    score = 0
    for count in CountsOfMotifs(motifs, 0): # true counts here
        score += sum(count) - max(count)
    return score

def probOfKmerGivenProfile(kmer, profile):
    "Return probability that kmer would appear, given a profile."
    assert len(kmer) == len(profile)
    prob = 1.0
    for i in range(len(kmer)):
        prob *= profile[i][indexOfAcid(kmer[i])]
    return prob

def bestKmerForProfile(s, profile):
    """
    Return a most probable k-mer in DNA string s, given a profile.
    If there is a tie, we return the leftmost such k-mer.
    """
    k = len(profile)
    assert len(s) >= k
    maxprob = -1
    maxkmer = ''
    for start in range(len(s)-k+1):
        kmer = s[start:start+k]
        prob = probOfKmerGivenProfile(kmer, profile)
        if prob > maxprob:
            maxprob = prob
            maxkmer = kmer
    return maxkmer

def bestMotifsForProfile(dnas, profile):
    """
    Given dnas and a profile, return the most probable Motifs.
    That is, pick a most probable k-mer from each DNA string.
    """
    return [bestKmerForProfile(dna, profile) for dna in dnas]

def readDnasAndK(file):
    """
    Read problem input from the open file, and return pair (Dnas, k).
    The input should be in the format of rosalind.info problems 3e and 3f.
    """
    tokens = file.readline().split()    # split up the first line
    k = int(tokens[0])                  # length of motif
    t = int(tokens[1])                  # number of DNA strings to read
    Dnas = []
    for i in range(t):
        Dnas.append(file.readline().strip())
    # Note we return a pair of values.
    # We do not need to return t, since t == len(Dnas).
    return (Dnas, k)

# Random functions:
import random

def randomKmer(dna, k):
    "Given dna and k, return a randomly chosen k-mer from the dna."
    start = random.randint(0, len(dna)-k)
    return dna[start:start+k]

def randomMotifs(Dnas, k):
    """
    Given Dnas and k, return a randomly chosen Motifs.
    That is, each kmer in Motifs is randomly chosen from
    the corresponding DNA string in Dnas.
    """
    return [randomKmer(dna, k) for dna in Dnas]

