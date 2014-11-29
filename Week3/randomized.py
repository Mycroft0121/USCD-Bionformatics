# You do not need to edit this for the lab.

# Randomized Motif Search.
# This is a solution to http://rosalind.info/problems/3f/
# However, we run only N=100 tries (rather than N=1000).
#
# Usage: python randomized.py INPUTFILE

import common

def RandomizedMotifSearch(Dnas, k,count):
    """
    Start with a random Motifs, then repeatedly update Profile
    and Motifs until the score does not improve.
    This function returns a pair: (bestMotifs, bestScore)
    """
    # Initial Motifs: pick each substring at random.
    # To save time, we are careful to compute each Score
    # only once. Returns a pair: (bestMotifs, bestScore)
    Motifs = common.randomMotifs(Dnas, k)
    bestMotifs = Motifs
    bestScore = common.Score(bestMotifs)
    while True:
        Profile = common.ProfileOfMotifs(Motifs, 1)
        Motifs = common.bestMotifsForProfile(Dnas, Profile)
        score = common.Score(Motifs)
        count +=1
        if score >= bestScore:
            return (bestMotifs, bestScore)
        else:
            bestMotifs = Motifs
            bestScore = score

def RepeatedRandomizedMotifSearch(Dnas, k, N):
    "Run RandomizedMotifSearch(Dnas, k) N times, and return best Motifs."
    # Some noisy print statements, so we see some progress
    count = 0
    print "Running RandomizedMotifSearch", N, "times"
    bestMotifs = common.randomMotifs(Dnas, k)
    bestScore = common.Score(bestMotifs)
    for i in range(1,N+1): # from 1 to N
        (Motifs, score) = RandomizedMotifSearch(Dnas, k, count)
        if score < bestScore:
            bestMotifs = Motifs
            bestScore = score
        # Print every score, and occasionally print bestScore
        print score,
        if i%10==0 or i==N:
            print "(after", i, "tries bestScore is", str(bestScore)+")"
    return bestMotifs, count

# If this module is the main program, run an example.
# We can optionally give the input filename on the comm