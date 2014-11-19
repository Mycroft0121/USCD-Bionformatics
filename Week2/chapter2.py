#create codon table as a dictionary keys: RNA 3-mers, value: Protein in one letter 
RNA_Codon ={}
for line in open("RNA_codon_table_1.txt"):
    RNA_Codon[line[0:3]]=line[4]

#create AminoAcidMass table as a dictionary, keys: Protein in one letter, value: mass integer 
AminoAcidMass = {}
for line in open("integer_mass_table.txt"):
    AminoAcidMass[line[0]]=int(line[2:5])

#input a string of RNA, output: translated protein sequence, break if seen stop codon
def RNAtoAminoAcid(text):
    num = len(text)/3
    aminoacid = "" 
    for i in range(num):
        if RNA_Codon[text[3*i:3*i+3]]!="\n":
            aminoacid += RNA_Codon[text[3*i:3*i+3]]
        else:
            break
    return aminoacid
#input: DNA seq in string, output: reverse complement seq in string    
def reversecomplement(pattern):
    complement = {"A":"T", "C":"G","G":"C","T":"A"}
    revcom = ""
    for char in pattern:
        revcom +=complement[char]
    return revcom[::-1]    
    
#input: DNA seq in string, amino acid sequence in string; output: substring of input DNA if that piece of DNA or its reverse complement encodes input aa seq   
def FindGeneFromAminoAcid(text, aminoacid):
    text1 = text.replace("T", "U")
    text2 = text1[1::]
    text3 = text1[2::]
    def FindGene(seq, aa):
        k = len(aminoacid)*3
        genes = []
        for i in xrange((len(seq)-k+1)/3):
            gene = seq[3*i:3*i+k]
            generc = reversecomplement(gene.replace("U", "T")).replace("T", "U")
            if RNAtoAminoAcid(gene)==aa:
                genes.append(gene)
            if RNAtoAminoAcid(generc)==aa:
                genes.append(gene)
        return genes
    totalSeq = [text1, text2, text3]
    Genes = []
    for ele in totalSeq:
        for gene in FindGene(ele, aminoacid):
            Genes.append(gene.replace("U", "T"))
    return len(Genes)     
 
#input: a string of aa seq, output: its molecular mass in integer       
def PeptideMass(peptide):
    return sum(AminoAcidMass[x] for x in peptide)
    
#input: a string or list of aa mass, output: its molecular mass in integer       
def PeptideMassIntString(peptide):
    if type(peptide)==list:
        return sum(x for x in peptide)
    else:
        peptide = peptide.replace("-", " ")
        peptide = peptide.split()
        peptide = [int(x) for x in peptide]
    return sum(x for x in peptide)    
              
#input: a string of aa seq, output: its full spectrum in string
def LinearSpectrum(Peptide):
    PrefixMass = [0]*(len(Peptide)+1)   
    for i in xrange(1,len(Peptide)+1):
        PrefixMass[i] = PeptideMass(Peptide[0:i])
    LinearSpectrum = [0]
    for i in xrange(len(Peptide)):
        for j in xrange(i+1,len(Peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    spectrum = sorted(LinearSpectrum)
    Spectrum = " ".join(str(x)for x in spectrum)
    return Spectrum
    
#input: a string of aa mass, output: its full spectrum in string
def LinearSpectrumIntString(peptide):
    peptide = peptide.replace("-", " ")
    peptide = peptide.split()
    peptide = [int(x) for x in peptide]
    PrefixMass = [0]*(len(peptide)+1)   
    for i in xrange(1,len(peptide)+1):
        PrefixMass[i] = PeptideMassIntString(peptide[0:i])
    LinearSpectrum = [0]
    for i in xrange(len(peptide)):
        for j in xrange(i+1,len(peptide)+1):
            LinearSpectrum.append(PrefixMass[j]-PrefixMass[i])
    spectrum = sorted(LinearSpectrum)
    Spectrum = " ".join(str(x)for x in spectrum)
    return Spectrum    
    
#input: a string of aa seq, output: its full spectrum in string
def CyclicSpectrum(Peptide):
    PrefixMass = [0]*(len(Peptide)+1) 
    for i in xrange(1,len(Peptide)+1):
        PrefixMass[i]= PeptideMass(Peptide[0:i])
    peptideMass = PrefixMass[len(Peptide)]
    CyclicSpectrum = [0]
    for i in xrange(len(Peptide)):
        for j in xrange(i+1,len(Peptide)+1):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i>0 and j < len(Peptide):
                CyclicSpectrum.append(peptideMass-PrefixMass[j]+PrefixMass[i])
    spectrum = sorted(CyclicSpectrum)
    Spectrum = " ".join(str(x)for x in spectrum)
    return Spectrum
    
#input: a string of aa mass, output: its full spectrum in string
def CyclicSpectrumIntString(peptide):
    peptide = peptide.replace("-", " ")
    peptide = peptide.split()
    peptide = [int(x) for x in peptide]
    PrefixMass = [0]*(len(peptide)+1) 
    for i in xrange(1,len(peptide)+1):
        PrefixMass[i]= PeptideMassIntString(peptide[0:i])
    peptideMass = PrefixMass[len(peptide)]
    CyclicSpectrum = [0]
    for i in xrange(len(peptide)):
        for j in xrange(i+1,len(peptide)+1):
            CyclicSpectrum.append(PrefixMass[j]-PrefixMass[i])
            if i>0 and j < len(peptide):
                CyclicSpectrum.append(peptideMass-PrefixMass[j]+PrefixMass[i])
    spectrum = sorted(CyclicSpectrum)
    Spectrum = " ".join(str(x)for x in spectrum)
    return Spectrum    
    
#input: a potential sublist and a target list. output: Boolean: whether the potential list is a sublist of the target list
def isSublist(sublist, targetlist):    
    sublistcc = [x for x in sublist]
    for ele in sublistcc:       
        try:
            targetlist.remove(ele)
            sublist.remove(ele)
        except ValueError:
            return False
        
    if len(sublist)==0:
        return True

     
#input: a full cyclospectrum in string, output: possible aa seq in string with duplicates   
def CyclopeptideSequencing(Spectrum):
    result = []
    peptides = [""]
    parentMass = int(Spectrum.split()[-1])
    AminoAcid = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    AminoAcidc = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    listspec = [int(i) for i in Spectrum.split()]
    for aa in AminoAcid:
        if AminoAcidMass[aa] not in listspec:
            AminoAcidc.remove(aa)
    def expand(peptides):
        length = len(peptides)
        oldpep = [x for x in peptides]
        peptides = []               
        for x in xrange(length):
            for aa in AminoAcidc:
                peptides.append(oldpep[x]+aa)      
        return list(set(peptides))
    while len(peptides)!=0:
        peptides = expand(peptides)
        peptidescc = [x for x in peptides]
        for peptide in peptidescc:

            mass = PeptideMass(peptide)
            if int(mass) == int(parentMass):
                peptides.remove(peptide)
                if CyclicSpectrum(peptide)==Spectrum:   
                    result.append(peptide)
            elif not isSublist(LinearSpectrum(peptide).split(), Spectrum.split()):
                peptides.remove(peptide)
    
    return list(set(result))

#convert result from CyclopeptideSequencing from a list of string of aa seq to a list of aa mass              
def CycSeqRes(Spectrum):
    aaseq = CyclopeptideSequencing(Spectrum)
    aamass = []
    for seq in aaseq:
        mass = [AminoAcidMass[x] for x in seq]
        mass = "-".join(str(x)for x in mass)
        aamass.append(mass)
    aamass = list(set(aamass))
    aamass = " ".join(str(x)for x in aamass)
    return aamass

#a peptite in string and a spectrum in string, output: integer score describing how overlapping the theoretical and experimental spectrum are              
def cyclopepScore(peptide, spectrum):
    theorectic = CyclicSpectrum(peptide).split()
    theorecticc = [x for x in theorectic]
    spectrum = spectrum.split()
    spectrumcc = [x for x in spectrum]
    score = 0
    for ele in theorectic:
        try:
            theorecticc.remove(ele)
            spectrumcc.remove(ele)
            score +=1
        except ValueError:
            score +=0
    return score

#input: a aa mass in string and a spectrum in string, output: integer score describing how overlapping the theoretical and experimental spectrum are              
def cyclopepScoreIntString(peptide, spectrum):
    theorectic = CyclicSpectrumIntString(peptide).split()
    theorecticc = [x for x in theorectic]
    spectrum = spectrum.split()
    spectrumcc = [x for x in spectrum]
    score = 0
    for ele in theorectic:
        try:
            theorecticc.remove(ele)
            spectrumcc.remove(ele)
            score +=1
        except ValueError:
            score +=0
    return score

#a peptite in string and a spectrum in string, output: integer score describing how overlapping the theoretical and experimental linear spectrum are              
def linearScore(peptide, spectrum):
    theorectic = LinearSpectrum(peptide).split()
    theorecticc = [x for x in theorectic]
    spectrum = spectrum.split()
    spectrumcc = [x for x in spectrum]
    score = 0
    for ele in theorectic:
        try:
            theorecticc.remove(ele)
            spectrumcc.remove(ele)
            score +=1
        except ValueError:
            score +=0
    return score
    
#a peptide mass seq in string and a spectrum in string, output: integer score describing how overlapping the theoretical and experimental linear spectrum are              
def linearScoreIntString(peptide, spectrum):
    theorectic = LinearSpectrumIntString(peptide).split()
    theorecticc = [x for x in theorectic]
    spectrum = spectrum.split()
    spectrumcc = [x for x in spectrum]
    score = 0
    for ele in theorectic:
        try:
            theorecticc.remove(ele)
            spectrumcc.remove(ele)
            score +=1
        except ValueError:
            score +=0
    return score

#input: leaderboard: a list of current high score peptides, N: integer of desired peptides
#spectrum: a string of experimental spectrum, 
#output: leaderboard: a list of N highest score peptides
def trim(Leaderboard, Spectrum, N):
    linearscore = []
    for peptide in Leaderboard:
        score = linearScore(peptide, Spectrum)
        linearscore.append([score, peptide])
    linearscore = sorted(linearscore)[::-1]
    for i in xrange(N,len(linearscore)):
        if linearscore[i][0] < linearscore[N-1][0]:
            linearscore = linearscore[0:i]
            break
    Leaderboard = [x[1] for x in linearscore]
    return Leaderboard 
                                       
#input: leaderboard: a list of current high score peptides mass string, N: integer of desired peptides
#spectrum: a string of experimental spectrum, 
#output: leaderboard: a list of N highest score peptides
def trimIntString(Leaderboard, Spectrum, N):
    linearscore = []
    for peptide in Leaderboard:
        score = linearScoreIntString(peptide, Spectrum)
        linearscore.append([score, peptide])
    linearscore = sorted(linearscore)[::-1]
    for i in xrange(N,len(linearscore)):
        if linearscore[i][0] < linearscore[N-1][0]:
            linearscore = linearscore[0:i]
            break
    Leaderboard = [x[1] for x in linearscore]
    return Leaderboard
                                                                                                            
#input a string of spectrum and a integer N, output N aa peptides strings with highest score
def leaderboardCyclopeptideSequencing(Spectrum, N):
    leaderboard = [""]
    leaderpeptide = ""
    parentMass = int(Spectrum.split()[-1])
    AminoAcid = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
    def expand(peptides):
        length = len(peptides)
        oldpep = [x for x in peptides]
        peptides = []               
        for x in xrange(length):
            for aa in AminoAcid:
                peptides.append(oldpep[x]+aa)      
        return list(set(peptides))
    while len(leaderboard):
        leaderboard = expand(leaderboard)
        leaderboardcc = [x for x in leaderboard]
        for peptide in leaderboardcc:
            mass = PeptideMass(peptide)
            if mass == parentMass:
                 if cyclopepScore(peptide, Spectrum) > cyclopepScore(leaderpeptide, Spectrum):
                    leaderpeptide =  peptide
            elif mass > parentMass:
                leaderboard.remove(peptide)
        leaderboard = trim(leaderboard, Spectrum, N)
    mass = [AminoAcidMass[x] for x in leaderpeptide]
    mass = "-".join(str(x)for x in mass)
    return mass

#Input: A collection of integers Spectrum(string).
#Output: The list of elements in the convolution of Spectrum. If an element has
#multiplicity k, it should appear exactly k times; you may return the elements in any order
def convolution(Spectrum):
    Spectrum = Spectrum.split()
    Spectrum = sorted([int(x) for x in Spectrum])
    convolution = []
    length = len(Spectrum)
    for i in xrange(length):
        for j in xrange(i+1,length):
            val = int(Spectrum[j])-int(Spectrum[i])
            if val >0:
                convolution.append(val)
    convolution = " ".join(str(x) for x in convolution)
    return convolution

#Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum in string.
# Output: A cyclic peptide LeaderPeptide with amino acids taken only from the top M elements
# (and ties) of the convolution of Spectrum that fall between 57 and 200, and where the size
# of Leaderboard is restricted to the top N (and ties).
def ConvolutionCyclopeptideSequencing(M, N, Spectrum):
    leaderboard = [""]
    leaderpeptide = ""
    Spectrum = Spectrum.split()
    Spectrum = sorted([int(x) for x in Spectrum ])
    parentMass =Spectrum[-1]
    Spectrum = " ".join(str(x) for x in Spectrum)
    convol = convolution(Spectrum)
    convol = convol.split()
    convol = [int(x) for x in convol if int(x)>=57 and int(x)<= 200]
    convoldict = {x:0 for x in set(convol)}
    for x in convol:
        convoldict[x] +=1
    convollist = sorted([count, mass] for mass, count in convoldict.items())[::-1]
    for i in xrange(M,len(convollist)):
        if convollist[i][0] < convollist[M-1][0]:
            convollist = convollist[0:i]
            break
    aminoacidmass = [int(x[1])for x in convollist]
    def expandIntString(peptides):
        length = len(peptides)
        oldpep = [x for x in peptides]
        peptides = []               
        for x in xrange(length):
            for aa in aminoacidmass:
                if len(oldpep[x]):
                    peptides.append(oldpep[x]+"-"+str(aa))
                else:
                    peptides.append(str(aa))      
        return list(set(peptides))
    while len(leaderboard):
        leaderboard = expandIntString(leaderboard)
        leaderboardcc = [x for x in leaderboard]
        for peptide in leaderboardcc:
            mass = PeptideMassIntString(peptide)            
            if mass == parentMass:
                 if cyclopepScoreIntString(peptide, Spectrum) > cyclopepScoreIntString(leaderpeptide, Spectrum):
                    leaderpeptide =  peptide
            elif mass > parentMass:
                leaderboard.remove(peptide)
        leaderboard = trimIntString(leaderboard, Spectrum, N)
    mass = leaderpeptide
    return mass 