import random
import copy
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = copy.deepcopy(Motifs)
    for j in range(N):
      i = random.randint(0, t-1)
      Motifs.pop(i)
      Profile = ProfileWithPseudocounts(Motifs)
      Motifi = ProfileGeneratedString(Dna[i], Profile, k)
      Motifs.insert(i,Motifi)
      if Score(Motifs) < Score(BestMotifs):
        BestMotifs = Motifs.copy()
    return BestMotifs

# place all subroutines needed for GibbsSampler below this line
def RandomMotifs(Dna, k, t):
  """Input is a dictionary of motifs (Dna), k is length of kmers and t is the 
     length of Dna. Output is a list of randomly generated kmers
     from each list in Dna""" 
  motifs = []
  t = len(Dna)
  j = len(Dna[0]) 
  for i in range(t):
    rand = random.randint(0, j-k)
    motifs.append(Dna[i][rand:rand+k])
  return motifs

def CountWithPseudocounts(Motifs):
  """Takes in a a list of kmers and outputs a dictionary of bases and with counts
     for each base in each location of the kmer"""
  t = len(Motifs)
  k = len(Motifs[0])
  count = {} # output variable
  k = len(Motifs[0])
  for symbol in "ACGT":
      count[symbol] = []
      for j in range(k):
              count[symbol].append(1)                
  t = len(Motifs)
  for i in range(t):
      for j in range(k):
          symbol = Motifs[i][j]
          count[symbol][j] += 1
  return count

def ProfileWithPseudocounts(Motifs):
  """Takes in a list of Dna motifs and outputs a dictionary profile with the
     percentage that the base pair appears in the index location"""
  t = len(Motifs)
  k = len(Motifs[0])
  profile = CountWithPseudocounts(Motifs) # returns a count profile as a dictionary
  for key in profile:
      for value in range(k):
          profile[key][value] = profile[key][value] / (t + 4) 
  return profile

def Normalize(Probabilities):
  """Input is probability dictionary and the output is a normalized dictionary"""
  d = {}
  sumvals = sum(Probabilities.values())
  for k,v in Probabilities.items():
    d[k] = Probabilities[k]/sumvals
  return d

def WeightedDie(Probabilities):
  """Input is a probability dictionary and the output is a randomly 
     generated kmer"""
  kmer = '' # output variable
  count = 0
  population = list(Probabilities.keys())
  weight = list(Probabilities.values())
  for i in range(1):
    kmer += random.choices(population, weights=weight, k=1)[0]
  return kmer


def ProfileGeneratedString(Text, profile, k):
  """Input is a string Text, a profile matrix profile, and an integer k and the
     output is a randomly generated kmer"""
  n = len(Text)
  probabilities = {}
  for i in range(0, n-k+1):
    probabilities[Text[i:i+k]] = pr(Text[i:i+k],profile)
  probabilities = Normalize(probabilities)
  return WeightedDie(probabilities)

def Score(Motifs):
  """Takes in a list of motifs and scores each motif combined against the
     consensus generated motif. Each time the base does not match the one in 
     in the consensus it increases the count by 1. The lower the score the 
     better"""
  count = 0
  consensus = Consensus(Motifs)
  t = len(Motifs)
  k = len(Motifs[0])
      # initialize loop to go through each row(list) of Motifs
  for i in range(t):
      # initialize loop to go through each index of each list
    for j in range(k):
      # use an if statement to check if character is the same character of the same index of consensus
      if Motifs[i][j] != consensus[j]:
        count += 1
  return count

def Consensus(Motifs):
  """Takes in a list of motifs and outputs a single motif string with the max
     value at each index"""
  k = len(Motifs[0])
  count = CountWithPseudocounts(Motifs)
  consensus = ""
  for j in range(k):
      m = 0
      frequentSymbol = ""
      for symbol in "ACGT":
          if count[symbol][j] > m:
              m = count[symbol][j]
              frequentSymbol = symbol
      consensus += frequentSymbol
  return consensus

def pr(gene, Profile):
  """Input is a string (gene) and a dictionary (Profile). Output is float type
     probability"""
  probability = 1
  for i, char in enumerate(gene):
    probability *= Profile[char][i]
  return probability
