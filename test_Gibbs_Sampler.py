import Gibbs_Sampler
import pytest
from collections import Counter
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
k = 6
t = len(Dna)

@pytest.fixture
def motif():
    motifs = ['ATGTTC', 'TCAAGC', 'CAGTAC']
    return motifs

@pytest.fixture
def counts(motif):
    return Gibbs_Sampler.CountWithPseudocounts(motif)

@pytest.fixture
def profile(motif):
    return Gibbs_Sampler.ProfileWithPseudocounts(motif)

@pytest.fixture
def probabilities(profile):
    n = len(Dna[0])
    text = Dna[0]
    probability_dict = {}
    for i in range(0, n-k+1):
        probability_dict[text[i:i+k]] = Gibbs_Sampler.pr(text[i:i+k],profile)
    return probability_dict

@pytest.fixture
def normalize(probabilities):
    d = Gibbs_Sampler.Normalize(probabilities)
    return d

@pytest.fixture
def consensus(motif):
    return Gibbs_Sampler.Consensus(motif)

@pytest.fixture
def score(motif):
    return Gibbs_Sampler.Score(motif)

def test_randommotifs_generate_different_motifs():
    results = []
    for i in range(2):
        results.append(Gibbs_Sampler.RandomMotifs(Dna, k, t))
    assert results[0] != results[1]

def test_random_motifs():
    motifs = Gibbs_Sampler.RandomMotifs(Dna, k, t)
    for i in range(t):
        assert motifs[i] in Dna[i]

def test_count_with_psuedocounts(motif, counts):
    assert counts['A'] == [2, 2, 2, 2, 2, 1] and\
           counts['T'] == [2, 2, 1, 3, 2, 1] and\
           counts['C'] == [2, 2, 1, 1, 1, 4] and\
           counts['G'] == [1, 1, 3, 1, 2, 1] 

def test_profile_with_pseudocounts(motif, profile, counts):
    t = len(motif)
    k = len(motif[0])
    assert len(profile['A']) == len(motif[0])
    for key in profile:
        for value in range(k):
            assert profile[key][value] == counts[key][value]/(t+4)

def test_Normalize(normalize):
    total = round(sum(normalize.values()), 6)
    assert total == 1

def test_weighted_die(normalize):
    rolls = []
    values = normalize.values()
    max_val = max(values)
    for i in range(1000):
        rolls.append(Gibbs_Sampler.WeightedDie(normalize))
    rolls_dict = Counter(rolls)
    best_val = 0
    for roll in rolls_dict:
        if rolls_dict[roll] > best_val:
            best_val = rolls_dict[roll]
            most_probable = roll
    assert normalize[most_probable] == max_val
    

def test_pr(motif, profile):
    kmer = motif[0]
    k = len(kmer)
    pr = round(Gibbs_Sampler.pr(kmer, profile), 6)
    test_probability = round(1*(2*2*3*3*2*4)/(7**k), 6)
    assert pr == test_probability

def test_consensus(consensus):
    assert consensus == 'AAGTAC'

def test_score(score):
    assert score == 8

