import alg_alignment
import alg_project4_solution
import random
import time
import math
from matplotlib import pyplot

PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"


# question 1
hep = alg_alignment.read_protein(HUMAN_EYELESS_URL)
fep = alg_alignment.read_protein(FRUITFLY_EYELESS_URL)

scoring_matrix = alg_alignment.read_scoring_matrix(PAM50_URL)
hep_fep_local_alignment = alg_project4_solution.compute_local_alignment(hep, fep, scoring_matrix,
                                                                        alg_project4_solution.compute_alignment_matrix(
                                                                            hep, fep, scoring_matrix, False))
human_eyeless_fruitfly_local_alignment_score = hep_fep_local_alignment[0]
# question 1 answer
print "local alignment for human and fruitfly eyeless genome: " + str(hep_fep_local_alignment)

# question 2
cpd = alg_alignment.read_protein(CONSENSUS_PAX_URL)

hep_local_alignment = hep_fep_local_alignment[1]
fep_local_alignment = hep_fep_local_alignment[2]

hep_local_alignment_no_dashes = hep_local_alignment.replace('-', '')

hep_no_dashes_cpd_global_alignment = alg_project4_solution.compute_global_alignment(hep_local_alignment_no_dashes, cpd, scoring_matrix, alg_project4_solution.compute_alignment_matrix(
                                                   hep_local_alignment_no_dashes, cpd, scoring_matrix, True))
fep_local_alignment_no_dashes = fep_local_alignment.replace('-', '')

fep_no_dashes_cpd_global_alignment = alg_project4_solution.compute_global_alignment(fep_local_alignment_no_dashes, cpd, scoring_matrix,
                                                                                    alg_project4_solution.compute_alignment_matrix(fep_local_alignment_no_dashes, cpd, scoring_matrix, True))

print hep_no_dashes_cpd_global_alignment
print fep_no_dashes_cpd_global_alignment

# compute the percentage of elements in these two sequences that agree
hndga = hep_no_dashes_cpd_global_alignment[1]
hndgacpd = hep_no_dashes_cpd_global_alignment[2]


human_consensus_match_count = 0
for idx, elem in enumerate(hndga):
    if hndgacpd[idx] == elem:
        human_consensus_match_count += 1


fndcga = fep_no_dashes_cpd_global_alignment[1]
fndcgacpd = fep_no_dashes_cpd_global_alignment[2]

fly_consensus_match_count = 0
for idx, elem in enumerate(fndcga):
    if fndcgacpd[idx] == elem:
        fly_consensus_match_count += 1


# question 2 answer
print "human consensus match percentage: " + str(human_consensus_match_count) + "/" + str(len(hndga)) + " = " + str(human_consensus_match_count/float(len(hndga)))
print "fly consensus match percentage: " + str(fly_consensus_match_count) + "/" + str(len(fndcga)) + " = " + str(fly_consensus_match_count / float(len(fndcga)))


# question 3
# the agreement would be unlikely, it suggests that the sequences are related.

# question 4
def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    """
    return a dictionary scoring_distribution that represents an un-normalized distribution
    """
    distribution = {}
    for trial in range(num_trials):
        start = time.time()
        rand_y = list(seq_y)
        random.shuffle(rand_y)
        rand_y = ''.join(rand_y)
        alignment_matrix = alg_project4_solution.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        alignment = alg_project4_solution.compute_local_alignment(seq_x, rand_y, scoring_matrix, alignment_matrix)
        score = alignment[0]
        if score in distribution:
            distribution[score] += 1
        else:
            distribution[score] = 1
    return distribution

#score_distribution = generate_null_distribution(hep, fep, scoring_matrix, 1000)

#print score_distribution

distribution = {38: 1, 39: 4, 40: 8, 41: 16, 42: 20, 43: 25, 44: 38, 45: 48, 46: 66, 47: 70, 48: 68, 49: 67, 50: 68,
                51: 54, 52: 60, 53: 51, 54: 39, 55: 46, 56: 40, 57: 29, 58: 29, 59: 25, 60: 19, 61: 10, 62: 15, 63: 11,
                64: 13, 65: 6, 66: 4, 67: 9, 68: 7, 69: 7, 70: 6, 71: 4, 73: 4, 74: 3, 75: 1, 76: 3, 78: 1, 82: 1,
                84: 1, 87: 1, 88: 1, 90: 1}

# question 4 answer
weighted_distribution = dict(map(lambda (k, v): (k, v / float(1000)), distribution.iteritems()))

# expanded distribution makes stats easier
expanded_distribution = list()
for score in distribution:
    expanded_distribution = expanded_distribution +[score] * distribution[score]

print "expanded distribution: " + str(expanded_distribution)
#pyplot.bar(range(len(weighted_distribution)), weighted_distribution.values())
#pyplot.xticks(range(len(weighted_distribution)), weighted_distribution.keys())
#pyplot.ylabel('Frequency')
#pyplot.xlabel('Alignment score')
#pyplot.title('Distribution of alignment scores for  HumanEyelessProtein and randomized FruitflyEyelessProtein')
#pyplot.show()

# question 5
# calculate the mean of the distribution
print "num scores: " + str(len(expanded_distribution))
print "sum of scores: " + str(sum(expanded_distribution))
mean_score = sum(expanded_distribution) / float(len(expanded_distribution))

def get_squared_difference(value, mean):
    return math.pow(value - mean, 2)

print distribution.values()
sum_of_squared_differences = sum(map(lambda score: get_squared_difference(score, mean_score), expanded_distribution))
std_deviation = math.sqrt(sum_of_squared_differences / float(len(expanded_distribution)))
# question 5 answers
print "mean: " + str(mean_score)
print "std deviation: " + str(std_deviation)
print "human/fruitfly local alignment score: " + str(human_eyeless_fruitfly_local_alignment_score)
z_score = (human_eyeless_fruitfly_local_alignment_score - mean_score) / std_deviation
print "z score: " + str(z_score)
# question 6 answer
# The score resulting from the local alignment is unlikely due to chance.
# odds of winning powerball are ~ 1 in 300 million
# odds of z score of 113 being due to chance, i'm not even sure how to calculate this.  it is very very unlikely
# based on this wikipedia article: https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule, the frequency of being outside of a 7 std deviation range is something like in in 390 billion
# our local alignment score is about 16 standard deviations from the mean


# question 7

print "abc/abc: " + str(alg_project4_solution.edit_distance('abc', 'abc'))
print "a/b: " + str(alg_project4_solution.edit_distance('a', 'b'))
print "kitten/sitting: " + str(alg_project4_solution.edit_distance('kitten', 'sitting'))
print "abc/def: " + str(alg_project4_solution.edit_distance('abc', 'def'))
print "ryan/lauren: " + str(alg_project4_solution.edit_distance('ryan', 'lauren'))
print "snargle/bargle: " + str(alg_project4_solution.edit_distance('snargle', 'bargle'))

# question 7 answer:
# diagonal = 2, off diagonal = 1, dash = 0

# question 8
def check_spelling(checked_word, dist, word_list):
    """
    iterates through word_list and returns the set of all words that are within dist of the string checked_word
    """
    words_within_edit_distance = set()
    for word in word_list:
        word_distance = alg_project4_solution.edit_distance(checked_word, word)
        if word_distance <= dist:
            words_within_edit_distance.add(word)
    return words_within_edit_distance

word_list = alg_alignment.read_words(WORD_LIST_URL)
humble_list = check_spelling("humble", 1, word_list)
firefly_list = check_spelling("firefly", 2, word_list)
print "within 1 of humble: " + str(humble_list) + ", len: " + str(len(humble_list))
print "within 2 of firefly: " + str(firefly_list) + ", len: " + str(len(firefly_list))

# question 8 answer
#within 1 of humble: set(['bumble', 'humbled', 'tumble', 'humble', 'rumble', 'humbler', 'humbles', 'fumble', 'humbly', 'jumble', 'mumble']), len: 11
#within 2 of firefly: set(['firefly', 'tiredly', 'freely', 'fireclay', 'direly', 'finely', 'firstly', 'liefly', 'fixedly', 'refly', 'firmly']), len: 11


