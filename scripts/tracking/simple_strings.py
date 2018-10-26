# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#      http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys 
import re

# Min. number of letters to be dropped at the end of a word
MIN_ABBREV = 2
# Min. number of letters to be skipped within a word
MIN_SKIP = 2
# If True prints out a bunch of debugging message along the way
DEBUG = True

# Compute length of overlap between sequences (strings or sequences of strings)
# Arg1: string or sequence of strings
# Arg2: string or sequence of strings
# Return type: int
def overlap(seq1, seq2):
	for i in range(min(len(seq1), len(seq2))):
		if seq1[i].upper() != seq2[i].upper():
			return i
	return min(len(seq1), len(seq2))

# Tests if all elements in string are characters 
# Arg1: string
# Return type: bool
def only_characters(seq):
	char_match = re.match('[a-zA-Z]+', seq)
	return char_match is not None and len(char_match.group(0)) == len(seq)

# Propose shortening of listed strings
# Arg1: list of strings
# Return type: dict(original string => shortened string)
def simplify_words(seqs):
	# Compute sorted list of unique input strings
	sorted_seqs = sorted(list(set(seqs)))
	if len(sorted_seqs) == 0:
		return dict()

	max_length = max(len(seq) for seq in sorted_seqs)
	# Final seqs to be returned
	new_seqs = list(sorted_seqs)

	# Compute length of overlap between each string and the previous one in the list
	overlaps = [0] + [overlap(X, Y) for X, Y in zip(sorted_seqs[1:], sorted_seqs[:-1])] + [0]

	for index in range(len(new_seqs)):
		# Remove any non-informative suffix by comparing to seqence to neigbhours in list
		non_unique_length = max(overlaps[index], overlaps[index + 1]) + 1
		for cutpoint in range(non_unique_length, len(sorted_seqs[index]) - MIN_ABBREV):
			# Do no cut on a linking character
			if sorted_seqs[index][cutpoint-1] not in '.-+':
				new_seqs[index] = sorted_seqs[index][:cutpoint] + '.'
				break

		# Remove non-informative middle bits:
		long_overlap = overlaps[index]
		# If you detect a previous overlaps which is longer than the next one
		while long_overlap > overlaps[index + 1] + MIN_SKIP:
			# See how many sequences in the list share that overlap
			for index2 in range(index, -1, -1):
				if overlaps[index2] < long_overlap:
					break
			# How long is their overlap
			common_overlap = min(overlaps[index2+1:index+1])
			# How much do they overlap with the two sequences before and after
			max_boundary_overlap = max(overlaps[index2], overlaps[index+1])
			# Remember for next round
			long_overlap = max_boundary_overlap
			# How many chars do you need to uniquely recognise the seqs shraing the prefix
			min_unique = max_boundary_overlap + 1
			# Check that the difference is large enough to justify abbreviation
			# Check that the skip only contains alphabet characters, not digits or funny chars
			if common_overlap - min_unique > MIN_SKIP and only_characters(new_seqs[index][min_unique:common_overlap]):
				for i in range(index2, index + 1):
					if len(new_seqs[i][common_overlap:]) > 0:
						new_seqs[i] = new_seqs[i][:min_unique] + "'" + new_seqs[i][common_overlap:]
					else:
						new_seqs[i] = new_seqs[i][:min_unique] + "." 

	res = dict(zip(sorted_seqs, new_seqs))
	
	if DEBUG:
		for overlp, sorted_seq, new_seq in zip(overlaps, sorted_seqs, new_seqs):
			print 'WORD', overlp, sorted_seq, '=>', new_seq

	return res

# Propose shortening of listed list of strings
# Arg1: list of list of strings
# Return type: dict(original list of string => shortened list of strings)
def simplify_phrases(phrases):
	# Sort unique input lists
	sorted_phrases = sorted(set(map(tuple, phrases)))
	if len(sorted_phrases) == 0:
		return dict()

	max_length = max(len(phrase) for phrase in sorted_phrases)
	new_phrases = map(list, sorted_phrases)

	# Compute overlap between each list and its next neighbour
	overlaps = [0] + [overlap(X, Y) for X, Y in zip(sorted_phrases[1:], sorted_phrases[:-1])] + [0]
	

	for index in range(len(new_phrases)):
		# Remove any non-informative suffix
		non_unique_length = max(overlaps[index], overlaps[index + 1]) + 1
		if non_unique_length < len(sorted_phrases[index]):
			new_phrases[index] = sorted_phrases[index][:non_unique_length] 

		# Remove non-informative middle bits:
		long_overlap = overlaps[index]
		# If you detect a previous overlaps which is longer than the next one
		while long_overlap > overlaps[index + 1] + 1:
			# See how many sequences in the list share that overlap
			for index2 in range(index, -1, -1):
				if overlaps[index2] < long_overlap:
					break
			# How long is their overlap
			common_overlap = min(overlaps[index2+1:index+1])
			# How much do they overlap with the two sequences before and after
			max_boundary_overlap = max(overlaps[index2], overlaps[index+1])
			# How many chars do you need to uniquely recognise the seqs shraing the prefix
			min_unique = max_boundary_overlap + 1
			# Check that the difference is large enough to justify abbreviation
			if common_overlap > min_unique:
				for i in range(index2, index + 1):
					new_phrases[i] = new_phrases[i][:min_unique] + new_phrases[i][common_overlap:]

			# Remember for next round
			long_overlap = max_boundary_overlap

	res = dict(zip(sorted_phrases, new_phrases))

	if DEBUG:
		for overlp, sorted_phrase, new_phrase in zip(overlaps, sorted_phrases, new_phrases):
			print 'PHRASE', overlp, " ".join(sorted_phrase), '=>', " ".join(new_phrase)

	return res

if __name__ == "__main__":
	# Read input
	seqs = [line.strip() for line in sys.stdin]

	# Decompose into lists of words
	phrases = [re.split(r'[\s,:;"()\[\]]+', seq) for seq in seqs]

	# Collect dictionary of words and propose unique abbreviations
	words = set.union(*[set(phrase) for phrase in phrases])
	simple_words = simplify_words(words)

	# Shorten the words 
	abbrev_phrases = [[simple_words[word] for word in phrase if len(word) > 0 ] for phrase in phrases]

	# Drop unecessary words
	simple_phrases = simplify_phrases(abbrev_phrases)
	result = [" ".join(simple_phrases[tuple(phrase)]) for phrase in abbrev_phrases]

	# Output
	for seq, res in zip(seqs, result):
		print seq, ">>", res

	if DEBUG:
		print "Avg. Input length: %f; Avg. output length: %f" % (1.0 * sum(len(seq) for seq in seqs) / len(seqs), 1.0 * sum(len(res) for res in result) / len(result))
		print "Max. input length: %i; Max. output length: %i" % (max(len(seq) for seq in seqs), max(len(res) for res in result))
		print "Max. length output sequences: "
		max_length = max(len(res) for res in result)
		print "\n".join(res for res in result if len(res) ==  max_length)
