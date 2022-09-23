#!/usr/bin/env python3
# USAGE: my-needleman-wunsch.py <fasta_file> <scoring_matrix> <gap_penalty> <output_filepath>

import numpy as np
import sys
from fasta import readFASTA

# read in sequences
input_sequences = readFASTA(open(sys.argv[1]))

seq1_id, sequence1 = input_sequences[0]
seq2_id, sequence2 = input_sequences[1]

# read in gap penalty and make it negative if needed
gap_penalty = float(sys.argv[3])
if gap_penalty >= 0:
	gap_penalty = -1 * gap_penalty

# print(gap_penalty)

#=======================#
# create scoring matrix #
#=======================#

# initialize dictionary to store row/col positions of nucleotides/AA
score_dict = {}

# initialize list of lists
list_of_rows = []

# to differentiate header row
row_count = 0

for row in open(sys.argv[2]):
	if row_count == 0:

		# store first row entries in dictionary
		for i, entry in enumerate(row.split()):
			score_dict[entry] = i
		row_count +=1 

	else:
		
		# initialize list to go in list of lists
		row_list = []
		for entry in row.split(): 
			
			# if the entry in the line is a number (not the row header)
			try: 
				entry = int(entry)

				# add it to the list
				row_list.append(entry)

			# otherwise leave it out
			except:
				pass

		# add the row to the list of lists
		list_of_rows.append(row_list)

# make list of lists into a numpy array
sm = np.array(list_of_rows)

# print(sm)

#=====================#
# initialize matrices #
#=====================#

# intialize F-matrix
Fm = np.zeros((len(sequence1) + 1, len(sequence2) + 1))

# initialize first column/row with gap penalty
for i in range(Fm.shape[0]):
	Fm[i,0] = i * gap_penalty

for j in range(Fm.shape[1]):
	Fm[0,j] = j * gap_penalty

# initialize traceback matrix
tbm = np.zeros((len(sequence1) + 1, len(sequence2) + 1))

# numerical codes for arrows
diag = 0
left = 1
up = 2

# set all first row to left arrows
tbm[0,1:] = left

# set all first col to up arrows
tbm[1:,0] = up

# print(Fm)
# print(tbm)

#===================#
# populate matrices #
#===================#

for i in range(1, Fm.shape[0]):
	for j in range(1, Fm.shape[1]):

		# determine match or mismatch score using the scoring dictionary
		match_mismatch = sm[score_dict[sequence1[i-1]],score_dict[sequence2[j-1]]]

		# calculate d, h and v
		d = Fm[i-1,j-1] + match_mismatch
		h = Fm[i,j-1] + gap_penalty
		v = Fm[i-1,j] + gap_penalty

		# calculate F
		F = max(d,h,v)

		# set F-matrix entry to F
		Fm[i,j] = F

		# set arrow in traceback matrix; if there is a tie, should set F to
		# d, then h, then v
		if F == v:
			tbm[i,j] = up
		if F == h:
			tbm[i,j] = left
		if F == d:
			tbm[i,j] = diag

# print(Fm)
# print(tbm)

#========================#
# find optimal alignment #
#========================#

# set starting position for while loop
i = len(sequence1)
j = len(sequence2)

# record the alignment score
align_score = Fm[i,j]

# initialize sequence alignment strings
align_1 = ''
align_2 = ''

# initialize gap counts
gaps_1 = 0
gaps_2 = 0

# as long as we're not at position 0,0...
while i > 0 or j > 0:

	# add the letter for both sequences to the beginning of the string if they
	# are aligned
	if tbm[i,j] == diag:
		align_1 = sequence1[i-1] + align_1
		align_2 = sequence2[j-1] + align_2

		# move up diagnoally 
		i -= 1
		j -= 1

	# add a gap in sequence 1
	elif tbm[i,j] == left:
		align_1 = '-' + align_1
		align_2 = sequence2[j-1] + align_2

		# move to the left
		j -= 1

		# increment gap count
		gaps_1 += 1

	# add a gap in sequence 2
	elif tbm[i,j] == up:
		align_1 = sequence1[i-1] + align_1
		align_2 = '-' + align_2

		# move to up
		i -= 1

		# increment gap count
		gaps_2 += 1

# Write the alignment to a file
file = open(sys.argv[4], 'w')
file.write('Sequence 1 alignment: ' + align_1 + '\n' + '\n')
file.write('Sequence 2 alignment: ' + align_2)
file.close()

# Print gap counts and alignment score
print('Gaps in sequence 1: ' + str(gaps_1))
print('Gaps in sequence 2: ' + str(gaps_2))
print('Final alignment score: ' + str(align_score))


