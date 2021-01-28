import pytest
import sys
import numpy as np
import pandas as pd
sys.path.append("..")
import align.algs as algs
import os
#@pytest.fixture
def some_relevant_data(dir1,dir2,dir3):
	'''
	I load the dir of two sequences, the dir of scores matrix, build a global alignment class and a local alignment class.
	'''
	seq1 = algs.fasta_reader(dir1)
	
	seq2 = algs.fasta_reader(dir2)
	mtx = algs.score_matrix(dir3)
	

	Align_global = algs.PairwiseAligner(dir1,dir2,dir3,11,3)
	Align_local  = algs.SmithWaterman(dir1,dir2,dir3,11,3)
	Align_global.compute_alignment_score()
	Align_local.compute_alignment_score()
	
	return seq1, seq2, mtx, Align_global, Align_local

dir1 = 'test/Test_sequences/Sequence-1.fa'
dir2 = 'test/Test_sequences/Sequence-1-paired.fa'
dir3 = 'scoring_matrices/PAM250.mat'

seq1, seq2, mtx, Align_global, Align_local = some_relevant_data(dir1, dir2, dir3)


def test_fasta_io():
	'''
	make sure it is the correct sequence
	'''
	#seq1 = algs.fasta_reader('test/Test_sequences/Sequence-1.fa')
	assert seq1.sequence_as_string() =='SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM', 'The sequence was read incorrectly'

def test_scoring_matrix_io():
	'''
	Make sure the matrix has the correct max, min values and correct deminsion
	'''
	assert np.max(mtx.to_numpy()) == 17, 'The matrix was not loaded correctly'
	assert np.min(mtx.to_numpy()) == -8, 'The matrix was not loaded correctly'
	assert sum(np.diag(mtx.to_numpy())) == sum([2,6, 2, 4,12, 4, 4, 5, 6, 5, 6, 5, 6, 9, 6, 2, 3, 17, 10, 4, 3, 3, -1, 1]), 'The matrix was not loaded correctly'

def test_identical():
	'''
	Here I use two exactly same sequences for comparision, you can see the align sequences are identical to each other.
	For both the global alignment and local alignment

	Notes: There is also a sequences with a large gap in the test folder. If someone use that in the implementation, you will say a nice gap.
	'''
	assert Align_local.align_seq1 == Align_local.align_seq2, 'Same sequence should give same alignment results'

def test_alignment_score():
	'''
	Since the two sequences are same, with the same paramemters, the global and local alignment will give you the same score.
	'''
	assert Align_global.get_alignment_score() == 46267, 'The score was computed incorrectly!'
	assert Align_global.get_alignment_score() == Align_local.get_alignment_score(), "They are same sequences, so the score should be same for global and local alignment"


