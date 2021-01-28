# Alignment
Implement alignment algorithms 


- [Background](#background)
- [Usage](#usage)
- [API](#api)
- [Example](#example)

## Background
In this assignment, I will implement two classical alignment algorithms and then evaluate each algorithm’s performance with a range of parameters. 

## Usage 
Common alignment methos (global alignment and local alignment) for comparing the similarity of two sequences and computed the similarity scores. 

## API
* Class PairwiseAligner(seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty)
* Class NeedlemanWunsch(seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty)
* Class SmithWaterman(seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty)

For algorthims initialization, we need the orignal sequence, the score_matrix, and defining the open gap penalty and extension penalty.
For excution, the align matrix is computed. And the "flow" when computing the alignment matrix was recorded. When needed to calculated the alignment score, the right coner (for NeedlemanWunsh) or the peak value (for SmithWaterman) was found and the alignment sequence was formed along the "flow"
For termination, once the point reach the left or top broader ( for NeedlemanWunsh) ; zero value shown up ( for SmithWaterman), the termination happens.

The initialzation was the same but the exuation and termination are defined differently.

Affined-gap based will punish the gap nonlinealy corelated with the length. For implmentation ,everytime we need to check if the gap has been open or not.


"""
A common alignment class for comparing the similarity of two sequences.

### Parameters:
   
seq1 : str, default = NA,
        The dir of sequence A
        
seq2 : str, default = NA,
        The dir of sequence B

mtx_dir : str, default = NA,
        The dir of score matrix
        
gap_open_penalty : int, default = NA,
        The value of gap open penalty
        
gap_extension_penalty, int, default = NA.
        The value of gap extension penalty
        
### Attributes

gop : int gap_open_penalty; 

gep : int gap_extension_penalty; 

sequence_1/sequence_2 str, the content of sequence A/B; 

ID_1/ID_2, the IDs of sequence A/B; 

score_matrix : pandas dataframe, the reference score matrix for amino acid comparision; 

align_mtx, pandas dataframe, the computed dataframe for score calculation and sequence alignment; 

align_score, int, the computed number representing the similarity of two sequences; 

align_optimal_track, list, a list of location (representing by [index, columns] used for computing the alignment sequences; 

align_seq1/align_seq2, str, the aligned sequences of sequence A/B; 

trace_mtx, pandas dataframe, it is a matrix representing the "flow" during the computation. When during the sequence alignment, once the start point was found, I can just trace back to find the whole sequence.

![alt text](https://i.stack.imgur.com/LDiz2.jpg)

As shown in this figure. align_mtx is the values dataframe. trace_mtx is the dataframe representing the arrows. Eg. For the df['T']['T']  the align_mtx is 6 align_mtx.iloc[2,3] =6
the trace_mtx of [2,3] is connected with [1.2], and also is linked to [2,4], [3,4], [3,3]. so trace_mtx.iloc[2,3] = [ [1,2] , [ [2,4], [3,4] , [3,3]]] 

### Methods

* Initialization

	The basic alignment class.
	The parent class of  NeedlemanWunsch (Global alignment) and SmithWaterman (Local Alignment) class.
   If you derectly use Pairwise class, the default method is global alignment. 

        Align = SmithWaterman('../sequences/prot-0091.fa', '../sequences/prot-0093.fa','../scoring_matrices/BLOSUM62.mat',2,1)
   Anyone of the three can be used in this way
   
* Calculating the alignment score

    Computing scores became eaiser with the help of track matrix, for global alignment, just starting from the right down conor, and always return to the previous location. And at the same time, the location was added into align_optimal_track for the further computing of the align sequence.
    ps: always run this function first prior to get the alignment scores and alignment sequences.
   

        Align.compute_alignment_score()
  
* Return the alignment sequence

        Align.get_alignment_sequence()
        
* Representing the alignment score
        
        Align.get_alignment_score()
        
*  penalty(i,j)

    Since we are doing afine gap penalty, then everytime when we want to have a gap we need to justify whether we should give a open penalty or externsion penalty. To decide on this, my method is to judge whether the previous point already began the gap or not. For this penalty funtion, [i,j] is the location of either the upper location or the left location. We judge whether this previous location is 1) linked to its left or upper 2) linked to the upper left (diagnoal) or don't have links. if 1), give the following one extension penalty. else, (Just started the sequence or the gap did not occur) give the open penalty
        
* alignment()

    This function is packaged in other function. No need to use it.
    Initialize align_mtx

* align()

    This function is packaged in other function. No need to use it.
    Compute the aligned sequences

## Example

In this turtorial, we will together initialize a global alignment class and a local alignment class, compute the alignment scores and alignment matrix.
Go to the folder test, open a jupyter notebook. ( Or other editors...)

* Initialization

Load packages

		import sys
		import numpy as np
		import pandas as pd
		sys.path.append("..")
		import align.algs as algs
		import os

dir1 is a fasta file, dir2 is the same fasta file but with a huge gap, dir 3 is the scoreing matrix

		dir1 = 'Test_sequences/Sequence-1.fa'
		dir2 = 'Test_sequences/Sequence-1-paired-GAP.fa'
		dir3 = '../scoring_matrices/PAM250.mat'

Initialize the class and calculate the alignment scores. (At the same time, the align_sequence is also computed.

		Align_global = algs.NeedlemanWunsch(dir1,dir2,dir3,11,3)
		Align_local  = algs.SmithWaterman(dir1,dir2,dir3,11,3)
		Align_global.compute_alignment_score()
		Align_local.compute_alignment_score()
		
Show the alignment scores:

		Align_global.get_alignment_score()
		
![alt text](https://user-images.githubusercontent.com/49840852/106100881-e75dc880-60f1-11eb-8c4b-87f5948714db.png)

		Align_local.get_alignment_score()
		
![alt text](https://user-images.githubusercontent.com/49840852/106100882-e75dc880-60f1-11eb-87b9-a850f2d3a6a0.png)		
		
Show the alignment sequences:

		Align_global.get_alignment_sequence()

![alt text](https://user-images.githubusercontent.com/49840852/106100876-e5940500-60f1-11eb-9646-b7fa02fc4335.png)

		Align_local.get_alignment_sequence()
		
![alt text](https://user-images.githubusercontent.com/49840852/106100880-e6c53200-60f1-11eb-93dc-df008b2746e8.png)
	
	




