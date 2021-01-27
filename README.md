# Alignment
Implement alignment algorithms 


- [Background](#background)
- [Usage](#usage)
- [API](#api)

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
        
### Atributes

gop : int gap_open_penalty; gep : int gap_extension_penalty; sequence_1 str, the content of sequence A; sequence_2 str, the content of sequence B; ID_1/ID_2, the IDs of sequence A/B; score_matrix : pandas dataframe, the reference score matrix for amino acid comparision; align_mtx, pandas dataframe, the computed dataframe for score calculation and sequence alignment; align_score, int, the computed number representing the similarity of two sequences; align_optimal_track, list, a list of location (representing by [index, columns] used for computing the alignment sequences; align_seq1/align_seq2, str, the aligned sequences of sequence A/B; trace_mtx, pandas dataframe, it is a matrix representing the "flow" during the computation. When during the sequence alignment, once the start point was found, I can just trace back to find the whole sequence.

### Methods

* Initialization

        Align = SmithWaterman('../sequences/prot-0091.fa', '../sequences/prot-0093.fa','../scoring_matrices/BLOSUM62.mat',2,1)
    
* Calculating the alignment score

        Align.compute_alignment_score()
       
* Return the alignment sequence

        Align.get_alignment_sequence()
        
* Representing the alignment score
        
        Align.get_alignment_score()
        
* alignment()

initialize align_mtx

* align()

Compute the aligned sequences




