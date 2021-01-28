import multiprocessing
import fastaparser
import glob
import os
import sys
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat
import sklearn.metrics as skl_metrics
import fastaparser
from algs import  SmithWaterman
import Compute_Scores
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gep',type=int, default= 0, help='Gap extension penalty')
parser.add_argument('--gop', type = int,default =0,help='Gap open penalty')
args = parser.parse_args()


ground_truth  = []
sequence_dirA = []
sequence_dirB = []
with open('../scoring_matrices/Negpairs.txt' , 'r') as f:
    
    
    for line in f:
        line = line.replace('\n','')
        seq1 = '../'+line.split(' ')[0]
        seq2 = '../'+line.split(' ')[1]
        sequence_dirA.append(seq1)
        sequence_dirB.append(seq2)
        ground_truth.append(0)

with open('../scoring_matrices/Pospairs.txt' , 'r') as f:
    
    
    for line in f:
     
        line = line.replace('\n','')
    
        seq1 = '../'+line.split(' ')[0]
        seq2 = '../'+line.split(' ')[1]
     
        sequence_dirA.append(seq1)
        sequence_dirB.append(seq2)
        ground_truth.append(1)
Align_scores = []
Indicators = []


label = str(args.gop) +'-'+str(args.gep)
for i in range(0,len(sequence_dirA)):
    print(i)
    
    Compute_Scores.compute_SmithWaterman_Score(sequence_dirA,sequence_dirB, 
                    ground_truth,i,'../scoring_matrices/BLOSUM62.mat',args.gop,args.gep, label)
print(i)    
print('what?')