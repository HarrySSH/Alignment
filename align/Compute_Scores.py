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



def compute_SmithWaterman_Score(list1, list2, truth,index,mtx_dir,gop,gep,filename):
    print(index)
    Align = SmithWaterman(list1[index], list2[index],mtx_dir,gop,gep)
    #Align_scores.append(Align.compute_alignment_score())
    #Indicators.append(ground_truth[index])
    score = Align.compute_alignment_score()
    with open(filename+".txt", "a") as text_file:
        text_file.write(str(index)+'\t'+str(truth[index])+'\t'+str(score)+'\n')