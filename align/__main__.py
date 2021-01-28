# All the codes for part 2


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
from algs import  SmithWaterman, score_matrix, NeedlemanWunsch
import seaborn as sns
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


def compute_SmithWaterman_Score(list1, list2, truth,index,mtx_dir,gop,gep,filename):
    print(index)
    Align = SmithWaterman(list1[index], list2[index],mtx_dir,gop,gep)
    #Align_scores.append(Align.compute_alignment_score())
    #Indicators.append(ground_truth[index])
    score = Align.compute_alignment_score()
    with open(filename+".txt", "a") as text_file:
        text_file.write(str(index)+'\t'+str(truth[index])+'\t'+str(score)+'\n')
        
def compute_NeedlemanWunsch_Score(list1, list2, truth,index,mtx_dir,gop,gep,filename):
    print(index)
    Align = SmithWaterman(list1[index], list2[index],mtx_dir,gop,gep)
    #Align_scores.append(Align.compute_alignment_score())
    #Indicators.append(ground_truth[index])
    score = Align.compute_alignment_score()
    with open(filename+".txt", "a") as text_file:
        text_file.write(str(index)+'\t'+str(truth[index])+'\t'+str(score)+'\n')



# Question 1:
'''
With the BLOSUM50 matrix and a gap opening cost of 11 and a gap extension
cost of 3, locally align these sequences and visualize the distribution of
alignment scores. How would you describe this distribution?
'''
# read the direction of all sequence pairs and record the ground truth. Pospairs as 1 and Negpairs as 0


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
        
for i in range(0,len(sequence_dirA)):
    compute_SmithWaterman_Score(sequence_dirA,sequence_dirB, ground_truth, i, '../scoring_matrices/BLOSUM50.mat',11,3,'BLOSUM50_11_3.txt')
    
table = pd.read_csv('BLOSUM50_11_3.txt.txt', sep ='\t', header = None, index_col= None)

sns.set(style ='whitegrid')

ax = sns.boxplot(table[1],table[2], showfliers=False )
#ax = sns.swarmplot(table[1], table[2], color = '.25', showfliers =False)


# Question 2:
'''
Generate a confusion matrix indicating the frequency of false positives, false
negatives, true positives, and true negatives when using the average alignment
score as a threshold. What is the threshold value, and how does the confusion
matrix suggest this algorithm performed?
'''

average = stat.mean(table[2])
Predict = [int(x) for x in [x >average for x in table[2]]]
cm = skl_metrics.confusion_matrix( ground_truth,Predict)
cm_display = skl_metrics.ConfusionMatrixDisplay(cm)
cm_display.plot()

# Questions 3:
'''
Create a ROC plot which shows the fraction of true positives on the Y axis and
the fraction of false positives on the X axis. Please take care to make your ROC
plots square, with both X and Y axes limited to the range [0:1].
'''

Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()

# Questions 4:
'''
Determine the area under the ROC curve (AUROC). What does this value
indicate about the performance of this algorithm? Can you confidently assess
the performance of this algorithm with this value alone? Why or why not?
'''

skl_metrics.roc_auc_score(ground_truth, Predict_ratio)

# Questions 5
'''
Once again, using local alignment, try a range of gap opening (1-20) and gap
extension (1-5) costs with the BLOSUM62 matrix. Using the AUROC of each
approach, determine which gap penalty performs the “best”. What does this pair
of values suggest about the evolution of these sequences and the likelihood of
insertions / deletions?
'''
# Don't run this ... it can take a while ... I did this on the server. Please read mult_task.py if you need more information,
# Thanks
for i in range(1,21):
    for j in range(1,6):
        label = str(i) + '+' + str(j)
        print('nohup python mult_task.py --gep '+str(i)+ ' --gop '+ str(j)+ ' >log/'+label+'.log &')
        os.system('nohup python mult_task.py --gep '+str(i)+ ' --gop '+ str(j)+ ' >log/'+label+'.log &')
        
roc_scores = []
IDs = []
GOP = range(1,21)
GEP = range(1,6)
for i in GOP:
    for j in GEP:
        ID = str(i)+'-'+str(j)
        IDs.append(ID)
        table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
        Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
        fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
        #roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
        roc_scores.append(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))
        
ID = IDs[int(np.where(roc_scores == max(roc_scores))[0])]
print(ID)
table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))

# Questions 6,7:
'''
Using the BLOSUM50, BLOSUM62, PAM100 and PAM250 scoring matrices,
evaluate the performance of the global alignment algorithm using the selected
pair of best performing gap penalties. Evaluate your optional additional
extension algorithm as well, using the parameters it requires.

For each algorithm, generate a ROC plot demonstrating performance using
each of the 4 matrices, given the fixed gap costs. Of all these algorithms, which
performs the best as measured by AUROC?
'''
# Qurstions 8:
'''
Comment qualitatively on the best algorithm. What does the best performing
algorithm indicate about the origin of these sequences?
'''

for i in range(0,100):
    compute_NeedlemanWunsch_Score(sequence_dirA,sequence_dirB, ground_truth, i ,'../scoring_matrices/BLOSUM50.mat',5,4,'N_BLO50_5_4')
    compute_NeedlemanWunsch_Score(sequence_dirA,sequence_dirB, ground_truth, i ,'../scoring_matrices/BLOSUM62.mat',5,4,'N_BLO62_5_4')
    compute_NeedlemanWunsch_Score(sequence_dirA,sequence_dirB, ground_truth, i ,'../scoring_matrices/PAM100.mat',5,4,'N_PAM100_5_4')
    compute_NeedlemanWunsch_Score(sequence_dirA,sequence_dirB, ground_truth, i ,'../scoring_matrices/PAM250.mat',5,4,'N_PAM250_5_4')                              

#N_PAM250         
table = pd.read_csv('N_PAM250_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))

#N_PAM100         
table = pd.read_csv('N_PAM100_5_4.txt', sep ='\t', header = None, index_col= None)


Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))

#N_BLOSUM62         
table = pd.read_csv('N_BLO62_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))

#N_BLOSUM50         
table = pd.read_csv('N_BLO50_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))



#BLOSUM50         
table = pd.read_csv('BLOSUM50_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))


#BLOSUM62        
table = pd.read_csv('5-4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))


#PAM250        
table = pd.read_csv('PAM250_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))


#PAM100        
table = pd.read_csv('PAM100_5_4.txt', sep ='\t', header = None, index_col= None)

#table = pd.read_csv(ID+'.txt', sep ='\t', header = None, index_col= None)
Predict_ratio = [(x - min(table[2])/(max(table[2]) - min(table[2]))) for x in table[2]]
fpr, tpr, _ = skl_metrics.roc_curve(ground_truth, Predict_ratio, pos_label=1)
roc_display = skl_metrics.RocCurveDisplay(fpr=fpr, tpr=tpr).plot()
print(skl_metrics.roc_auc_score(ground_truth, Predict_ratio))



            








    
        
