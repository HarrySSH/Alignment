import os
import sys
import glob
import pandas as pd
import numpy as np
import fastaparser


def score_matrix(dir):
	'''
	With the direction of the score matrix files, this function build a pandas dataframe and return it.
	'''
	with open(dir , 'r') as f:
		p = 0
		for line in f:
			if line.startswith('#'):
				continue
			else:
		
				string = line.strip('\n')
				string = string.strip(' ') 
				string = string.split(' ')
				string = [x for x in string if len(x)>=1]
				
				if p == 0:
					df = pd.DataFrame(columns=string)
					index = string
					p = p+1
				else:
					string = [int(x) for x in string]

					df_ = pd.DataFrame(string, index = index).T
					df = df.append(df_, ignore_index=True)
					#p = p+1
					#df = pd.DataFrame(columns=string, )
					#df = df.append(string, )
		
		df.index = index
		return df

def fasta_reader(dir):
	'''
	With the direction of fasta files, this function return a sequence class. 
	For using it:
		seq = fasta_reader(dir)
		The_actual_sequence = seq.sequence_as_string()
	'''
	with open(dir) as fasta_file:
		parser = fastaparser.Reader(fasta_file)
		
		for seq in parser:
			return seq
class PairwiseAligner:
	'''
	The basic alignment class.
	The parent class of  NeedlemanWunsch (Global alignment) and SmithWaterman (Local Alignment) class.

	If you derectly use Pairwise class, the default method is global alignment.
	'''
	def __init__(self, seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty):
		'''
		Parameters: 
		
		seq1 : str, default = NA, The dir of sequence A

		seq2 : str, default = NA, The dir of sequence B

		mtx_dir : str, default = NA, The dir of score matrix

		gap_open_penalty : int, default = NA, The value of gap open penalty

		gap_extension_penalty, int, default = NA. The value of gap extension penalty
		'''
		ID_1 = seq1.split('/')[-1].split('.')[0]    # Initialize the ID
		ID_2 = seq2.split('/')[-1].split('.')[0]
		self.ID_1 = ID_1
		self.ID_2 = ID_2


		seq1 = fasta_reader(seq1)    # Initialize the sequence
		seq2 = fasta_reader(seq2)
		self.sequence_1 = seq1.sequence_as_string()
		self.sequence_2 = seq2.sequence_as_string()
		self.sequence_1 = "*" + self.sequence_1	
		self.sequence_2 = "*" + self.sequence_2


		#gap__open_penalty
		self.gop = gap_open_penalty
		#gap_extension_panalty
		self.gep = gap_extension_penalty
		
		# Initialize the score matrix
		self.score_matrix = score_matrix(mtx_dir)

		# Defualt is zero
		self.align_score = 0

		# The track during alignment score computing, for recording the movement of location. 
		# In this way, it would be a lot easier when return the alignment sequences
		self.align_optimal_track = []
		

		# The aligned sequences, default is ''.
		self.align_seq1 =''
		self.align_seq2 =''

		# The alignment matrix, first intialize by a matrix with the dimension equal to te length of seq1 and seq2.
		# Later this matrix will be filled according to the score matrix and the standard of walk (different between global and local)
		
		# The trace matrix, this matrix representing the "flow" during the alignment matrix initialization. 
		# First build a pandas dataframe, each elment is intialize by a two element list. 
		# eg for a element in the location [4,6] on alignment matrix,  its corresponding [[4,5], [[5,7],[5,6]]]
		# The first element : [4,5] representing the previous location that pointing at this location. 
		# Initialize by int(zero), demostrating no previous location can point to this one. 
		# The second element : [[5,7],[5,6]] is a list  that point at the next element/elements since a element on score matrix 
		# can point to more than one element.
		# Initilize by [0], means that there is no element that this location point to.
		self.trace_mtx = pd.DataFrame(np.zeros((len(self.sequence_1), len(self.sequence_2))), index =list(self.sequence_1), columns = list(self.sequence_2))
		self.align_mtx = pd.DataFrame(np.zeros((len(self.sequence_1), len(self.sequence_2))), index =list(self.sequence_1), columns = list(self.sequence_2))
		self.trace_mtx = self.trace_mtx.astype('object')
		h = self.align_mtx.shape[0]
		w = self.align_mtx.shape[1]	
		# initialize
		for i in range(0,h):
			for j in range(0,w):
				self.trace_mtx.iloc[i,j] = [0,[0]]
		
    # Since we are doing afine gap penalty, then everytime when we want to have a gap we need to justify whether we
    # should give a open penalty or externsion penalty. To decide on this, my method is to judge whether the previous point
    # already began the gap or not. For this penalty funtion, [i,j] is the location of either the upper location or the left location.
    # We juade whether this previous location is 1) linked to its left or upper 2) linked to the upper left (diagnoal) or don't have links. if 1), give the 
    # following one extension penalty. else, (Just started the sequence or the gap did not occur) give the open penalty
	def penalty(self,i,j):	
		
		if (self.trace_mtx.iloc[i,j][0] == 0) :
			return self.gop
		elif ([i,j] in self.trace_mtx.iloc[i-1,j][1]) or ([i,j] in self.trace_mtx.iloc[i,j-1][1]) :
			return self.gep

		elif ([i,j] in self.trace_mtx.iloc[i-1,j-1][1]):
			return self.gop
		else:
			assert 2==1, 'SOmething went wrong!'


    # Alignment is the function to fill the align_mtx, the default id global alignment. 

    # What is different is, when I computng the matrix, I also track the links or tragetories and store it in the track_mtx
	def alignment(self):
		h = self.align_mtx.shape[0]
		w = self.align_mtx.shape[1]	
		# compute the matrix
		self.align_mtx.iloc[0,0] = self.score_matrix[self.sequence_1[0]][self.sequence_2[0]]
		
		for i in range(1,h):
			self.align_mtx.iloc[i,0]= self.score_matrix[self.sequence_1[i]][self.sequence_2[0]] + self.score_matrix[self.sequence_1[0]][self.sequence_2[0]]
			
		for j in range(1,w):
			self.align_mtx.iloc[0,j]= self.score_matrix[self.sequence_1[0]][self.sequence_2[j]] + self.score_matrix[self.sequence_1[0]][self.sequence_2[0]]
			
		
		for i in range(1,h):
			for j in range(1,w):

					self.align_mtx.iloc[i,j] = max(
						
							self.align_mtx.iloc[i-1,j] - self.penalty(i-1,j),
							self.align_mtx.iloc[i-1,j-1] + self.score_matrix[self.sequence_1[i]][self.sequence_2[j]],
							self.align_mtx.iloc[i,j-1] - self.penalty(i,j-1)
							
						
					)
					# After decide which direction to go. For example, if in align_matrix, [4,5] ---> [5,6] then in the 
                    # track_matrix, 
                    # For the previous location [4,5] : self.trace_mtx.iloc[4,5][1].append([5,6]) , add one more link (Also remove 0 if 0 is the only element in)
                    # The the later location [5,6] : self.trace_mtx.iloc[5,6][0] = [4,5]
                    # In this way, the track chain can be built, which would be super helpful for further alignment.
					if self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i-1,j] - self.penalty(i-1,j):
						if self.trace_mtx.iloc[i-1,j][1] == [0]:
							self.trace_mtx.iloc[i-1,j][1].append([i,j])
							self.trace_mtx.iloc[i-1,j][1] = self.trace_mtx.iloc[i-1,j][1][self.trace_mtx.iloc[i-1,j][1]!= 0]
							self.trace_mtx.iloc[i-1,j][1] = [self.trace_mtx.iloc[i-1,j][1]]
						else:
							assert 	self.trace_mtx.iloc[i-1,j][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i-1,j][1].append([i,j])


						self.trace_mtx.iloc[i,j][0] = [i-1,j]
					elif self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i-1,j-1] + self.score_matrix[self.sequence_1[i]][self.sequence_2[j]]:
						if self.trace_mtx.iloc[i-1,j-1][1] == [0]:
							self.trace_mtx.iloc[i-1,j-1][1].append([i,j])
							self.trace_mtx.iloc[i-1,j-1][1] = self.trace_mtx.iloc[i-1,j-1][1][self.trace_mtx.iloc[i-1,j-1][1] != 0]
							self.trace_mtx.iloc[i-1,j-1][1] = [self.trace_mtx.iloc[i-1,j-1][1]]
						else:
							assert 	self.trace_mtx.iloc[i-1,j-1][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i-1,j-1][1].append([i,j])
						
						self.trace_mtx.iloc[i,j][0]= [i-1, j-1]
					elif self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i,j-1] - self.penalty(i,j-1):
						if self.trace_mtx.iloc[i,j-1][1] == [0]:
							self.trace_mtx.iloc[i,j-1][1].append([i,j])
							self.trace_mtx.iloc[i,j-1][1] = self.trace_mtx.iloc[i,j-1][1][self.trace_mtx.iloc[i,j-1][1] != 0]
							self.trace_mtx.iloc[i,j-1][1] = [self.trace_mtx.iloc[i,j-1][1]]
						else:
							assert 	self.trace_mtx.iloc[i,j-1][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i,j-1][1].append([i,j])
						self.trace_mtx.iloc[i,j][0]= [i, j-1]

					else:
						assert 1 ==2 , "Something went wrong"

	'''
	Computing scores became eaiser with the help of track matrix, for global alignment, just starting from the right down conor, 
	and always return to the previous location. And at the same time, the location was added into align_optimal_track for the further computing 
	of the align sequence.
	'''
	def compute_alignment_score(self):
		self.alignment() # Call alignment to perform the calculation



		# start from the biggest value
		
		h = int(self.align_mtx.shape[0]) -1
		w = int(self.align_mtx.shape[1]) -1
		
		self.align_score = self.align_mtx.iloc[h,w]

		self.align_optimal_track.append((h,w))
		point = (h,w)
		h = point[0]
		w = point[1]			

		while (h!= 0 and w!=0):

				
			point = self.trace_mtx.iloc[h,w][0]
			h = point[0]
			w = point[1]
			self.align_score = self.align_score + self.align_mtx.iloc[h,w]
			self.align_optimal_track.append((h,w))
				
		self.align()
		return self.align_score
    

    # In the end of function compute_alignment_score, this function was performed to put the resulting aligned sequences into self.align_sequence
	def align(self):
		p = 0
		self.align_seq1 ==''
		self.align_seq1 ==''
		#assert self.align_seq1 !='', 'Something went wrong!'
		for loc in self.align_optimal_track:
            
			#assert self.align_seq1 =='', 'Something went wrong!'
			if p == 0:
				h = loc[0]
				w = loc[1]
				self.align_seq1 = self.sequence_1[h] + self.align_seq1
				self.align_seq2 = self.sequence_2[w] + self.align_seq2
				p = p +1
			else:
				_h = loc[0]
				_w = loc[1]
      
				if h == _h +1:
					self.align_seq1 = self.sequence_1[_h] + self.align_seq1
				elif h == _h:
					self.align_seq1 = '_' + self.align_seq1
				else:
					assert 1 ==2, 'Something went wrong'
				if w == _w +1:
					self.align_seq2 = self.sequence_2[_w] + self.align_seq2
				elif w == _w:
					self.align_seq2 = '_' + self.align_seq2
				else:
					assert 1 ==2, 'Something went wrong'
				h = _h
				w = _w
    # We can call this to function to return the sequcens or scores, but do make sure that you run compute_alignment_score function first!
	def get_alignment_score(self):
		print('----------------------')
		##assert self.align_score == 100000000, "The alignment did not happen!"
		print('The alignment score between '+str(self.ID_1) + ' and ' + str(self.ID_2) +' is: ' + str(self.align_score))
		return  self.align_score



	def get_alignment_sequence(self):
		
		print('The alignment results:')
		print('----------------------')
		print(self.ID_1 +':')
		print(self.align_seq1)
		print(self.ID_2 +':')
		print(self.align_seq2)
		return self.align_seq1, self.align_seq2



	
		



class SmithWaterman(PairwiseAligner):
	'''
	Similar with global alignment but this time the alignemtn starts from the highest value and stop when there is a zero value
	'''
	
	
	
	def __init__(self, seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty):
		# Initialize attributes of the parent class
		super().__init__(seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty)
	'''
	alignment is a little different here because whenever there is a zero, the function stop

	'''
	def alignment(self):
		h = self.align_mtx.shape[0]
		w = self.align_mtx.shape[1]	
		# compute the matrix

		self.align_mtx.iloc[0,0] = max(self.score_matrix[self.sequence_1[0]][self.sequence_2[0]],0)
		
		for i in range(1,h):
			self.align_mtx.iloc[i,0]= max(self.score_matrix[self.sequence_1[i]][self.sequence_2[0]] + self.score_matrix[self.sequence_1[0]][self.sequence_2[0]],0)
			
		for j in range(1,w):
			self.align_mtx.iloc[0,j]= max(self.score_matrix[self.sequence_1[0]][self.sequence_2[j]] + self.score_matrix[self.sequence_1[0]][self.sequence_2[0]],0)
			
		
		for i in range(1,h):
			for j in range(1,w):
				    
					self.align_mtx.iloc[i,j] = max(
						
							self.align_mtx.iloc[i-1,j] - self.penalty(i-1,j),
							self.align_mtx.iloc[i-1,j-1] + self.score_matrix[self.sequence_1[i]][self.sequence_2[j]],
							self.align_mtx.iloc[i,j-1] - self.penalty(i,j-1),
							0
						
					)
					

					#assert self.trace_mtx.iloc[i-1,j][1] == 0, 'haha I got you bitch'
					#assert self.trace_mtx.iloc[i,j-1][1] == 0, 'haha I got you bitch'
					if self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i-1,j] - self.penalty(i-1,j):
						if self.trace_mtx.iloc[i-1,j][1] == [0]:
							self.trace_mtx.iloc[i-1,j][1].append([i,j])
							self.trace_mtx.iloc[i-1,j][1] = self.trace_mtx.iloc[i-1,j][1][self.trace_mtx.iloc[i-1,j][1]!= 0]
							self.trace_mtx.iloc[i-1,j][1] = [self.trace_mtx.iloc[i-1,j][1]]
						else:
							assert 	self.trace_mtx.iloc[i-1,j][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i-1,j][1].append([i,j])


						self.trace_mtx.iloc[i,j][0] = [i-1,j]
					elif self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i-1,j-1] + self.score_matrix[self.sequence_1[i]][self.sequence_2[j]]:
						if self.trace_mtx.iloc[i-1,j-1][1] == [0]:
							self.trace_mtx.iloc[i-1,j-1][1].append([i,j])
							self.trace_mtx.iloc[i-1,j-1][1] = self.trace_mtx.iloc[i-1,j-1][1][self.trace_mtx.iloc[i-1,j-1][1] != 0]
							self.trace_mtx.iloc[i-1,j-1][1] = [self.trace_mtx.iloc[i-1,j-1][1]]
						else:
							assert 	self.trace_mtx.iloc[i-1,j-1][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i-1,j-1][1].append([i,j])
						
						self.trace_mtx.iloc[i,j][0]= [i-1, j-1]
					elif self.align_mtx.iloc[i,j] == self.align_mtx.iloc[i,j-1] - self.penalty(i,j-1):
						if self.trace_mtx.iloc[i,j-1][1] == [0]:
							self.trace_mtx.iloc[i,j-1][1].append([i,j])
							self.trace_mtx.iloc[i,j-1][1] = self.trace_mtx.iloc[i,j-1][1][self.trace_mtx.iloc[i,j-1][1] != 0]
							self.trace_mtx.iloc[i,j-1][1] = [self.trace_mtx.iloc[i,j-1][1]]
						else:
							assert 	self.trace_mtx.iloc[i,j-1][1][0]!= 0, "Something went wrong!"
							self.trace_mtx.iloc[i,j-1][1].append([i,j])
						self.trace_mtx.iloc[i,j][0]= [i, j-1]
					elif self.align_mtx.iloc[i,j] == 0:
						pass
					else:
						assert 1 ==2 , "Something went wrong"
					
						



	def compute_alignment_score(self):

		self.alignment()


		# start from the biggest value, the only different with the global alignment here.

		start = np.max(self.align_mtx.to_numpy())
		loc = np.where(self.align_mtx == start)

		if (len(loc) ==2 and (len(loc[0]) ==1)) :
			self.align_score = start

			h = int(loc[0])
			w = int(loc[1])
			self.align_optimal_track.append((h,w))
			point = (h,w)
			h = point[0]
			w = point[1]			

			while ((self.align_mtx.iloc[h,w]!= 0 and w != 0)and h != 0):

				
				point = self.trace_mtx.iloc[h,w][0]

				h = point[0]
				w = point[1]

				self.align_score = self.align_score + self.align_mtx.iloc[h,w]
				self.align_optimal_track.append((h,w))
			self.align()
				
		else:
			score_list = []
			track_list = []
			num = len(loc[0])
			for chop in range(0,num):
			
				self.align_score = start
				h = int(loc[0][chop])
				w = int(loc[1][chop])
				self.align_optimal_track.append((h,w))
				point = (h,w)
				
				h = point[0]
				w = point[1]
							
				while ((self.align_mtx.iloc[h,w]!= 0 and w != 0)and h != 0):


					point = self.trace_mtx.iloc[h,w][0]
					h = point[0]
					w = point[1]
					self.align_score = self.align_score + self.align_mtx.iloc[h,w]
					self.align_optimal_track.append((h,w))
                    
				score_list.append(self.align_score)
				track_list.append(self.align_optimal_track)
				self.align_optimal_track=[]
			index = score_list.index(max(score_list))
			self.align_optimal_track = track_list[index]
			self.align_score = score_list[index]
			self.align()
		return self.align_score

							







class NeedlemanWunsch(PairwiseAligner):
	def __init__(self, seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty):
		# Initialize attributes of the parent class
		super().__init__(seq1, seq2, mtx_dir, gap_open_penalty, gap_extension_penalty)


    






