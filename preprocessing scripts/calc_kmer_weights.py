#generate normalised k-mer weights (2:7)
#all possible k-mer combinations are generated
#followed by mapping them onto the k-mers of individual sequence files
#a tuple for every k-mer:transcript pair is created
#computationally extremely intensive

import os
from os.path import join, abspath, isfile
import pandas as pd
import numpy as np
from itertools import product, cycle
from multiprocess import Pool

def readf(k):
	flanks_down = pd.read_csv("fs1/w_fickett/down_fickett")#only fs1, as kmer tuples will be reused
#for the sequences are the same across fs1, fs2, fs3
	print (flanks_down.head())
	flanks_down.drop(["Unnamed: 0"], 1, inplace=True)
	genkmer(flanks_down, k)
	#fetchfiles(flanks_down,k)

def genkmer(flanks_down,k):
	kmer_insert=pd.DataFrame(columns=flanks_down.index)
	nt = ['A','T','G','C']
	#k = 7
	kmer3 = [''.join(p) for p in product(nt, repeat=k)]
	kmer_insert.insert(loc=0, column="kmer", value="NNN")
	for i in range(len(kmer3)):
		kmer_insert.loc[i, "kmer"] = kmer3[i]
	#kmer_insert.to_csv(join("tuples", "k7"))
	#print (kmer_insert.head())
	fetchfiles(flanks_down, kmer_insert,k)

def fetchfiles(flanks_down, kmer_insert, k):
	kmer_insert=kmer_insert.set_index('kmer')
	print (kmer_insert.head())
	#k=7
	cols = flanks_down['gene_id']
	runs = True
	col_cycle = cycle(cols)
	col_val = next(col_cycle)
	startc = col_val
	col_insert = 0
	col_present = []
#while runs:
	for col_val in cols:
    #print(col_val)
		if col_insert<464:
			if isfile(abspath(os.path.join("sno", "kmer"+str(k), col_val+".count"))):
				indi_t = pd.read_csv(abspath(join("sno", "kmer"+str(k), col_val+".count")), 
                            header=None, names=['kmer', 'freq'], delim_whitespace=True)
				indi_t=indi_t.set_index('kmer')
				for index in indi_t.index:
					for indekx in kmer_insert.index:
						if index==indekx:
							kmer_insert.at[indekx, col_insert]=indi_t.loc[index, "freq"] 
				col_insert += 1
            #col_present.append(col_val)        
		if col_insert>463 and col_insert<928:            
			if isfile(abspath(os.path.join("mir", "kmer"+str(k), col_val+".count"))):
				indi_t = pd.read_csv(abspath(join("mir", "kmer"+str(k), col_val+".count")), 
                            header=None, names=['kmer', 'freq'], delim_whitespace=True)
				indi_t=indi_t.set_index('kmer')            
				for index in indi_t.index:
					for indekx in kmer_insert.index:
						if index==indekx:
							kmer_insert.at[indekx, col_insert]=indi_t.loc[index, "freq"] 
				col_insert += 1
            #col_present.append(col_val)        
		if col_insert>927:                      
			if isfile(abspath(os.path.join("lnc", "kmer"+str(k), col_val+".count"))):
				indi_t = pd.read_csv(abspath(join("lnc", "kmer"+str(k), col_val+".count")), 
                            header=None, names=['kmer', 'freq'], delim_whitespace=True)
				indi_t=indi_t.set_index('kmer')            
				for index in indi_t.index:
					for indekx in kmer_insert.index:
						if index==indekx:
							kmer_insert.at[indekx, col_insert]=indi_t.loc[index, "freq"] 
				col_insert += 1
            #col_present.append(col_val)            
    #col_val = next(col_cycle)    
	print ("stop")
	print (kmer_insert.head())
	maketuples(kmer_insert, k)

def maketuples(kmer_insert, k):
	#k=7
	path="fs1/tuples"
	if not os.path.isdir(join(path, "kmer"+str(k))):
		os.makedirs(join(path, "kmer"+str(k)))
	kmer_insert.to_csv(abspath(join(path, "kmer"+str(k), "only_freq")))
	print ("saving file... done")

#calculate weight
#get every freq from every transcript, get the sum of freqs for that transcript
#divide freq/sum(freq)
	kmer_trans_w1 = pd.DataFrame(index=kmer_insert.index, columns=kmer_insert.columns)
	print ("Weights for kmer%d"%(k))
	for t in kmer_insert.columns:
		for index in kmer_insert.index:
			kmer_trans_w1.at[index, t]=kmer_insert[t].sum()
#kmer_trans_w=kmer_trans_w1.copy()
	for t in kmer_insert.columns:
		kmer_trans_w1[t]=kmer_insert[t]/kmer_trans_w1[t]

	kmer_trans_w1.to_csv(abspath(join(path, "kmer"+str(k), "only_weight")))
	print("saving file... done")

#assign key to every (kmer, trans) pair
#if present 1, otherwise 0/NaN?
	kmer_trans_key = pd.DataFrame(index=kmer_insert.index, columns=kmer_insert.columns)
	print ("Keys for kmer%d"%(k))
	for t in kmer_insert.columns:
		kmer_trans_key.at[kmer_insert[t].notnull(), t]=1
#kmer_trans_key

	kmer_trans_key.to_csv(abspath(join(path,  "kmer"+str(k), "only_key")))
	print ("saving file... done")

#join all the three dataframes
#create a tuple for each (kmer, trans) pair
#(key, weight, freq)
	kmer_trans_tup = pd.DataFrame(index=kmer_insert.index, columns=kmer_insert.columns)
	print ("Tuples for kmer%d"%(k))
	for t in kmer_insert.columns:
		for index in kmer_insert.index:
			tup = str(kmer_trans_key.loc[index, t]), str(kmer_trans_w1.loc[index, t]), str(kmer_insert.loc[index, t])        
			kmer_trans_tup.at[index, t] = tup
#kmer_trans_tup

	kmer_trans_tup.to_csv(abspath(join(path, "kmer"+str(k), "tuple_train")))
	print ("Saving master file... done")

def savetuples():
	path="fs1/tuples"
#	if os.path.isdir(join(path, "kmer"+str(k))):
#		os.makedirs(join(path, "kmer"+str(k)))
	for k in range(2,8):
		kmer = "kmer"+str(k)
		print (kmer)
		kmer_trans_tup = pd.read_csv(abspath(join(path, kmer, "tuple_train")), index_col="kmer", low_memory=False)
		tupt=kmer_trans_tup.T
		for ci in tupt.columns:    
			tupt.loc[(tupt[ci]=="('nan', 'nan', 'nan')"), ci] = ('0')
            #for i in range(len(ci)):
		if k==2:
			tupdf = tupt
		else:
			tupdf = pd.concat([tupdf, tupt], axis=1, copy=False)
    
	print (tupdf.info())
	tupdf.to_hdf(abspath(join(path, "tuple_cols_27")), key="k37", mode="w")

if __name__=="__main__":
	os.chdir("../training_1")
#	with Pool() as p:
	with Pool(16) as p: 
       #print(p.map(f, [1, 2, 3]))
#		x=p.apply(savetuples)
#		x.get()
		p.map(readf, range(2,8))
		#savetuples()
		y=p.apply(savetuples)
		#y.get()

