#bins for conservation scores and pairing probabilities

import numpy as np
from itertools import cycle
import os
from subprocess import call
from scipy import stats

os.chdir(os.path.join("..","training_3","lnc", "pp_60"))#or cons
print (os.getcwd())
#call("ls > reads", shell=True)
files = [x for x in os.listdir() if os.path.isfile(x)]
for pp in files:
#with open ("reads", "r") as f:
	if not os.path.isdir(os.path.join("pp_60_bins")):
		os.mkdir(os.path.join("pp_60_bins"))
	#while (pruns):
	with open(os.path.join("pp_60_bins", pp.rstrip()+".bins"), "w") as b:
		with open(pp.rstrip(), "r") as p:
			lines = p.readlines()
			prob_list = []
				#print (pp.rstrip())
			for x in lines:
				prob_list.append(x.split("\t")[1].rstrip())				
			prob_arr = np.asarray(prob_list).astype(np.float)
			#print (np.size(prob_arr))
				#prob_arr_sorted = np.sort(prob_arr)
				#print (prob_arr_sorted)
				#prob_arr_split = np.array_split(prob_arr_sorted, 20)
			edges = []
			for y in np.arange(0., 1.05, 0.05):
				edges.append(round(y,2))
				#print (edges)
				#print (prob_arr_split)
			bins = stats.binned_statistic(prob_arr, prob_arr, 'count', bins=edges)
				#print (bins)
				#print (len(bins[1]))
				#print (len(bins[0]))
			for each in range(len(bins[0])):
				if (each+1)<len(bins[0]):				
					b.write(f'{bins[0][each]}\t[{bins[1][each]} - {bins[1][each+1]})\n')
				else:
					b.write(f'{bins[0][each]}\t[{bins[1][each]} - {bins[1][each+1]}]\n')	

