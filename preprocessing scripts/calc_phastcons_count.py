#conservation scores

from os.path import abspath, join, isdir, isfile
from itertools import cycle
import os
from Bio import SeqIO
from collections import Counter
from statistics import mean

path = "../Dataset4"
fa_name=[]#""
fa_real=""
fa_c=0
for fa in SeqIO.parse(open(join(path, "Merged_PhastCons_exons_without_payload_two_four_windows.fa"), "r"), "fasta"):
	#name_ex=fa.id.rstrip()
	line = fa.id.split("_")
	#name_ex = line[1][:line[1].index("_")]#mir
	name_ex = line[3]#for lnc
	if not name_ex in fa_name:
		fa_name.append(name_ex)
		fa_c=0
		fa_real=name_ex+".0"		
	else:
		c = Counter(fa_name)
		fa_c = c[name_ex]+1
		fa_name.append(name_ex)		
		fa_real=name_ex+"."+str(fa_c)
	seq = str(fa.seq)
	cons = seq.rstrip().split(";")
	#print (cons[:-1])
#only for windows check if length>100	
#	if len(cons)>100:
	if not isdir(join(path, "training_4", "cons")):
		os.makedirs(join(path, "training_4", "cons"))
	with open(join(path, "training_4", "cons", fa_real+".cons"), "w") as w:						
		#cons = [float(x) for x in cons]
		icons = []
#only for windows cons[:100]
		for sc in cons[:-1]:
			if sc=="0":
				icons.append(0.0)
			else:
				icons.append(float(sc))
						#print(facons.id)	
		for i in range(len(icons)):
			w.write(f'{i+1}\t')		
			w.write(f'{icons[i]}\n')

print(f'fin')
