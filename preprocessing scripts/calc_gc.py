#calculates GC content 

import os
from decimal import Decimal
from Bio import SeqIO
from collections import Counter

os.chdir("../training_1/sno")
fa_name=[]#""
fa_real=""
fa_c=0
files = [x for x in os.listdir() if os.path.isfile(x)]
#f="Merged_Filtered_Closest_Exon_sno.fa"
for f in files:
	for fa in SeqIO.parse(open (f, "r"), "fasta"):
	#print (fa.id)	
	#extracts sequences, avoids duplicates in names
	#depending upon the gene class (lnc are named differently than snhg and mirhg) 
	#the following section may be uncommented
	
		'''line = fa.id.split(":")
		name_ex = line[0][line[0].rindex("|")+1:]
		#name_ex = line[0]#for lnc
		if not name_ex in fa_name:
			fa_name.append(name_ex)
			fa_c=0
			fa_real=name_ex+".0"		
		else:
			c = Counter(fa_name)
			fa_c = c[name_ex]+1
			fa_name.append(name_ex)		
			fa_real=name_ex+"."+str(fa_c)'''
	#print (fa_real)
		fa_real = fa.id
		g, c = 0, 0
		seq = str(fa.seq).upper()
		for s in seq:
			if s=="G":
				g+=1
			elif s=="C":
				c+=1
		gc = (g+c)/len(seq)
		with open(os.path.abspath(os.path.join("..", "gc.lnc")), "a") as w:
			w.write(f'{fa_real}\t{Decimal(gc).quantize(Decimal("1.00"))}\n')#\t{len(seq)}\n')
print ("fin")
