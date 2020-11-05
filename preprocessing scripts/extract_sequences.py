#extracting sequences from datasets

from os.path import abspath, join, isdir, isfile
from itertools import cycle
import os
from Bio import SeqIO
from collections import Counter

if not isdir("../training_1"):
        os.makedirs("../training_1")
path = "../training_1"
fa_name=[]#""
fa_real=""
fa_c=0
count=0#for lncs
for fa in SeqIO.parse(open(join("../Dataset1/", "sno_exons.fa"), "r"), "fasta"):
	#name_ex=fa.id.rstrip()
	#line = fa.id.split(":")
	#name_ex = line[0][line[0].rindex("|")+1:]
	name_ex = fa.id#for lnc
	if not name_ex in fa_name:
		fa_name.append(name_ex)
		fa_c=0
		fa_real=name_ex+".0"		
	else:
		c = Counter(fa_name)
		fa_c = c[name_ex]+1
		fa_name.append(name_ex)		
		fa_real=name_ex+"."+str(fa_c)
	if not isdir(join(path, "sno")):# mir/sno/lnc
		os.makedirs(join(path, "sno"))
	count+=1
	#if count<751:#only for lnc
	with open(join(path, "sno", fa_real), "w") as w:
		w.write(f'>{fa_real}\n{str(fa.seq).upper()}')

print(f'fin')
