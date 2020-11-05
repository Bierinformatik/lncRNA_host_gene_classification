#calculate fickett score for all individual sequences
from fickett import *
from itertools import cycle
from os.path import abspath, join
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import os

#set the path, since all sequences are in individual files
path = abspath("../training_1/sno")
files = [f for f in os.listdir(path) if os.path.isfile(join(path, f))]
for x in files:
	with open (abspath(join(path, x)), "r") as f:
		flines = f.readlines()
		flen = len(flines)
		fcycle = cycle(flines)
		read = next(fcycle)
		#print (read[1:].rstrip())
		start = read
		runs = True
		fa_name = [] 
		fa_real = ""
		fa_c=0
		line = 1
		while runs:				
			if read.startswith(">"):
			#to extract the name of the entry and avoid having duplicates
			#since multiple extracted sequences can come from a single transcript/read
			
				#fa_id = read.split()[0]
				#x= (read[1:].rstrip())
				#if not read.find('|')==-1:
				#name_ex = fa.id[fa.id.index('|')+1:fa.id.index(':')]
				#	name_ex = read[:read.index('|')]
				#else:
				#	name_ex = read[:read.index(':')]
				name_ex = read[1:].rstrip()
				'''if not name_ex in fa_name:
					fa_name.append(name_ex)
					fa_c=0			
					fa_real=name_ex+".0"
				else:
					c = Counter(fa_name)
					fa_c = c[name_ex]+1
					fa_name.append(name_ex)
					fa_real=name_ex+"."+str(fa_c)'''
				read = next(fcycle).rstrip().upper()			
				fs = fickett.fickett_value(read)
				print (f'{name_ex}\t{fs}')
			read = next(fcycle)
			line+=2
			if line>flen:
				runs = False

