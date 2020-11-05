import numpy as np
import pandas as pd
import os
import argparse
from training_1 import *
from training_2 import *
from training_3 import *
from training_4 import *
import rf

def choose_training(*t):
	if args.t==1:
		print ("Dataset 1 selected")
		os.chdir("training_1")
		if args.fs==1:
			features_combined=training1_fs1(args.fi)
			return (features_combined)
		elif args.fs==2:
			features_combined=training1_fs2(args.fi)
			return (features_combined)
		elif args.fs==3:
			features_combined=training1_fs3(args.fi)
			return (features_combined)
	elif args.t==2:
		print ("Dataset 2 selected")
		os.chdir("training_2")
		if args.fs==1:
			features_combined=training2_fs1(args.fi)
			return (features_combined)
		elif args.fs==2:
			features_combined=training2_fs2(args.fi)
			return (features_combined)
		elif args.fs==3:
			features_combined=training2_fs3(args.fi)
			return (features_combined)
	elif args.t==3:
		print ("Dataset 3 selected")
		os.chdir("training_3")
		if args.fs==1:
			features_combined=training3_fs1(args.fi)
			return (features_combined)
		elif args.fs==2:
			features_combined=training3_fs2(args.fi)
			return (features_combined)
		elif args.fs==3:
			features_combined=training3_fs3(args.fi)
			return (features_combined)
	elif args.t==4:
		print ("Dataset 4 selected")
		os.chdir("training_4")
		if args.fs==1:
			features_combined=training4_fs1(args.fi)
			return (features_combined)
		elif args.fs==2:
			features_combined=training4_fs2(args.fi)
			return (features_combined)
		elif args.fs==3:
			features_combined=training4_fs3(args.fi)
			return (features_combined)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description="performs supervised training to classify the three classes: NoHG, SNHG, MIRHG")
	parser.add_argument("-t", "-training_set", help="specify data set to train upon",type=int,choices=[1,2,3,4],default=1)
	parser.add_argument("-fs","-feature_set", help="specify feature set: 1 for fs1 OR 2 for fs2 OR 3 for fs3", type=int,choices=[1,2,3],default=1)
	parser.add_argument("-fi","-fickett_score",help="specify inclusion of fickett score: 1 to enable, 0 to disable",type=int,choices=[0,1],default=0)
	args = parser.parse_args()
	#print (args)
	features_combined=choose_training(args.t,args.fs,args.fi)
	print (features_combined.info())
	rf.preprocess(features_combined,args.fs)
