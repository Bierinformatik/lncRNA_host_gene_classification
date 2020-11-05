import os
import pandas as pd

def training1_fs1(fickett):
	if fickett==1:
		flanks_down = pd.read_csv("fs1/w_fickett/down_fickett", index_col=0)
		print (f'Feature set 1 selected with fickett score:\n{flanks_down.info()}')
	elif fickett==0:
		flanks_down = pd.read_csv("fs1/wo_fickett/down_no_fickett", index_col=0)
		print (f'Feature set 1 selected without fickett score:')
	print (flanks_down.info())
	print (os.getcwd())
	cmlb=pd.read_hdf("fs1/cmlb5", key="mlb")
	features_combined = pd.concat([flanks_down, cmlb], axis=1, copy=False)
	return features_combined

def training1_fs2(fickett):
	if fickett==1:
		flanks_down = pd.read_csv("fs2/w_fickett/down_pp_fickett", index_col=0)
		print (f'Feature set 2 selected with fickett score:\n{flanks_down.info()}')
	elif fickett==0:
		flanks_down = pd.read_csv("fs2/wo_fickett/down_pp_no_fickett", index_col=0)
		print (f'Feature set 2 selected without fickett score:')
		print (flanks_down.info())
	cmlb=pd.read_hdf("fs1/cmlb5", key="mlb")
	features_combined = pd.concat([flanks_down, cmlb], axis=1, copy=False)
	return features_combined

def training1_fs3(fickett):
	if fickett==1:
		flanks_down = pd.read_csv("fs3/w_fickett/down_cons_fickett", index_col=0)
		print (f'Feature set 3 selected with fickett score:\n{flanks_down.info()}')
	elif fickett==0:
		flanks_down = pd.read_csv("fs3/wo_fickett/down_cons_no_fickett", index_col=0)
		print (f'Feature set 3 selected without fickett score:')
		print (flanks_down.info())
	cmlb=pd.read_hdf("fs1/cmlb5", key="mlb")
	features_combined = pd.concat([flanks_down, cmlb], axis=1, copy=False)
	return features_combined

