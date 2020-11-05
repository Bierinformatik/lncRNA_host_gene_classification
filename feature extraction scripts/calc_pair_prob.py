#calculates pairing probability 
#sets sliding window lengths

import sys,os,inspect
import traceback as tb
import logging
import RNA
import argparse
from io import StringIO
import time
import math
import gzip
import importlib
import traceback as tb
#from Randseq import createrandseq
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq
#numpy and matplolib and pyplot
import numpy as np
from collections import Counter

def fold_up():
    try:
        logid = 'pair_prob.fold: '
        #set path for VRNA lib if necessary
        #if 'vrna' in kwargs:
            #sys.path=[kwargs['vrna']] + sys.path
        try:
            global RNA
            RNA = importlib.import_module('RNA')
            globals().update(
            {n: getattr(RNA, n) for n in RNA.__all__}
                if hasattr(RNA, '__all__')
                else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')
                })

            md = RNA.md()
            md = None

        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
            log.error(logid+''.join(tbe.format()))
            sys.exit()

        md = RNA.md()
        #set the sliding window length
        window = 60
        span = window
        region = 1
        md.max_bp_span = span # kwargs['span']
        md.window_size = window #kwargs['window']
        #for the negative lncs
        out = os.path.abspath(os.path.join("..", "training_1","sno","pp_60"))
        print (out)
        if not os.path.isdir(out):
            os.mkdir(out)
        #f = os.path.join("..", "fas", "Merged_Filtered_Closest_Exon_lnc.fa")
        files = [f for f in os.listdir(os.path.join("..", "training_1","sno")) 
                 if os.path.isfile(os.path.join("..", "training_1","sno",f))]
        fa_name = [] 
        fa_real = ""
        fa_c=0
        lnc=0
        #print (files[0])
        for f in files:
            for fa in SeqIO.parse(open (os.path.join("..", "training_1","sno",f), "r"), "fasta"):
                print (fa.id)
                #line = fa.id.split(":")
            	#name_ex = line[0][line[0].rindex("|")+1:]
                #name_ex = line[0]#for lnc    
                name_ex=fa.id #for pp from windows             
                if not name_ex in fa_name:
                    fa_name.append(name_ex)
                    fa_c=0
                    fa_real=name_ex+".0"		
                else:
                    c = Counter(fa_name)
                    fa_c = c[name_ex]+1
                    fa_name.append(name_ex)		
                    fa_real=name_ex+"."+str(fa_c)
                fa_real=fa.id #for pp from windows                
                print (fa_real)			
                seq = str(fa.seq).upper()
                #print(seq)
            	
		#for whole exons, check >60, for windows everything below happens if >100
		#to avoid overlaps
         
                if len(seq) > 99:
                    lnc+=1
                    print (len(seq))
                    window = 60
                else:
                    window = len(seq)
                window = 60
                span = window
                region = 1            
                md.max_bp_span = span # kwargs['span']
                md.window_size = window #kwargs['window']            
            # create new fold_compound object
                fc = RNA.fold_compound(str(seq), md, RNA.OPTION_WINDOW)
            # call prop window calculation
                data = { 'up': [] }
                fc.probs_window(region, RNA.PROBS_WINDOW_UP, up_callback, data)
            
                write_up(seq=str(seq), data=data, region=int(region), window=str(window), span=str(span), outdir=out, name=fa_real)
        print(lnc)  

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        #log.error(logid+''.join(tbe.format()))
        sys.exit()

def up_callback(v, v_size, i, maxsize, what, data):

    logid = 'pair_prob.up_callback: '
    try:
        if what & RNA.PROBS_WINDOW_UP:
            data['up'].extend([v])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        #log.error(logid+''.join(tbe.format()))

#def fetch():
def print_up(data=None, seqlength=None, region=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    #logid = scriptname+'.print_up: '
    try:
        if data:
            ups=''
            for i in range(int(seqlength)):
                for x in range(1,region+1):
                    if isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
            return ups
        else:
            print('No up data to print')
            return ups
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        #clog.error(logid+''.join(tbe.format()))

def isinvalid(x=None):
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def write_up(seq, data, region, window, span, outdir, name):

    try:
        #print(print_up(data['up'],len(seq),region))
        with open(os.path.join(outdir, name+".pp"), "w") as w:
            w.write(print_up(data['up'],len(seq),region))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
    return 1

if __name__ == '__main__':

    logid = 'pair_prob.main: '
    try:
        #args=parseargs()
        #if args.loglevel != 'WARNING':
          #log = setup_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)

        #log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores')
        #preprocess(sequence=args.sequence, window=args.window, span=args.span, region=args.region, outdir=args.outdir, fold=args.fold)
        fold_up()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        #log.error(logid+''.join(tbe.format()))
