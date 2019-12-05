import argparse,collections,copy,csv,datetime,filetype
import glob,gzip,logging,os,pysam,random,string,sys
import numpy as np
from Bio import SeqIO, Seq, AlignIO,Alphabet
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline, MafftCommandline
import Bio.Align.AlignInfo
from Bio import pairwise2

#=====================================================================================        
# Parameters

tsvfile="ncbi_contamination/gibbon.13.ncbi.contamination.scaffolds.totrim.txt"
genomefile="assemblies/gibbon.15.fasta"

#=====================================================================================        
# read assembly
filekind=filetype.guess(genomefile)
handle=None
if filekind and filekind.extension in ['gz','GZ']:
    handle=gzip.open(genomefile, 'rb') 
else:
    handle=open(genomefile, "rt")

genome=SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
handle.close()
# len(genome) # 18834

#=====================================================================================        
# Read totrim file
totrim=[]
handle=None
if filekind and filekind.extension in ['gz','GZ']:
    handle=gzip.open(tsvfile, 'rb') 
else:
    handle=open(tsvfile, "rt")

for line in handle.readlines():
    l=line.strip().split()
    if not l[0] =="Sequence":
        totrim+=[l]

handle.close()
#=====================================================================================        
sequences=[]

#=====================================================================================        

for line in range(0, len(totrim)):
    scaffold=totrim[line][0]
    # col1= len
    start=int(totrim[line][2])
    end=int(totrim[line][3])
    sequences+=[genome[scaffold][start:end]]

with open("ncbi_contamination/gibbon.13.sequences.to.trim.fasta", "wb") as handle:
    for i in range(0,len(sequences)):
        _=SeqIO.write(sequences[i], handle, "fasta")