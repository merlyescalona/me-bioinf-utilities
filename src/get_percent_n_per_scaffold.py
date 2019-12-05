#! /usr/bin/python
"""
extract.flank.sequences.from.assembly
Flanks around gaps and beginning and end of scaffolds
Merly Escalona <mmescalo@ucsc.edu>
Tue Oct  1 15:18:59 PDT 2019
"""
# Extract flank sequences to aligned them to the Human genome.
import pysam,gzip,  datetime, collections, csv, os, sys, copy, logging,argparse, filetype
import numpy as np
from Bio import SeqIO, Seq
from StringIO import StringIO

PROGRAM_NAME="get_percent_n_per_scaffold"
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
# ========================================================================================
APPLOGGER = logging.getLogger(PROGRAM_NAME)
APPLOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
loggerFormatter = logging.Formatter(\
    fmt = "[%(asctime)s]  %(levelname)s: %(message)s",\
    datefmt = "%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)

##########################################################################################
def run(REFFILE):
  """ Get reference file info"""
  APPLOGGER.info("Reading scaffolds from reference file...")
  # print("Reading scaffolds from reference file...")
  start=datetime.datetime.now()
  scaffolds=dict()
  counter=0
  # ScpqPdg
  filekind=filetype.guess(REFFILE)
  handle=None
  if filekind and filekind.extension in ['gz','GZ']:
      handle=gzip.open(REFFILE, 'rb') 
  else:
      handle=open(REFFILE, "rt")
  for record in SeqIO.parse(handle, "fasta"):
    c=collections.Counter(str(record.seq).upper())
    print("{}\t{}".format(record.id, (c["N"]*1.0)/len(record.seq) ))
  handle.close()
     

def main():
    parser = argparse.ArgumentParser(\
        prog=PROGRAM_NAME,\
        description='',\
            add_help=False\
        )
    requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
    requiredArgs.add_argument(\
            '-g','--genome',\
            metavar = '<assembly.fasta>',\
        type = str,\
        required = True,\
        help = 'Filepath of the reference genome file to be used.')
    informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
    informationGroup.add_argument('-v', '--version',\
    action='version',\
    version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
    help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
    action='store_true',\
    help="Show this help message and exit")
    args = parser.parse_args()  
    if (os.path.exists(os.path.abspath(args.genome))):
        REFFILE=os.path.abspath(args.genome)
    else:
        message="Reference file does not exist. Please verify. Exiiting."
        return (-1, message)
    return run(REFFILE)

if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)
