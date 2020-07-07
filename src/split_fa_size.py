
#!/usr/bin/python
import pysam,gzip,  argparse,datetime, collections, csv, os, sys, filetype, logging
import numpy as np
from Bio import SeqIO

# ========================================================================================
PROGRAM_NAME="split_fa_size"
VERSION=1
MIN_VERSION=0
FIX_VERSION=0
# ========================================================================================
APPLOGGER = logging.getLogger(PROGRAM_NAME)
APPLOGGER.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
loggerFormatter = logging.Formatter(\
  fmt = "[%(asctime)s]  %(levelname)s (%(funcName)s|%(lineno)d): %(message)s",\
  datefmt = "%d/%m/%Y %I:%M:%S %p")
ch.setFormatter(loggerFormatter)
APPLOGGER.addHandler(ch)

def run(parameters):
  message="All correct."
  APPLOGGER.info("Reading file...")
  handle=None
  filekind=filetype.guess(parameters['input'])
  if filekind and filekind.extension in ['gz','GZ','gZ','Gz']:
    APPLOGGER.info("Running gzip...")
    handle=gzip.open(parameters['input'], "rt")
  else:
    handle=open(parameters['input'], "rt")
  APPLOGGER.info("Creating output file...")
  out=open(parameters['output'], "wb")
  APPLOGGER.info("Parsing FASTA file...")
  for record in SeqIO.parse(handle, "fasta"):
    chunks, chunk_size = len(record.seq), len(record.seq)/parameters['size']
    print(chunks,chunk_size)
    subseqs=[ record.seq[i:i+chunk_size] for i in range(0, chunks, chunk_size) ]
    num_digits=len(str(len(subseqs)))
    for index in range(0,len(subseqs)):
        seq=SeqIO.SeqRecord(\
            seq=Seq.Seq(subseqs[index]),\
            id="{0}_{1:0{2}d}".format(record.id, index, num_digits),\
            description=""
        )
        SeqIO.write(seq,out, "fasta")
  handle.close()
  out.close()
  APPLOGGER.info("Closing files...")
  return True, message



def main():
  parser = argparse.ArgumentParser(\
  prog=PROGRAM_NAME,\
    description='Splits sequences of a FASTA file in sequences of a given length',\
    add_help=False\
  )
  requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
  requiredArgs.add_argument(\
    '-i','--input',\
    metavar = '<sequences.fasta>',\
    type = str,\
    required = True,\
    help = 'Filepath of the reference genome file to be used.')
  requiredArgs.add_argument(\
    '-o','--output',\
    metavar = '<file_path>',\
    type = str,\
    required = True,\
    help = 'Output file path.')
  requiredArgs.add_argument(\
    '-s','--size',\
    metavar = 'INT',\
    type = int,\
    help = 'Length of final subsequences.')
  optionalArts=parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
  optionalArts.add_argument(\
    '-l','--log',\
    metavar = '<log_level>',\
    type = str,\
    choices = ['INFO','DEBUG'],\
    default= 'INFO',\
    help = 'Verbosity levels.')
  informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
  informationGroup.add_argument('-v', '--version',\
    action='version',\
    version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
    help="Show program's version number and exit")
  informationGroup.add_argument('-h', '--help',\
    action='store_true',\
    help="Show this help message and exit")
  try:
    args = parser.parse_args()  
    APPLOGGER.debug(args)
    parameters=dict()
    if args.help:
      return(0,"Finished")
    # Checking files
    if (os.path.exists(os.path.abspath(args.input))):
      REFFILE=os.path.abspath(args.input)
    else:
      message="Input file does not exist. Please verify. Exiiting."
      return (-1, message)
    OUTPUT=os.path.abspath(args.output)
    if (os.path.exists(OUTPUT)):
      message="Output file exist. Not overwriting. Please verify. Exiiting."
      return (-1, message)
    parameters['size']=int(args.size)
    parameters['output']=OUTPUT
    parameters['input']=REFFILE
    return run(parameters)
  except:
    message="Parsing error. Exiting"
    parser.print_help()
    return (-1, message)

if __name__ == "__main__":
  SIGNAL_EXIT,MESSAGE=main()
  if SIGNAL_EXIT == -1:
    APPLOGGER.error(MESSAGE)
  else:
    APPLOGGER.info(MESSAGE)
  sys.exit(SIGNAL_EXIT)


