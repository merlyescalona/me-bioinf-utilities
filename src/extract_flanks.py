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

PROGRAM_NAME="extract_flanks"
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

def getBed(bedfile):
    """ Reading bedfile to extract gap info"""
    APPLOGGER.info("Reading BED file ({})".format(bedfile))
    # print("Reading BED file ({})".format(bedfile))
    start=datetime.datetime.now()
    bed=dict()
    # BED Dictionary:
    # BED[SCAFFOLD]: GAPID,START,END,SIZE
    bedCounter=1
    with open(bedfile) as f:
        for line in f:
            bedline=line.strip().split()
            try:
                bed[bedline[0]]+=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
            except:
                bed[bedline[0]]=[(bedCounter,int(bedline[1]), int(bedline[2]), int(bedline[2])-int(bedline[1]))]     
            bedCounter+=1
        end=datetime.datetime.now()
    APPLOGGER.info("Done reading BED file > {}".format( end-start))
    # print("Done reading BED file > {}".format( end-start))
    return bed

##########################################################################################

def getRefereceScaffolds(REFFILE,bed):
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
        if (counter%10)==0:APPLOGGER.debug("Scaffold:\t{}".format(counter))
        # if (counter%10)==0: print("Scaffold:\t{}".format(counter))
        newid=record.id.split("|")[0]
        try:
            scaffolds[newid]={'id':counter,'seq':str(record.seq),'size':len(record.seq), 'flanks':[], 'numgaps':len(bed[newid])}  
        except:
            scaffolds[newid]={'id':counter,'seq':str(record.seq),'size':len(record.seq), 'flanks':[], 'numgaps':0}
        counter+=1
    handle.close()
    end=datetime.datetime.now()
    APPLOGGER.info("Done reading reference file ({}) > {}".format(REFFILE, end-start))
    # print("Done reading reference file ({}) > {}".format(REFFILE, end-start))
    return scaffolds
# ----------------------------------------------------------------------------------------

def getFlankTable(scaffolds, bed, READFLANK):
    APPLOGGER.info("Extracting flank regions...")
    # print("Extracting flank regions...")
    start=datetime.datetime.now()
    flankTable=[]
    flankCounter=1
    for g in bed:
        for gap in bed[g]:
            gapID=gap[0]
            scaffold=g
            GS=gap[1]
            GE=gap[2]
            FSGS=GS-READFLANK
            FEGS=GS
            FSGE=GE
            FEGE=GE+READFLANK
            if FSGS<0: FSGS=1
            if FEGE>scaffolds[scaffold]['size']: FEGE=scaffolds[scaffold]['size']-1
            flankTable+=[\
            (\
                gapID,\
                flankCounter,\
                scaffold,\
                FSGS,\
                FEGS,\
                GS,\
                GE\
            ),(\
                gapID,\
                flankCounter+1,\
                scaffold,\
                FSGE,\
                FEGE,\
                GS,\
                GE\
            )]
            flankCounter+=2
    end=datetime.datetime.now()
    APPLOGGER.info("Done with flanks extraction: > {}".format(end-start))
    # print("Done with flanks extraction: > {}".format(end-start))
    return flankTable

# ----------------------------------------------------------------------------------------

def run(REFFILE, BEDFILE,OUTPUT,READFLANK):
    # Get the gap info
    bed=getBed(BEDFILE)
    numGaps=len([j for i in bed.keys() for j in range(0,len(bed[i]))])
    scaffolds=getRefereceScaffolds(REFFILE,bed)
    flankTable=getFlankTable(scaffolds,bed,READFLANK)
    # Updating scaffolds data structure
    for f in flankTable:
        scaffolds[f[2]]['flanks']+=[f[1]]
    for f in flankTable:
        scaffolds[f[2]]['numgaps']=len(scaffolds[f[2]]['flanks'])/2
    ###############################################################################  
    for sc in scaffolds:
      FS=0;  FE=READFLANK
      if (READFLANK > scaffolds[sc]['size']): 
        FE=scaffolds[sc]['size']
      TAG="SS"
      SEQUENCE= SeqIO.SeqRecord(\
          seq= Seq.Seq(str(scaffolds[sc]['seq'][FS:FE]).upper()),\
          id="{}:{}-{}:{}_0".format(sc, FS,FE,TAG),\
          description="(Starting scaffold {})".format(READFLANK)\
      )
      FS=scaffolds[sc]['size']-READFLANK
      FE=scaffolds[sc]['size']
      TAG="SE"
      if (READFLANK > scaffolds[sc]['size']): 
        FS=0
      SEQUENCE2= SeqIO.SeqRecord(\
        seq= Seq.Seq(str(scaffolds[sc]['seq'][FS:FE]).upper()),\
        id="{}:{}-{}:{}_INF".format(sc, FS,FE,TAG),\
        description="(Ending flank {})".format(READFLANK)\
      )
      with open(OUTPUT,"a") as handle:
            _=SeqIO.write(SEQUENCE,handle,"fasta")
            _=SeqIO.write(SEQUENCE2,handle,"fasta")
    ###############################################################################
    # Generating a FASTA file with all the Flanks
    TAG=['GS', 'GE']
    for findex in range(0,len(flankTable)):
        FID=flankTable[findex][1]
        s=flankTable[findex][2]
        FS=flankTable[findex][3]
        FE=flankTable[findex][4]
        GS=flankTable[findex][5]
        GE=flankTable[findex][6]
        TAG="GS"

        if FID % 2 ==0: TAG="GE"
        APPLOGGER.debug("{}:{}-{}:{}-F{} (Flank={})".format(s, FS,FE,TAG,FID, READFLANK))
    
        if findex % 10 ==0 : APPLOGGER.debug("[{}] {}".format(datetime.datetime.now(),findex))
        SEQUENCE= SeqIO.SeqRecord(\
            seq= Seq.Seq(str(scaffolds[s]['seq'][FS:FE]).upper()),\
            id="{}:{}-{}:{}-F{} (Flank={})".format(s, GS,GE,TAG,FID, READFLANK),\
            description="{}-{}".format(FS,FE)\
        )
        with open(OUTPUT,"a") as handle:
            _=SeqIO.write(SEQUENCE,handle,"fasta")
    return 1, "Finished"

def main():
    parser = argparse.ArgumentParser(\
        prog=PROGRAM_NAME,\
        description='This scripts generates a FASTA file with multiple sequences.'
            'These sequences correspond to the flanking regions before and after a gap and the flanking zones at the beginning and end of a scaffold '
            'in a reference assembly file. The size of the flank is give by the user, as'
            ' well as the BED file describing the gaps and the reference assembly file '
            'from which the final FASTA output file will be generated.',\
            add_help=False\
        )
    requiredArgs=parser.add_argument_group("{0}Required arguments{1}".format("\033[1m","\033[0m"))
    requiredArgs.add_argument(\
            '-g','--genome',\
            metavar = '<assembly.fasta>',\
        type = str,\
        required = True,\
        help = 'Filepath of the reference genome file to be used.')
    requiredArgs.add_argument(\
        '-b','--bed',\
        metavar = '<gaps.bed>',\
        type = str,\
        required = True,\
        help = 'Filepath of the bed file describing the gaps of the genome')
    requiredArgs.add_argument(\
        '-o','--OUTPUT',\
        metavar = '<file_path>',\
        type = str,\
        required = True,\
        help = 'Output file path.')
    requiredArgs.add_argument(\
        '-f','--readflank',\
        metavar = 'FLANK_SIZE',\
        type = int,\
        help =  'Flank size to be used to select the reads that are in the '
                'surroundings of the gap and determine whether there are '
                'reads that span the g.ap or not.',\
        required = True,\
    )
    optionalArts=parser.add_argument_group("{0}Optional arguments{1}".format("\033[1m","\033[0m"))
    optionalArts.add_argument(\
        '-l','--log',\
        metavar = '<log_level>',\
        type = str,\
        choices = ["INFO","DEBUG"],\
        default= "INFO",\
        help = 'Verbosity levels.')
    informationGroup= parser.add_argument_group("{0}Information arguments{1}".format("\033[1m","\033[0m"))
    informationGroup.add_argument('-v', '--version',\
    action='version',\
    version='{0}: Version {1}.{2}.{3}'.format(PROGRAM_NAME,VERSION,MIN_VERSION,FIX_VERSION),\
    help="Show program's version number and exit")
    informationGroup.add_argument('-h', '--help',\
    action='store_true',\
    help="Show this help message and exit")
    args = parser.parse_args()  
    # Checking files
    if args.log == 'DEBUG':
        APPLOGGER.setLevel(logging.DEBUG)
    else:
        APPLOGGER.setLevel(logging.INFO)
    if (os.path.exists(os.path.abspath(args.genome))):
        REFFILE=os.path.abspath(args.genome)
    else:
        message="Reference file does not exist. Please verify. Exiiting."
        return (-1, message)
    if (os.path.exists(os.path.abspath(args.bed))):
        BEDFILE=os.path.abspath(args.bed)
    else:
        message="BED file does not exist. Please verify. Exiiting."
        return (-1, message)
    OUTPUT=os.path.abspath(args.OUTPUT)
    if (os.path.exists(OUTPUT)):
        message="Output file exist. Not overwriting. Please verify. Exiiting."
        return (-1, message)
    READFLANK=int(args.readflank)
    return run(REFFILE, BEDFILE,OUTPUT,READFLANK)

if __name__ == "__main__":
    SIGNAL_EXIT,MESSAGE=main()
    if SIGNAL_EXIT == -1:
        APPLOGGER.error(MESSAGE)
    else:
        APPLOGGER.info(MESSAGE)
    sys.exit(SIGNAL_EXIT)
