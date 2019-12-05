import pysam,gzip,  datetime, collections, csv, os, sys, copy, logging,argparse, filetype
import numpy as np
from Bio import SeqIO, Seq
from StringIO import StringIO
with open("gibbon.duplicates.fasta", "rb") as handle, open("gibbon.duplicates.n.txt", "wb") as out:
    for record in SeqIO.parse(handle, "fasta"):
            c=collections.Counter(str(record.seq).upper())
