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
tsvfile="ncbi_contamination/gibbon.13.ncbi.contamination.scaffolds.duplicates.txt"
genomefile="assemblies/gibbon.14.fasta"
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
# Read duplicates file
duplicates=[]
handle=None
if filekind and filekind.extension in ['gz','GZ']:
    handle=gzip.open(tsvfile, 'rb') 
else:
    handle=open(tsvfile, "rt")

for line in handle.readlines():
    duplicates+=[line.strip().split()]

handle.close()
#=====================================================================================        
to_remove=[]
to_keep=[]
#=====================================================================================        
for line in range(0, len(duplicates)):
    scaffolds=duplicates[line]
    sequences=[ genome[sc] for sc in scaffolds]
    if len(sequences)==2:
        aln=pairwise2.align.globalxx(sequences[0].seq, sequences[1].seq)
        lenseq=len(sequences[0].seq)
        score=int(aln[0][2])
        matches=0
        for index in range(0,lenseq):
            if sequences[0].seq[index]==sequences[1].seq[index]:
                matches+=1
        if matches == lenseq:
            to_remove+=[ sequences[1].id ]
            to_keep+=[ sequences[0].id ]
            # del genome[sequences[1].id]
        print("[{}] {} vs {}\tScore:{}\t Length:{}\tMatches:{}".format(\
            line,\
            sequences[0].id,\
            sequences[1].id,\
            score,\
            lenseq,\
            matches))
    if len(sequences)>2:
        alnfile="/scratch7/Gibbon/ncbi_contamination/msa/line_{}.fasta".format(line)
        msafile="/scratch7/Gibbon/ncbi_contamination/msa/line_{}.msa.fasta".format(line)
        # Write sequences to file 
        with open(alnfile,  "wb") as out:
            _=SeqIO.write(sequences, out, "fasta")
        # align
        aligner_cline = MuscleCommandline(input=alnfile)
        try:
            stdout, stderr = aligner_cline()
            aligner="MUSCLE"
        except:
            aligner_cline=MafftCommandline(input=alnfile)
            stdout, stderr = aligner_cline()
            aligner="MAFFT"
        with open(msafile,"w") as handle: 
            handle.write(stdout)
        # get similarity
        align=AlignIO.read(msafile, "fasta")
        summary=Bio.Align.AlignInfo.SummaryInfo(align)
        cnss=Bio.Align.AlignInfo.SummaryInfo.gap_consensus(summary,\
            threshold=0.1,\
            ambiguous='N',\
            require_multiple=2)
        ncount=np.sum(np.array([i for i in cnss])=='N')
        gapcount=np.sum(np.array([i for i in cnss])=='-')
        if ncount > 0 or gapcount >0:
            to_keep+=[ i.id for i in sequences]
        else:
            to_keep+=[sequences[0].id]
            to_remove+=[ i.id for i in sequences if not i.id in to_keep]
        print("[{}] {}\tLength:{}\tN:{}\tGaps:{}".format(\
            line,\
            [sequences[i].id for i in range(0,len(sequences))],\
            len(cnss),\
            ncount,\
            gapcount))
         

with open("ncbi_contamination/gibbon.13.ncbi.contamination.scaffolds.duplicates_removed.txt", "wb") as handle:
    for i in to_remove:
        handle.write("{}\t{}\n".format(i, len(genome[i].seq)))

with open("ncbi_contamination/gibbon.13.ncbi.contamination.scaffolds.duplicates_kept.txt", "wb") as handle:
    for i in to_keep:
        handle.write("{}\t{}\n".format(i, len(genome[i].seq)))


newversion={i:genome[i] for i in genome if not i in to_remove}

with open("assemblies/gibbon.15.fasta", "wb") as handle:
    for i in genome:
        if not i in to_remove:
            _=SeqIO.write(genome[i], handle, "fasta")