#!/bin/bash/python3

from Bio import SeqIO
import argparse,numpy, os, sys, copy

# ========================================================================================

def read_genome_sizes(filename):
    scaffold_sizes=dict()
    with open(filename, "r") as handle:
        for line in handle.readlines():
            scaff=line.strip().split()[0]
            size=int(line.strip().split()[1])
            scaffold_sizes[scaff]=size
    return scaffold_sizes

# ========================================================================================

def bin_scaffold(scaffold_sizes, binsize):
    bin_scaffolds=dict()
    for scaffold in scaffold_sizes.keys():
        bin_scaffolds[scaffold]=range(0,scaffold_sizes[scaffold], binsize)
    return bin_scaffolds

# ========================================================================================

def read_tad(filename,binsize):
    tads=dict()
    with open(filename, "r") as handle:
        for line in handle.readlines():
            scaff=line.strip().split()[0]
            start=int(line.strip().split()[1])
            end=int(line.strip().split()[2])
            try:
                if not newboundary in tmp[scaff]:
                    tads[scaff]+=[start,end]
            except:
                tads[scaff]=[start, end]
    for scaffold in tads.keys():
        tads[scaffold]=[val for val in numpy.unique(tads[scaffold])]
    return tads

# ========================================================================================

def translate_boundary(binned_scaffold_sizes, scaffold, bound):
    val=0
    for index in range(0,len(binned_scaffold_sizes[scaffold])): 
        x=binned_scaffold_sizes[scaffold][index]
        if bound < x: 
            val=binned_scaffold_sizes[scaffold][index-1]
            break
    return val

# ========================================================================================

allboundaries=dict()

with open("mmul10.TADs.all.bed", "r") as handle:
    for line in handle.readlines():
        scaffold=line.strip().split()[0]
        start=int(line.strip().split()[1])
        end=int(line.strip().split()[2])
        try:
            if not newboundary in tmp[scaffold]:
                allboundaries[scaff]+=[start,end]
        except:
            allboundaries[scaffold]=[start, end]

for s in allboundaries.keys():
    allboundaries[s]=[ ss for ss in numpy.unique(allboundaries[s])]

with open("mmul10.boundaries.all.bed", "w") as f:
    for s in allboundaries.keys():
        for i in range(0, len(allboundaries[s])):
            f.write("{}\t{}\t{}\n".format(s,allboundaries[s][i],allboundaries[s][i]+1))

# ========================================================================================

scaffold_sizes_filename="/scratch7/RheMac/assemblies/mmul10/mmul.chrom.size"
scaffold_sizes=read_genome_sizes(scaffold_sizes_filename)

# ========================================================================================

tads_filename_1="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.1K.bedpe"
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.2K.bedpe"
highRes=2000
tad1=read_tad(tads_filename_1,1000)
tad2=read_tad(tads_filename_2,2000)
tmp=tad2
bins_lower_res=bin_scaffold(scaffold_sizes, 1000)
bins_higher_res=bin_scaffold(scaffold_sizes, 2000)
resolutions={ s:[highRes]*len(tad2[s]) for s in tad2.keys() }
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]




print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.2500.bedpe"
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,2500)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, 2500)
highRes=2500
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.3K.bedpe"
highRes=3000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.4K.bedpe"
highRes=4000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.7K.bedpe"
highRes=7000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.8K.bedpe"
highRes=8000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.9K.bedpe"
highRes=9000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.10K.bedpe"
highRes=10000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.20K.bedpe"
highRes=20000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.25K.bedpe"
highRes=25000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.30K.bedpe"
highRes=30000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.40K.bedpe"
highRes=40000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.50K.bedpe"
highRes=50000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.60K.bedpe"
highRes=60000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.70K.bedpe"
highRes=70000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.80K.bedpe"
highRes=80000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.90K.bedpe"
highRes=90000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]




print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.100K.bedpe"
highRes=100000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]





print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.200K.bedpe"
highRes=200000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.250K.bedpe"
highRes=250000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.300K.bedpe"
highRes=300000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.400K.bedpe"
highRes=400000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]


print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.500K.bedpe"
highRes=500000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.600K.bedpe"
highRes=600000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]

print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.700K.bedpe"
highRes=700000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.800K.bedpe"
highRes=800000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]



print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.900K.bedpe"
highRes=900000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]




print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.1M.bedpe"
highRes=1000000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]




print(len([b for t in tmp for b in tmp[t]]))
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.2M.bedpe"
highRes=2000000
tad1=copy.copy(tmp)
tad2=read_tad(tads_filename_2,highRes)
tmp=tad2
bins_higher_res=bin_scaffold(scaffold_sizes, highRes)
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]

with open("mmul10.highres.filtered.bed", "w") as f:
    for s in tmp.keys():
        for i in range(0, len(tmp[s])):
            f.write("{}\t{}\t{}\n".format(s,tmp[s][i],tmp[s][i]+1))





with open("mmul10.highres.txt", "w") as f:
    for s in tmp.keys():
        for i in range(0, len(tmp[s])):
            f.write("{}\t{}\t{}\n".format(s,tmp[s][i], resolutions[s][i]))


print(len([b for t in tmp for b in tmp[t]]))

tmp2=dict()
for scaffold in tmp:
    for bound in tmp[scaffold]:
        tmp2[scaffold]=[ b for b in numpy.unique(tmp[scaffold])]



with open("mmul10.highres.filtered.bed", "w") as f:
    for s in tmp2.keys():
        for i in range(0, len(tmp2[s])):
            f.write("{}\t{}\t{}\n".format(s,tmp2[s][i],tmp2[s][i]+1))


print(len([b for t in tmp for b in tmp2[t]]))


######################################################################################################


tads_filename_1="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.1K.bedpe"
tads_filename_2="/scratch7/RheMac/tad_calling/bed_contacts/mmul10.TADs.2K.bedpe"
lowhRes=1000
highRes=2000
tad1=read_tad(tads_filename_1,1000)
tad2=read_tad(tads_filename_2,2000)
tmp=tad1
bins_lower_res=bin_scaffold(scaffold_sizes, 1000)
bins_higher_res=bin_scaffold(scaffold_sizes, 2000)
resolutions={ s:[highRes]*len(tad2[s]) for s in tad2.keys() }
for scaffold in tad1.keys():
    for boundary in tad1[scaffold]:
        newboundary=translate_boundary(bins_higher_res, scaffold, boundary)
        print("({}) {} => {}".format(scaffold, boundary, newboundary))
        try:
            if not newboundary in tmp[scaffold]:
                tmp[scaffold]+=[newboundary]
                resolutions[scaffold]+=[highRes]
        except:
            tmp[scaffold]=[newboundary]
            resolutions[scaffold]=[highRes]




# # 12734
# if __name__ == "__main__":
#     SIGNAL_EXIT, MESSAGE=main()
#     if SIGNAL_EXIT == -1:
#         APPLOGGER.error(MESSAGE)
#     else:
#         APPLOGGER.info(MESSAGE)
#     sys.exit(SIGNAL_EXIT)

