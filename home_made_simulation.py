import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from random import randint

# sys.argv[1] -> counts
# sys.argv[2] -> transcript fasta
# sys.argv[3] -> simulated fasta

def seqread(record):
    start = randint(0,len(record.seq) -33)
    end = start + (26 + randint(0,6))
    seq = record.seq[start:end]
    return "N" not in str(seq),seq

dictcount = {}

with open(sys.argv[1]) as fi:
    for line in fi:
        if not line.startswith("#") and not line.startswith("Geneid"):
            key = line.split("\t")[0].rstrip() #TRANSCRIPT ID 
            count = line.split("\t")[-1].rstrip() #COUNT
            dictcount[key] = int(count)

simulation = []
seq = ""
name = "AB"
c = 0
for record in SeqIO.parse(sys.argv[2], "fasta"):
    if record.id.split(":")[0] in dictcount:
        c+=1
        i = 0
        while i < (dictcount[record.id.split(":")[0]]):
            flag = False
            while not flag:
                flag,seq = seqread(record)
            read = SeqRecord(seq,name+"_" +str(c)+"_"+str(i),"","")
            simulation.append(read)
            i += 1
SeqIO.write(simulation, sys.argv[3], "fasta")

