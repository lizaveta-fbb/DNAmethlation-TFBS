from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
from FisherExact import fisher_exact
from itertools import izip
import random
import sys
TF = sys.argv[1]
type_PWM = sys.argv[2]
quality = sys.argv[3]
#"HT-SELEX_1"
#type_PWM = "H11MO.2.C"

GFF_annotation = open("/home/besedina/data/genomes/Annotation/Homo_sapiens.GRCh38.78.gtf", 'r')
gene_start_dict = {}
for line in GFF_annotation:
    line = line.strip().split('\t')
    if line[0] in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]:
        if line[1] == "ensembl_havana" or "ensembl" or "havana":
            if line[2] == "gene":
                chrom = "chr"+line[0]
                if chrom not in gene_start_dict:
                    gene_start_dict[chrom] = {}
                    gene_start_dict[chrom]["+"] = []
                    gene_start_dict[chrom]["-"] = []
                if line[6] =="+":
                    if int(line[3]) not in gene_start_dict[chrom]["+"]:
                        gene_start_dict[chrom]["+"].append(int(line[3]))
                else:
                    if int(line[4]) not in gene_start_dict[chrom]["-"]:
                        gene_start_dict[chrom]["-"].append(int(line[4]))
    
import os
if not os.path.exists("/home/besedina/analysis/TFs_Separated_MotifLists/"+TF):
    os.makedirs("/home/besedina/analysis/TFs_Separated_MotifLists/"+TF)
inner_directory = "/home/besedina/analysis/TFs_Separated_MotifLists/"+TF+"/"+type_PWM
if not os.path.exists(inner_directory):
    os.makedirs(inner_directory)


if "H11MO." in type_PWM:
    ##HOCOMOCO
    lengthPWM = sum(1 for line in open('/home/besedina/data/PWMs/HOCOMOCO/human/HOCOMOCOv11_full_pwm_HUMAN_mono/'+TF+'_HUMAN.'+type_PWM+'.pwm'))-1
    log_file = open("/home/besedina/analysis/TFs_logs/"+TF+"_HUMAN."+type_PWM+".full_"+quality+".log", 'r')
    pval_file = open("/home/besedina/analysis/TFs_pvals/"+TF+"_HUMAN."+type_PWM+".full_"+quality+".pval", 'r')
else:
    ##SELEX
    lengthPWM = sum(1 for line in open('/home/besedina/data/PWMs/in_Vitro/matrices/'+TF+'.'+type_PWM+'.pwm'))
    log_file = open("/home/besedina/analysis/TFs_logs/"+TF+"."+type_PWM+"."+quality+".log", 'r')
    pval_file = open("/home/besedina/analysis/TFs_pvals/"+TF+"."+type_PWM+"."+quality+".pval", 'r')



inMotif_CpG_distances = []


methyl_values = open("/home/besedina/data/cistrome_methylome_intersection_whole/"+TF+"_HUMAN."+quality+".methylome_threshold19.bed", 'r')
methyl_values_dict = {}
for line in methyl_values:
    line = line.strip().split('\t')
    methyl_values_dict[line[0]+"."+line[1]] = line[4]

    

coordinates = []

for line_log, line_pval in izip(log_file, pval_file):
    pval = float(line_pval.strip())
    line_log = line_log.strip().split("\t")
    #if pval <= 0.0005:
    strand = line_log[2]
    start = int(line_log[1])
    end = start+lengthPWM-1
    coordinates.append([start, end, strand, pval])
    
i = 0
if quality == "C":
    fasta_file = "/home/besedina/data/cistrome/fasta/medium/human/"+TF+"_HUMAN."+quality+".fasta"
    bed_file = open("/home/alla/GTRD/12_07_17_GTRD/small3/medium/"+TF+"_HUMAN."+quality+".bed",'r')
elif quality == "B":
    fasta_file = "/home/besedina/data/cistrome/fasta/high/human/"+TF+"_HUMAN."+quality+".fasta"
    bed_file = open("/home/alla/GTRD/12_07_17_GTRD/small3/high/"+TF+"_HUMAN."+quality+".bed",'r')
else:
    fasta_file = "/home/besedina/data/cistrome/fasta/highest/human/"+TF+"_HUMAN."+quality+".fasta"
    bed_file = open("/home/alla/GTRD/12_07_17_GTRD/small3/highest/"+TF+"_HUMAN."+quality+".bed",'r')
    
print fasta_file
    
for record in SeqIO.parse(fasta_file, "fasta"):
    if coordinates[i][3] <= 0.0005:
        record.seq = record.seq.upper()
        if coordinates[i][2] == "+":
            word = record.seq[coordinates[i][0]:coordinates[i][1]+1]
            unword1 = record.seq[0:coordinates[i][0]]
            unword2 = record.seq[coordinates[i][1]+1:]
        else:
            word = record.seq[coordinates[i][0]:coordinates[i][1]+1].reverse_complement()
            unword1 = record.seq[0:coordinates[i][0]]
            unword2 = record.seq[coordinates[i][1]+1:]
        if word.count("CG") != 0:
            chr_ = record.id.split(":")[0]
            chr_pos = int(record.id.split(":")[1].split("-")[0])
            values = []
            if coordinates[i][2] == "+":
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(word))]
                for pos in positions_of_CpG:
                    if chr_+"."+str(pos) in methyl_values_dict:
                        starts = []
                        for start in gene_start_dict[chr_]["+"]:
                            if start-10000 <= pos <= start:
                                starts.append(abs(start-pos))
                                if abs(start-pos) >10000:
                                    print "+"
                        if len(starts) >=1:
                            inMotif_CpG_distances.append(random.choice(starts))
                        else:
                            inMotif_CpG_distances.append("NA")
                        del methyl_values_dict[chr_+"."+str(pos)]
            elif coordinates[i][2] == "-":
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(record.seq[coordinates[i][0]:coordinates[i][1]+1]))] #looking on sense strand
                for pos in positions_of_CpG:
                    if chr_+"."+str(pos) in methyl_values_dict:
                        starts = []
                        for start in gene_start_dict[chr_]["-"]:
                            if start+10000 >= pos >=start:
                                starts.append(abs(start-pos))
                                if abs(start-pos) >10000:
                                    print "-"
                        if len(starts) >=1:
                            inMotif_CpG_distances.append(random.choice(starts))
                        else:
                            inMotif_CpG_distances.append("NA")         
                        del methyl_values_dict[chr_+"."+str(pos)]
    i+=1



inMotif_CpG_distances_file = open(inner_directory+"/inMotif_CpG_distances"+"."+quality+".txt", 'w')
for value in inMotif_CpG_distances:
    inMotif_CpG_distances_file.write(str(value)+'\n')
inMotif_CpG_distances_file.close()


outMotif_CpG_distances = []
i=0
for line in bed_file:
    line = line.strip().split('\t')
    chr_ = line[0]
    start = int(line[1])
    end = int(line[2])
    for key in methyl_values_dict.keys():
        key_=key.split('.')
        chrom=key_[0]
        pos=int(key_[1])
        if chrom == chr_:
            if start<=pos<=end:
                if coordinates[i][2] == "+":
                    dists = []
                    for TSS in gene_start_dict[chr_]["+"]:
                        if TSS-10000 <= pos <= TSS:
                            dists.append(abs(TSS-pos))
                            if abs(TSS-pos) >10000:
                                print "+"
                    if len(dists) >=1:
                            outMotif_CpG_distances.append(random.choice(dists))
                    else:
                            outMotif_CpG_distances.append("NA")
                    del methyl_values_dict[chr_+"."+str(pos)]
                elif coordinates[i][2] == "-":
                    dists = []
                    for TSS in gene_start_dict[chr_]["-"]:
                        if TSS+10000 >= pos >=TSS:
                            dists.append(abs(TSS-pos))
                            if abs(TSS-pos) >10000:
                                print "-"
                    if len(dists) >=1:
                            outMotif_CpG_distances.append(random.choice(dists))
                    else:
                            outMotif_CpG_distances.append("NA")         
                    del methyl_values_dict[key]
    i+=1
                    

outMotif_CpG_distances_file = open(inner_directory+"/outMotif_CpG_distances"+"."+quality+".txt", 'w')
for value in outMotif_CpG_distances:
    outMotif_CpG_distances_file.write(str(value)+'\n')
outMotif_CpG_distances_file.close()

