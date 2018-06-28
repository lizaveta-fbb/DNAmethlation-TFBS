from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
from FisherExact import fisher_exact
from itertools import izip
import random
import sys
import pybedtools
from pybedtools import BedTool
import os
hg38 = pybedtools.example_filename('/home/besedina/data/genomes/hg38.fa')
one_side_window = 10000
methyl_values = open("/home/ivavlakul/_asc_mutations_methylation.2018/_hydroxm_data/liver_T4.bed", 'r')
methyl_values_dict = {}
for line in methyl_values:
    line = line.strip().split('\t')
    #methyl_values_dict[line[0]+"."+line[1]] = line[4]
    methyl_values_dict[line[0]+"."+line[1]] = line[3]
    print line[0]
    
CGI = pybedtools.BedTool('/home/besedina/data/CpGIslands/cpgIslandExt.bed')

def neighbour_methylation_intersect_with_CGI(TF,type_PWM,quality):
    inner_directory = "/home/besedina/analysis/TFs_Separated_MotifLists/"+TF+"/"+type_PWM
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
    outfile_CGI = open(inner_directory+"/summary_hydroxymethylation_window_WholeGenome"+"."+quality+".CGI.liver_T4.txt", 'w')
    outfile_notCGI = open(inner_directory+"/summary_hydroxymethylation_window_WholeGenome"+"."+quality+".notCGI.liver_T4.txt", 'w')
    print outfile_CGI
    print outfile_notCGI
    for record in SeqIO.parse(fasta_file, "fasta"):
        if coordinates[i][3] <= 0.0005:
            record.seq = record.seq.upper()
            word = record.seq[coordinates[i][0]:coordinates[i][1]+1]
            if word.count("CG") != 0:
                chr_ = record.id.split(":")[0]
                chr_pos = int(record.id.split(":")[1].split("-")[0])
                word_bed = pybedtools.BedTool(' '.join([chr_, str(chr_pos+coordinates[i][0]), str(chr_pos+coordinates[i][1])]), from_string=True)
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(word))] #looking on sense strand in ny case
                default_meth_value_none = True
                for pos in positions_of_CpG:
                    if chr_+"."+str(pos) in methyl_values_dict:
                        IsInCGI = str(CGI.intersect(word_bed)) !=""
                        if IsInCGI:
                            outfile_CGI.write("0"+"\t"+"0"+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        else:
                            outfile_notCGI.write("0"+"\t"+"0"+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        default_meth_value_none = False
                if default_meth_value_none:
                    i+=1
                    continue
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(word))]
                if coordinates[i][2] == "+":
                    unword1_coord = ' '.join([chr_, str(chr_pos+coordinates[i][0]-one_side_window), str(chr_pos+coordinates[i][0])])
                    a = pybedtools.BedTool(unword1_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword1 = open(a.seqfn).read().split("\n")[1]
                    unword2_coord = ' '.join([chr_, str(chr_pos+coordinates[i][1]+1), str(chr_pos+coordinates[i][1]+1+one_side_window)])
                    a = pybedtools.BedTool(unword2_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword2 = open(a.seqfn).read().split("\n")[1]
                else:
                    unword2_coord = ' '.join([chr_, str(chr_pos+coordinates[i][0]-one_side_window), str(chr_pos+coordinates[i][0])])
                    a = pybedtools.BedTool(unword2_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword2 = open(a.seqfn).read().split("\n")[1]
                    unword1_coord = ' '.join([chr_, str(chr_pos+coordinates[i][1]+1), str(chr_pos+coordinates[i][1]+1+one_side_window)])
                    a = pybedtools.BedTool(unword1_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword1 = open(a.seqfn).read().split("\n")[1]
                if unword2.count("CG") != 0:
                    if coordinates[i][2] == "+":
                        positions_of_CpG = [chr_pos+coordinates[i][1]+1+m.start() for m in re.finditer("CG", str(unword2))]
                        positions_of_CpG = sorted(positions_of_CpG)
                    else:
                        positions_of_CpG = [chr_pos+coordinates[i][0]-one_side_window+m.start() for m in re.finditer("CG", str(unword2))]
                        positions_of_CpG = sorted(positions_of_CpG, reverse=True)
                    cpg_num = 1
                    for pos in positions_of_CpG:
                        distance_to_motif = min(abs(pos - (chr_pos+coordinates[i][0])), abs(pos - (chr_pos+coordinates[i][1])))
                        if chr_+"."+str(pos) in methyl_values_dict:
                            if IsInCGI:
                                outfile_CGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                            else:
                                outfile_notCGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        else:
                            if IsInCGI:
                                outfile_CGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                            else:
                                outfile_notCGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                        cpg_num +=1                        
                if unword1.count("CG") != 0:
                    if coordinates[i][2] == "+":
                        positions_of_CpG = [chr_pos+coordinates[i][0]-one_side_window+m.start() for m in re.finditer("CG", str(unword1))]
                        positions_of_CpG = sorted(positions_of_CpG, reverse=True)
                    else:
                        positions_of_CpG = [chr_pos+coordinates[i][1]+1+m.start() for m in re.finditer("CG", str(unword1))]
                        positions_of_CpG = sorted(positions_of_CpG)
                    cpg_num = -1
                    for pos in positions_of_CpG:
                        distance_to_motif = min(abs(pos - (chr_pos+coordinates[i][0])), abs(pos - (chr_pos+coordinates[i][1])))
                        if chr_+"."+str(pos) in methyl_values_dict:
                            if IsInCGI:
                                outfile_CGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                            else:
                                outfile_notCGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        else:
                            if IsInCGI:
                                outfile_CGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                            else:
                                outfile_notCGI.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                        cpg_num -=1
        i+=1
    outfile_CGI.close()
    outfile_notCGI.close()
TF_methylation_intersect_with_CGI(TF)

def TF_methylation_intersect_with_CGI(TF):
    for matrix_type in ["H11MO.0.A", "H11MO.1.A", "H11MO.2.A", "H11MO.0.B",
                            "H11MO.1.B", "H11MO.2.B", "H11MO.0.C", "H11MO.1.C",
                            "H11MO.2.C", "H11MO.0.D", "H11MO.1.D", "H11MO.2.D"]:
        if os.path.isfile("/home/besedina/data/PWMs/HOCOMOCO/human/HOCOMOCOv11_full_pwm_HUMAN_mono/"+TF+"_HUMAN."+matrix_type+".pwm"):
            for quality in ["A", "B", "C"]:
                neighbour_methylation_intersect_with_CGI(TF,matrix_type,quality)
    




def neighbour_methylation(TF,type_PWM,quality):
    inner_directory = "/home/besedina/analysis/TFs_Separated_MotifLists/"+TF+"/"+type_PWM
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
    outfile = open(inner_directory+"/summary_hydroxymethylation_window_WholeGenome"+"."+quality+".txt", 'w')
    print outfile
    for record in SeqIO.parse(fasta_file, "fasta"):
        if coordinates[i][3] <= 0.0005:
            record.seq = record.seq.upper()
            if coordinates[i][2] == "+":
                word = record.seq[coordinates[i][0]:coordinates[i][1]+1]
            else:
                word = record.seq[coordinates[i][0]:coordinates[i][1]+1]
            if word.count("CG") != 0:
                chr_ = record.id.split(":")[0]
                chr_pos = int(record.id.split(":")[1].split("-")[0])
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(word))] #looking on sense strand in ny case
                default_meth_value_none = True
                for pos in positions_of_CpG:
                    if chr_+"."+str(pos) in methyl_values_dict:
                        outfile.write("0"+"\t"+"0"+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        default_meth_value_none = False
                if default_meth_value_none:
                    i+=1
                    continue
                positions_of_CpG = [chr_pos+coordinates[i][0]+m.start() for m in re.finditer("CG", str(word))]
                if coordinates[i][2] == "+":
                    unword1_coord = ' '.join([chr_, str(chr_pos+coordinates[i][0]-one_side_window), str(chr_pos+coordinates[i][0])])
                    a = pybedtools.BedTool(unword1_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword1 = open(a.seqfn).read().split("\n")[1]
                    unword2_coord = ' '.join([chr_, str(chr_pos+coordinates[i][1]+1), str(chr_pos+coordinates[i][1]+1+one_side_window)])
                    a = pybedtools.BedTool(unword2_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword2 = open(a.seqfn).read().split("\n")[1]
                else:
                    unword2_coord = ' '.join([chr_, str(chr_pos+coordinates[i][0]-one_side_window), str(chr_pos+coordinates[i][0])])
                    a = pybedtools.BedTool(unword2_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword2 = open(a.seqfn).read().split("\n")[1]
                    unword1_coord = ' '.join([chr_, str(chr_pos+coordinates[i][1]+1), str(chr_pos+coordinates[i][1]+1+one_side_window)])
                    a = pybedtools.BedTool(unword1_coord, from_string=True)
                    a = a.sequence(fi=hg38)
                    unword1 = open(a.seqfn).read().split("\n")[1]
                if unword2.count("CG") != 0:
                    if coordinates[i][2] == "+":
                        positions_of_CpG = [chr_pos+coordinates[i][1]+1+m.start() for m in re.finditer("CG", str(unword2))]
                        positions_of_CpG = sorted(positions_of_CpG)
                    else:
                        positions_of_CpG = [chr_pos+coordinates[i][0]-one_side_window+m.start() for m in re.finditer("CG", str(unword2))]
                        positions_of_CpG = sorted(positions_of_CpG, reverse=True)
                    cpg_num = 1
                    for pos in positions_of_CpG:
                        distance_to_motif = min(abs(pos - (chr_pos+coordinates[i][0])), abs(pos - (chr_pos+coordinates[i][1])))
                        if chr_+"."+str(pos) in methyl_values_dict:
                            outfile.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        else:
                            outfile.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                        cpg_num +=1                        
                if unword1.count("CG") != 0:
                    if coordinates[i][2] == "+":
                        positions_of_CpG = [chr_pos+coordinates[i][0]-one_side_window+m.start() for m in re.finditer("CG", str(unword1))]
                        positions_of_CpG = sorted(positions_of_CpG, reverse=True)
                    else:
                        positions_of_CpG = [chr_pos+coordinates[i][1]+1+m.start() for m in re.finditer("CG", str(unword1))]
                        positions_of_CpG = sorted(positions_of_CpG)
                    cpg_num = -1
                    for pos in positions_of_CpG:
                        distance_to_motif = min(abs(pos - (chr_pos+coordinates[i][0])), abs(pos - (chr_pos+coordinates[i][1])))
                        if chr_+"."+str(pos) in methyl_values_dict:
                            outfile.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+str(methyl_values_dict[chr_+"."+str(pos)])+"\n")
                        else:
                            outfile.write(str(cpg_num)+"\t"+str(distance_to_motif)+"\t"+"NA"+"\n")
                        cpg_num -=1
        i+=1
    outfile.close()






    
for matrix_type in ["H11MO.0.A", "H11MO.1.A", "H11MO.2.A", "H11MO.0.B",
                        "H11MO.1.B", "H11MO.2.B", "H11MO.0.C", "H11MO.1.C",
                        "H11MO.2.C", "H11MO.0.D", "H11MO.1.D", "H11MO.2.D"]:
        if os.path.isfile("/home/besedina/data/PWMs/HOCOMOCO/human/HOCOMOCOv11_full_pwm_HUMAN_mono/"+TF+"."+matrix_type+".pwm"):
            print matrix_type
import os.path
TF_file = open("/home/besedina/data/cistrome/ABC.human.txt", 'r')
TF_list = []
for line in TF_file:
    TF = line.strip().split("_")[0]
    TF_list.append(TF)


def TF_methylation(TF):
    for matrix_type in ["H11MO.0.A", "H11MO.1.A", "H11MO.2.A", "H11MO.0.B",
                            "H11MO.1.B", "H11MO.2.B", "H11MO.0.C", "H11MO.1.C",
                            "H11MO.2.C", "H11MO.0.D", "H11MO.1.D", "H11MO.2.D"]:
        if os.path.isfile("/home/besedina/data/PWMs/HOCOMOCO/human/HOCOMOCOv11_full_pwm_HUMAN_mono/"+TF+"_HUMAN."+matrix_type+".pwm"):
            for quality in ["A", "B", "C"]:
                neighbour_methylation(TF,matrix_type,quality)
    for matrix_type in ["HT-SELEX_1", "Methyl-HT-SELEX_1", "HT-SELEX_2", "Methyl-HT-SELEX_2", "HT-SELEX_3",
                            "Methyl-HT-SELEX_3", "HT-SELEX_4", "Methyl-HT-SELEX_4", "HT-SELEX_5", "Methyl-HT-SELEX_5"]:
        if os.path.isfile("/home/besedina/data/PWMs/in_Vitro/matrices/"+TF+"."+matrix_type+".pwm"):
            for quality in ["A", "B", "C"]:
                neighbour_methylation(TF,matrix_type,quality)



from joblib import Parallel, delayed
import multiprocessing

num_cores = multiprocessing.cpu_count()

results = Parallel(n_jobs=num_cores)(delayed(TF_methylation)(TF) for TF in TF_list)
