import os
import itertools
import subprocess


filenames = os.listdir("/home/besedina/data/PWMs/in_Vitro/matrices")
filenames = sorted(filenames)
new_filenames = []
d= {}
#i = 0
for filenam in filenames:
    filename = filenam.split(".")
    name = filename[0]
    selex_type = filename[1][0:-2]
    if name not in d:
        d[name] = {}
    if selex_type not in d[name]:
        d[name][selex_type] = []
        d[name][selex_type].append(1)
    else:
        d[name][selex_type].append(d[name][selex_type][-1]+1)



EvalSimilarity = 'java -cp /home/besedina/programs/macroape/ape-19aug2017.jar ru.autosome.macroape.di.EvalSimilarity'
Info_TF = {}
MethPairNonmeth = {}
for TF in d:
    print TF
    m = 0
    n = 0
    if "Methyl-HT-SELEX" in d[TF]:
        meth = d[TF]["Methyl-HT-SELEX"]
        m = 1
    if "HT-SELEX" in d[TF]:
        non = d[TF]["HT-SELEX"]
        n = 1
    if m == 1 and n == 1:
        if len(meth)==1 and len(non)==1:
            print 'mn=1'
            pwm_meth = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.Methyl-HT-SELEX_1.pwm'
            pwm_non = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.HT-SELEX_1.pwm'
            result = subprocess.check_output([EvalSimilarity+' '+pwm_meth+' '+pwm_non+' --first-from-mono --second-from-mono'], shell=True)
            result = result.split('\n')
            distance = float(result[22].split('\t')[1])
            max_meth_non_distance = distance
            min_meth_non_distance = distance
            max_meth_meth_distance = 0
            max_non_non_distance = 0
        else:
            if len(meth)>1:
                m = len(meth)
                methset = True
                distances_meth = []
                for pair in itertools.combinations(meth,2):
                    pwm1 = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.Methyl-HT-SELEX_'+str(pair[0])+'.pwm'
                    pwm2 = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.Methyl-HT-SELEX_'+str(pair[1])+'.pwm'
                    result = subprocess.check_output([EvalSimilarity+' '+pwm1+' '+pwm2+' --first-from-mono --second-from-mono'], shell=True)
                    result = result.split('\n')
                    distance = float(result[22].split('\t')[1])
                    distances_meth.append(distance)
                    max_meth_meth_distance = max(distances_meth)
            else:
                max_meth_meth_distance = 0
            if len(non) >1:
                n = len(non)
                nonset = True
                distances_non = []
                for pair in itertools.combinations(non,2):
                    pwm1 = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.HT-SELEX_'+str(pair[0])+'.pwm'
                    pwm2 = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.HT-SELEX_'+str(pair[1])+'.pwm'
                    result = subprocess.check_output([EvalSimilarity+' '+pwm1+' '+pwm2+' --first-from-mono --second-from-mono'], shell=True)
                    result = result.split('\n')
                    distance = float(result[22].split('\t')[1])
                    distances_non.append(distance)
                    max_non_non_distance = max(distances_non)
            else:
                max_non_non_distance = 0
            dist_meth_non_list = []
            for i in range(1, m+1):
                pwm_meth = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.Methyl-HT-SELEX_'+str(i)+'.pwm'
                for j in range(1, n+1):
                    pwm_non = '/home/besedina/data/PWMs/in_Vitro/matrices/'+str(TF)+'.HT-SELEX_'+str(j)+'.pwm'
                    result = subprocess.check_output([EvalSimilarity+' '+pwm_meth+' '+pwm_non+' --first-from-mono --second-from-mono'], shell=True)
                    result = result.split("\n")
                    distance_meth_non = float(result[22].split('\t')[1])
                    dist_meth_non_list.append(distance_meth_non)
            max_meth_non_distance = max(dist_meth_non_list)
            min_meth_non_distance = min(dist_meth_non_list)
    Info_TF[TF] = [m,n, max_meth_non_distance,min_meth_non_distance,max_meth_meth_distance,max_non_non_distance]                   
                    
            

pairs_meth_non = open('/home/besedina/data/PWMs/in_Vitro/pairs_maxBetween_minBetween_maxMeth_maxNon_distance.txt', 'w')
for TF in Info_TF:
    pairs_meth_non.write(TF+"\t"+str(Info_TF[TF][0])+"\t"+str(Info_TF[TF][1])+'\t'+str(Info_TF[TF][2])+'\t'+
                         str(Info_TF[TF][3])+'\t'+str(Info_TF[TF][4])+'\t'+str(Info_TF[TF][5])+"\n")
pairs_meth_non.close()
