import numpy as np
import pandas as pd
import itertools
from math import log
table_parentness = pd.read_csv('Parentness_of_filtered_alg2.txt', header=None, sep='\t')
contigs = np.array(list(table_parentness[0]))

plants = {}
for i in range(1,51):
    plants["plant"+str(i)] = table_parentness[i]


def pair_rec_dist(pair_cont):
    AABB = 0.0
    AABb = 0.0
    AAbb = 0.0
    AaBB = 0.0
    AaBb = 0.0
    Aabb = 0.0
    aaBB = 0.0
    aaBb = 0.0
    aabb = 0.0
    f_cont = pair_cont[0]
    s_cont = pair_cont[1]
    f_ind = np.where(contigs==f_cont)[0][0]
    s_ind = np.where(contigs==s_cont)[0][0]
    for plant in plants:
        if plants[plant][f_ind] == 1:
            if plants[plant][s_ind] == 1:
                AABB +=1
            if plants[plant][s_ind] == 0:
                AABb +=1
            if plants[plant][s_ind] == -1:
                AAbb +=1
        if plants[plant][f_ind] == 0:
            if plants[plant][s_ind] == 1:
                AaBB +=1
            if plants[plant][s_ind] == 0:
                AaBb +=1
            if plants[plant][s_ind] == -1:
                Aabb +=1
        if plants[plant][f_ind] == -1:
            if plants[plant][s_ind] == 1:
                aaBB +=1 
            if plants[plant][s_ind] == 0:
                aaBb +=1
            if plants[plant][s_ind] == -1:
                aabb +=1
    SummaryCount =  float(2*(AABB+AABb+AAbb+AaBB+AaBb+Aabb+aaBB+aaBb+aabb))
    hetero = AaBb - (AABB + aabb)
    if hetero < 0:
        hetero = 0
    recombDistance = float(hetero+2*aaBB+2*AAbb+AaBB+Aabb+AABb+aaBb)/SummaryCount
    if recombDistance < 0.5:
        recombDistCorr = 0.25*log((1+2* recombDistance)/(1-2*recombDistance))
    else:
        recombDistCorr = 'INF'
    return recombDistCorr

def pair_dist(pair_list):   
    f_element = pair_list[0].split('\t')
    s_element = pair_list[1].split('\t')
    if len(f_element) > 1 or len(s_element) > 1:
        recdistlist = []
        for pair_cont in itertools.product(f_element, s_element):
            if pair_cont[0] != pair_cont[1]:
                recombDistCorr = pair_rec_dist(pair_cont)
            else:
                recombDistCorr = 0.0
            recdistlist.append(recombDistCorr)

        if 'INF' not in recdistlist: 
            distance = np.mean(recdistlist)
        else:
            distance = 'INF'
    else:
        distance = pair_rec_dist(pair_list)
    return distance


def pair_from_index(list, index):
    m = n = len(list) - 1
    while index >= n:
        index -= n
        n -= 1
    m -= n
    return list[m], list[m + index + 1]


def isTriangle(sides):
    smallest,medium,biggest = sorted(sides)
    return smallest+medium>=biggest

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
