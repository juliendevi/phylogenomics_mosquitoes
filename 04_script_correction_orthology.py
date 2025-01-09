#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:19:51 2020

@author: mg16464
"""

import os 
from ete3 import Tree
from collections import defaultdict
import statistics
from Bio import SeqIO
import sys

my_cdir = os.chdir(sys.argv[1])
my_dir = os.listdir(my_cdir)

my_trees = [file for file in my_dir if file.endswith('.rooted')]  ### Change here the suffix

### Create a dictionary with the long branches

dict_long_branches = defaultdict(list)

for trees in my_trees:
    current_tree = Tree(trees)
    branch_lengths = [node.dist for node in current_tree.traverse('preorder')]
    mean_bl = statistics.mean(branch_lengths)
    sd_bl = statistics.stdev(branch_lengths)
    threshold = (sd_bl * 2)+ mean_bl             ### Change here the SD treshold!!!

    for node in current_tree.traverse('preorder'):
        if node.dist > threshold:
            leaves = node.get_leaf_names()
            if len(leaves) < 0.4 * len(current_tree.get_leaf_names()):
                dict_long_branches[trees.replace('.fasta.fasta.treefile.rooted','.fasta')].append(leaves) ### Change here the suffix
            elif len(leaves) > 0.4 * len(current_tree.get_leaf_names()):
                dict_long_branches[trees.replace('.fasta.fasta.treefile.rooted','.fasta')].append(['Problematic'])

### break the lists of list in dict values (try to come up with a anew method)

dict_long_branches_single =  defaultdict(list)

for k,v in dict_long_branches.items():
    for j in v:
        for m in j:
            dict_long_branches_single[k].append(m)

### open the fasta using the keys in the dict

                                
for k,v in dict_long_branches_single.items():
    if 'Problematic' in v:
        print(k)
    else:
        print(v, file=open(k + '_sequences_removed.txt', 'a'))
        with open(k, 'r') as my_fasta:
            for record in SeqIO.parse(my_fasta, "fasta"):
                if record.id not in v:
                    print('>' + record.id + '\n' + record.seq, file=open(k + '.SD2.new', 'a'))
