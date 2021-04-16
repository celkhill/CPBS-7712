import numpy as np
import re
import pandas as pd
from os import path,makedirs
try:
    import ujson as json
except:
    import json
import pdb
from DBGContigMaker import ParseJson, ParseReads, ReadFasta

def GetHammingDistance(string1,string2):
    return np.sum([int(x!=y) for x,y in zip(string1,string2)])

def CorrectErrors(node_table):
    alphabet = ['G','C','T','A']
    #find the degree of each node in the table
    low_degree_nodes = []
    threshold = 4
    for seq in node_table:
        # find where these are less than or equal to the threshold
        if len(node_table[seq].next) <= threshold:
            low_degree_nodes+=[seq]
            for idx in range(len(seq)-1,0,-1):
                options = [seq[0:idx] + a + seq[idx+1:] for a in alphabet]
                found_options = [o for o in options if o in node_table]
                found_options_better_than_thresh = np.where(np.array([len(node_table[x].next) for x in found_options])> threshold)[0]
                if len(found_options_better_than_thresh) > 0:
                    corrected_kmer = found_options[np.random.choice(found_options_better_than_thresh)]
                    pdb.set_trace()
        #look for another node close in hamming distance
        
    pdb.set_trace() 

existing_graph_filename = '../analysis/kmer_13/DeBruijnGraph_kmer13.json'
if not path.isfile(existing_graph_filename):
    raise ValueError('Existing De Bruijn Graph file not found. Please check the following argument:\n\t%s'%existing_graph_filename)
print('\nLoading existing graph found in %s\n\n'%existing_graph_filename)
#loading in our node table (DBG)
if existing_graph_filename.endswith('.json'):
    with open(existing_graph_filename, 'r') as fp:
        node_table,kmersize = ParseJson(json.load(fp))#we can't use object hook here because ujson doesn't support it
elif existing_graph_filename.endswith('.dat'):
    with open(existing_graph_filename, 'rb') as fp:
        node_table = pickle.load(fp)
        kmersize = int(node_table['kmersize'])
node_table_kmersize = len(node_table[list(node_table.keys())[0]].seq)+1

reads_fasta = '../data/READS.fasta'
reads = ParseReads(ReadFasta(reads_fasta))

CorrectErrors(node_table)
    

