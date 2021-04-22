import numpy as np
from tqdm import tqdm
import pandas as pd
from os import path, makedirs
import pdb
from .FileIOUtils import ReadFasta, SaveTableToFile, LoadTableFromFile
from .DeBruijnGraph import *


def CreateKmersFromReads(reads, kmersize):
    """
    Function to create kmers of a given size based on a dictionary of reads
    Inputs;
        reads (dict): reads
        kmersize (int): The size a kmer will be
    Output:
        kmers (list): All kmers (non-unique) existing in the reads.
    """

    kmers = []
    for read in reads.values():
        idx = 0;
        for idx in range(len(read)-(kmersize-1)):
            #don't care about dupes, just want a list of all k-1 mers
            kmers += [read[idx:kmersize+idx]]
    return kmers

def GetReverseComplement(seq, pairing_dict={'G':'C','C':'G','A':'T','T':'A'}):
    """
    Function to extract the reverse compliment from a DNA sequence based on a standardized pairing dictionary.
    Inputs:
        seq (str): A DNA sequence
        pairing_dict (dict): A dictionary of complementary base pairs.
    Outputs:
        new_seq (str): The reverse complement to the input sequence.
    """
    new_seq = ''
    for s in seq[::-1]:
        if s in pairing_dict:
            new_seq += pairing_dict[s]
        else:
            print('Letter not found in alphabet!')
    return new_seq

def CreateDBGraph(kmer_table, kmersize):
    """
    Function to generate a De Bruijn Graph and save it to a file.
    Inputs:
        kmer_table (list): All kmers found within the reads, as created by CreateKmersFromReads
        kmersize (int): Value to use as the kmersize in the graph.
    Output:
        node_table (dict): The constructed node_table as a dictionary with nodes of k-1 mer and edges connecting them.
    """
    
    print('Creating graph with the following attributes:\n\t# of kmers: %d\n\tkmersize: %d\n\n'%(len(kmer_table),kmersize))
    node_table = {}
    #read and parse all our k-mers
    for direction in ['forward','reverse']:
        for kmer in tqdm(kmer_table):
            #account for the reverse direction!
            kmer_degree = kmer_table[kmer]
            if direction == 'reverse':
                kmer = GetReverseComplement(kmer)
            kmer_1 = kmer[:-1]
            kmer_2 = kmer[1:]
            #draw all our edges!
            if kmer_1 in node_table:
                currentNode = node_table[kmer_1]
            else:
                currentNode = GraphNode(seq = kmer_1)
            #add our next node
            append_to_node(currentNode,kmer[-1],'next',count=kmer_degree)
            #add the previous node!
            if kmer_2 in node_table:
                nextNode = node_table[kmer_2]
            else:
                nextNode = GraphNode(seq = kmer_2)
            append_to_node(nextNode,kmer[0],'prev',count=kmer_degree)
    
            #update them in our node table!
            node_table[kmer_1] = currentNode
            node_table[kmer_2] = nextNode
    return node_table

def CorrectErrors(kmer_table, threshold = 2): 
    """
    Function to correct errors in a kmer table based on a specified threshold. Will change each kmer-1 with a degree < threshold to the closest kmer with degree > threshold.
    Inputs:
        kmer_table (dict): A dictionary with key=kmer sequence and value = degree
        threshold (int): A threshold of degree each kmer should be found in the reads to be considered without errors.
    Outputs:
        kmer_table (dict): A kmer_table with each node < threshold degree having been corrected
    """

    #find the degree of each node in the table
    low_degree_nodes = []
    seq_loop = 0
    for seq in kmer_table:
        seq_loop+=1
        if seq_loop % int(len(kmer_table)/10) == 0:
            print('Successfully analyzed %d out of %d kmers.'%(seq_loop, len(kmer_table)))
        # find where these are less than or equal to the threshold
        if kmer_table[seq] <= threshold:
            low_degree_nodes+=[seq]
            corrected_kmer = FindCorrectedKmer(seq,kmer_table, threshold = threshold)
            #delete only if we can correct it
            if corrected_kmer is not None:
                kmer_table[corrected_kmer] += kmer_table[seq]
    for seq in low_degree_nodes:
        #delete the node from our table
        del kmer_table[seq]
    return kmer_table

def FindCorrectedKmer(seq, kmer_table, threshold = 4, alphabet = ['G','C','A','T'], depth = True):
    """
    Function to find the corrected Kmer based on the threshold and given kmer.
    Inputs:
        seq (str): The sequence to correct.
        kmer_table (dict): A dictionary of kmers and their degrees to find a correction in.
        threshold (int): The degree threshold that a correction is found.
        alphabet (list): All possible letters within the sequence.
    Output:
        corrected_kmer (str): The corrected kmer
    """
    for idx in range(len(seq)-1,0,-1):
        options = [seq[0:idx] + a + seq[idx+1:] for a in alphabet]
        found_options = [o for o in options if o in kmer_table]
        found_options_better_than_thresh = np.where(np.array([kmer_table[x] for x in found_options])> threshold)[0]
        if len(found_options_better_than_thresh) > 0:
            #randomly correct here
            corrected_kmer = found_options[np.random.choice(found_options_better_than_thresh)]
            return corrected_kmer
        #hamming distance of 2 now -- using the "depth" parameter so we don't recurse forever
        elif len(found_options)>0 and depth:
            corrected_kmer_dist2_options = [FindCorrectedKmer(opt,kmer_table, depth = False) for opt in found_options]
            if len(corrected_kmer_dist2_options)>0:
                corrected_kmer = np.random.choice(corrected_kmer_dist2_options)
                return corrected_kmer

    #no fix found...
    return None

def GenerateHistograms(reads, kmer_table,output_folder):
    """
    Function to generate the histograms of kmer degree and read length.
    Inputs:
        reads (dict): Dictionary of reads.
        kmer_table (dict): Dictionary of kmer sequences and degrees.
        output_folder (str): Folder to write the histograms to.
    Outputs:
        Saved dataframes of read length and kmer degrees.
    """
    if not path.isdir(output_folder):
        makedirs(output_folder)
    #now let's make some histogram data!
    #read length histograms
    read_lengths = [len(x) for x in reads.values()]
    pd.DataFrame(np.histogram(read_lengths,bins=1000),index=['Counts','Lengths']).T.to_csv(path.join(output_folder,'ReadLengthHistogram.csv'),index=False)
    #kmer histograms
    kmer_degrees = [x for x in kmer_table.values()]
    pd.DataFrame(np.histogram(kmer_degrees,bins=1000),index=['Counts','Degree']).T.to_csv(path.join(output_folder,'KmerDegreeHistogram.csv'),index=False)

def CreateDeBruijnGraphAndSave(reads_fasta, kmersize, de_bruijn_graph_filename, existing_kmer_table = None, do_error_correction = True, generate_histograms = False, save_kmer_table = True):
    """
    Function to generate a De Bruijn graph from a reads fasta file with a given kmer size.
    Inputs:
        reads_fasta (str): FASTA file containing the reads to generate a graph from.
        kmersize (int): The size to use for kmers.
        de_bruijn_graph_filename (str): The filename to write the DBG to.
        existing_kmer_table (None, str): If not None, will generate a graph using an existing dict of kmer sequence and degrees.
        do_error_correction (bool): To perform error correction.
        generate_histograms (bool): To generate histograms of read length and kmer degree.
        save_kmer_table (bool): To save the kmer_table used to generate the De Bruijn Graph.
    Output:
        node_table (dict): Generated De Bruijn graph.
    """

    reads = ReadFasta(reads_fasta)
    smallest_read_length = np.min([len(r) for r in reads.values()])

    if kmersize > smallest_read_length:
        raise ValueError('Chosen kmer size is larger than the smallest read. Please lower the kmer size to at most %d'%smallest_read_length)

    if existing_kmer_table is not None:
        print('Loading in existing kmer table...')
        kmer_table = LoadTableFromFile(existing_kmer_table)
    else:
        print('Parsing kmers of size %d from reads.\n'%kmersize)
        kmers = CreateKmersFromReads(reads,kmersize)
        
        kmers_un, degrees = np.unique(kmers,return_counts = True)
        kmer_table = dict(zip(kmers_un, degrees))
        if do_error_correction is True:
            print('Correcting errors in kmers. This will take a while (multiple minutes)...')
            kmer_table = CorrectErrors(kmer_table)
        
        if save_kmer_table is True:
            ext = 'dat' #must always be dat
            kmer_table_path = path.join(path.dirname(de_bruijn_graph_filename),'kmer_table_%s%s'%('corrected' if do_error_correction else 'uncorrected',ext))
            SaveTableToFile(kmer_table, kmer_table_path)
            print('Successfully wrote kmer table to file: %s'%kmer_table_path)

    if generate_histograms:
        GenerateHistograms(reads, kmer_table, path.dirname(de_bruijn_graph_filename))
        
    # now let's construct the graph from the table!
    node_table = CreateDBGraph(kmer_table, kmersize)
    SaveTableToFile(node_table, de_bruijn_graph_filename)
    print('\nDone creating graph!\nSuccessfully saved to file %s'%de_bruijn_graph_filename)

    return node_table
