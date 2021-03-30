import numpy as np
import re
import pandas as pd
from os import path,makedirs
import matplotlib.pyplot as plt
try:
    import ujson as json
except:
    import json
import pickle
import pdb
import argparse
from types import SimpleNamespace


def ReadFasta(filename):
    with open(filename,'r') as f:
        data = f.readlines()
    return data

def ParseReads(data):
    reads = {}
    for num in range(len(data)):
        if data[num].startswith('>'):
            key = data[num].rstrip() # remove newline char
            val = data[num+1].rstrip() # remove newline char
            reads[key] = val
    return reads

class GraphNode:
    def __init__(self, seq='', next = [], prev = [], visited = []):
        self.seq = seq
        self.next = np.array(next)
        self.prev = np.array(prev)
        self.visited = np.array(visited)

def construct_sequence_from_nodes(node_list):
    out_sequence = ''.join([n.seq[0] for n in node_list[0:-1]])+node_list[-1].seq
    return out_sequence

def RemoveEdge(edge_list,edge_to_remove):
    ind = np.where(edge_list == edge_to_remove)[0]
    if len(ind)==0:
        raise ValueError('Edge %s not found in list! Please check.')
    else:
        edge_list = np.delete(edge_list,ind[-1])
    return edge_list

def SelectNextNode(current_node,node_table, mode = 'forwards'):
    if mode == 'forwards':
        options = np.unique(current_node.next)
    elif mode == 'reverse':
        options = np.unique(current_node.prev)
    if len(options)>1:
        #node degree of all our possible next nodes
        possible_next_node_degrees = []
        for x in options:
            try:
                if mode == 'forwards':
                    next_node = node_table[current_node.seq[1:]+x]
                if mode =='reverse':
                    next_node = node_table[x+current_node.seq[:-1]]
            except KeyError:
                raise KeyError('Node not found in node table! You messed something up somewhere when making the DBG.')
            if mode == 'forwards':
                possible_next_node_degrees.append(len(next_node.next))
            if mode == 'reverse':
                possible_next_node_degrees.append(len(next_node.prev))
        #if none of them have any edges, just pick one at random, we will stop here anyways!
        possible_next_node_degrees = np.array(possible_next_node_degrees)
        if not np.any(possible_next_node_degrees>0):
            next_na = np.random.choice(options)
        else:
            #try to find ones that aren't bridges--have a degree of more than 1
            try:
                next_na = np.random.choice(options[np.where(possible_next_node_degrees>1)])
            #if they are all bridges, then just pick any that aren't zero
            except ValueError:
                next_na = np.random.choice(options[np.where(possible_next_node_degrees>0)])
    elif len(options) == 1:
        next_na = options[0]
    else:
        next_na = False
    # print('%s out of %s'%(next_na,' '.join(options)))
    return next_na

def WriteFasta(sequence, filename):
    with open(filename,'w') as f:
        f.write('>FinalContig:Length:%s'%len(sequence))
        f.write('\n')
        f.write(sequence)
        f.close()
    return

def CreateGraphAndSave(reads,kmersize,outname = '../data/DeBruijneGraph'):
    
    #read and parse all our k-mers
    node_table = {}
    read_idx = 0
    for read in reads.values():
        idx = 0;
        previousNode = None
        while idx<len(read)-kmersize+1:
            kmer = read[idx:kmersize+idx]
            kmer_1 = kmer[0:-1]
            kmer_2 = kmer[-1]
            if kmer_1 in node_table:#grab the existing node
                currentNode = node_table[kmer_1]
            else:#create new node
                currentNode = GraphNode(kmer_1)
            #add our next node
            currentNode.next = np.append(currentNode.next, kmer_2)
            #add the previous node if its here
            if previousNode:
                currentNode.prev = np.append(currentNode.prev, previousNode)
            #save it to our node table!
            node_table[kmer_1] = currentNode
            #set our previous node for reversal of the graph
            previousNode = kmer_1[0]
            idx+=1
        #now create the last node!
        last_kmer = kmer_1[1:]+ kmer_2
        if not last_kmer in node_table:
            ending_node = GraphNode(last_kmer)
            ending_node.prev = np.append(ending_node.prev, kmer_1[0])
            node_table[last_kmer] = ending_node
    
    node_table['kmersize']= kmersize
    print('Writing out DBG to file %s'%outname)
    if outname.endswith('.json'):
        with open(outname, 'w') as fp:
            json.dump(node_table,fp)
    elif outname.endswith('.dat'):
        with open(outname, 'wb') as fp:
            pickle.dump(node_table,fp)

# solve for our query sequence
def FindQueryInGraph(query_seq,kmersize,node_table):
    idx = 0
    solved_nodes = []
    removed_edges = []
    while idx<len(query_seq)-kmersize+2:#+2 because we need the end of the sequece, i.e. we have to dial ONE. NUMBER. HIGHER.
        kmer_1 = query_seq[idx:kmersize+idx-1]
        try:
            cn = node_table[kmer_1]
        except KeyError:
            raise ValueError('Kmer not found %s'%kmer_1)
        #check that our next node solution is valid!
        if len(query_seq)>(kmersize+idx):
            next_na = query_seq[kmersize+idx-1] 
            if not next_na in cn.next:
                raise ValueError('No solution found! Stopping at this sequence: %s'%construct_sequence_from_nodes(solved_nodes))
            else:
                #remove the edge!
                RemoveEdge(cn.next,next_na)
        solved_nodes.append(cn)    
        idx+=1
    return solved_nodes

#extend our solution!
def CreateContigFromGraph(solved_nodes,node_table,direction = 'forwards'):
    if len(solved_nodes)>0:
        cn = solved_nodes[-1]
        # print('Extending solution %s...'%direction)
    else:#pick a random starting node TODO: make this better
        cn = np.random.choie(list(node_table.keys()))
    keep_looping = True
    while keep_looping:
        next_na = SelectNextNode(cn, node_table, mode = direction)
        #if we can't find a next node, stop looping
        if next_na is False:
            keep_looping = False
        #if we can, look for it in the node table
        else:
            try:
                if direction == 'forwards':
                    next_node = node_table[cn.seq[1:]+next_na]
                    cn.next = RemoveEdge(cn.next,next_na)
                    cn.visited = np.append(cn.visited,next_na)
                elif direction == 'reverse':
                    next_node = node_table[next_na + cn.seq[:-1]]
                    cn.prev = RemoveEdge(cn.prev,next_na)
                    cn.visited = np.append(cn.visited,'-%s'%next_na)
                solved_nodes.append(next_node)
                cn = next_node
            except KeyError:
                if direction == 'forwards':
                    seq = cn.seq[1:]+next_na
                elif direction == 'reverse':
                    seq = next_na + cn.seq[:-1]
                print('Sequence %s not found in the node table!'%seq)
                keep_looping = False
        if len(solved_nodes) == 5000:
            keep_looping = False
    return solved_nodes

def ConstructArgs():
    defaults = {'output_folder':'../analysis/',
            'output_filetype': 'json',
            'kmersize':27,
            'plot_histograms':False,
            'existing_graph_filename':None,
            'solving_iterations':10,
            'existing_allele': None}

    parser = argparse.ArgumentParser(description='Here are the arguments for the query sequence finder tool:')

    req_group = parser.add_argument_group('Required Arguments','You must specify these arguments.')
    opt_group = parser.add_argument_group('Optional Arguments')

    req_group.add_argument('-reads_fasta',action='store',dest='reads_fasta',
            help='Fasta containing the reads to construct our De Bruijn Graph')
    req_group.add_argument('-query_fasta',action='store',dest='query_fasta',
            help='Fasta containing the sequence to search for within the De Bruijn Graph.')
    opt_group.add_argument('-output_folder',action='store',dest = 'output_folder',
            help="Output folder where the results will be saved to.",default = defaults['output_folder'])
    opt_group.add_argument('-output_filetype',action='store',dest = 'output_filetype',
            help = 'Filetype to save out the graph and winning contig node objects to. Options are json or dat', default = defaults['output_filetype'])
    opt_group.add_argument('-kmersize',action = 'store',dest = 'kmersize',
            type=int,help = 'Kmersize the De Bruijn Graph will be constructed with.', default = defaults['kmersize'])
    opt_group.add_argument('-plot_histograms',action='store_true', dest = 'plot_histograms',
            help = 'Option to generate histograms of the read length and kmer edge degree.', default= defaults['plot_histograms'])
    opt_group.add_argument('-existing_graph_filename',action = 'store',dest = 'existing_graph_filename',
            help='Argument to use an existing graph if one has already been created.',default= defaults['existing_graph_filename'])
    opt_group.add_argument('-solving_iterations',action = 'store',dest='solving_iterations',
            type = int, help='Number of solving iterations to try before selecting the longest contig found.', default = defaults['solving_iterations'])
    opt_group.add_argument('-existing_allele', action = 'store',dest = 'existing_allele',
            help='Argument to use an existing solved allele that has already been created.', default = defaults['existing_allele'])

    return parser

def ParseJson(jdata):
    if type(jdata) == dict:
        #delete the kmersize so we can just get the node table
        kmersize = int(jdata['kmersize'])
        del jdata['kmersize']
        #convert to GraphNode class from jdata dict
        return {x: GraphNode(**jdata[x]) for x in jdata}, kmersize
    elif type(jdata)==list:
        return [GraphNode(**x) for x in jdata]

def GenerateOutputTable(final_sequence, reads, query_seq):
    output_table = []
    for read in reads:
        for match in re.finditer(reads[read],final_sequence):
            query_match = re.search(reads[read],query_seq)
            if query_match:
                qstart = query_match.start()
                qend = query_match.end()
            else:
                qstart = np.nan
                qend = np.nan
            output_table.append({'sseqid':read,
                'sstart':match.start(),
                'send':match.end(),
                'qstart': qstart,
                'qend': qend})

    return output_table

def ResetNodeTable(node_table):
    for key in node_table:
        for edge in node_table[key].visited:
            if '-' in edge:
                node_table[key].prev = np.append(node_table[key].prev,edge.replace('-',''))
            else:
                node_table[key].next = np.append(node_table[key].next,edge)
        node_table[key].visited = np.array([])
    return node_table

if __name__ == "__main__":
    #create the arg parser
    parser = ConstructArgs()

    #read in the results from parse args, pass into a dict containing our variables
    settings = parser.parse_args()

    reads = ParseReads(ReadFasta(settings.reads_fasta))

    #parse our reads and figure out our kmer size
    kmersize = settings.kmersize
    smallest_read_length = np.min([len(r) for r in reads.values()])
    if kmersize > smallest_read_length:
        raise ValueError('Chosen kmer size is larger than the smallest read. Please lower the kmer size to at most %d'%smallest_read_length)
    query_seq = list(ParseReads(ReadFasta(settings.query_fasta)).values())[0]
    if not path.isdir(settings.output_folder): makedirs(settings.output_folder)

    #if we just want to use the existing allele, do it right now!
    if settings.existing_allele is not None:
        print('Solved allele file found! Generating output table...')
        final_sequence = list(ParseReads(ReadFasta(settings.existing_allele)).values())[0]
    else:

        #load an existing graph if we have one, otherwise make a new one!
        if settings.existing_graph_filename is None:
            print('Creating graph with the following attributes:\n\t# of reads: %d\n\tkmersize: %d\n\n'%(len(reads),kmersize))
            outname = path.join(settings.output_folder,'DeBruijneGraph.'+ settings.output_filetype)
            CreateGraphAndSave(reads,kmersize, outname)
            print('\nDone creating graph!\n')
        else:
            outname = settings.existing_graph_filename
            if not path.isfile(outname):
                raise ValueError('Existing De Bruijn Graph file not found. Please check the following argument:\n\t%s'%outname)
            print('\nLoading existing graph found in %s\n\n'%outname)
            #loading in our node table (DBG)
            if outname.endswith('.json'):
                with open(outname, 'r') as fp:
                    node_table,kmersize = ParseJson(json.load(fp))#we can't use object hook here because ujson doesn't support it
            elif outname.endswith('.dat'):
                with open(outname, 'rb') as fp:
                    node_table = pickle.load(fp)
                    kmersize = int(node_table['kmersize'])
            node_table_kmersize = len(node_table[list(node_table.keys())[0]].seq)+1
            if kmersize != node_table_kmersize:
                if idx == 0:
                    print('Specified kmer size does not agree with the node table sequence lengths. Using kmer size of %d'%node_table_kmersize)
                kmersize = node_table_kmersize
            print('Graph loaded!\n\n')
    
        #now let's make some histogram data!
        if settings.plot_histograms:
            #read length histograms
            read_lengths = [len(x) for x in reads.values]
            pd.DataFrame(zip(np.histogram(read_lengths,bins=100)),columns=['Counts','Lengths']).to_csv(path.join(settings.output_folder,'ReadLengthHistogram.csv'),index=False)
            #kmer histograms
            kmer_degrees = [len(x.next) for x in node_table.values()]
            pd.DataFrame(zip(np.histogram(kmer_degrees,bins=100)),columns=['Counts','Degree']).to_csv(path.join(settings.output_folder,'KmerDegreeHistogram.csv'),index=False)
    
        lengths = []
        best_solved_nodes = ['One']
        max_iter = settings.solving_iterations
        #now let's solve!
        for idx in range(max_iter):
            if idx >0:
                node_table = ResetNodeTable(node_table)
            print('Solving iteration %d out of %d'%(idx+1,max_iter))
            # read in our query
            solved_nodes = FindQueryInGraph(query_seq,kmersize,node_table)
            forwards_contig = CreateContigFromGraph(solved_nodes,node_table,direction = 'forwards')
            #now let's go in reverse!
            all_solved_nodes = CreateContigFromGraph(forwards_contig[::-1],node_table,direction = 'reverse')[::-1] #reverse the node list back!
            if len(all_solved_nodes)> len(best_solved_nodes):
                best_solved_nodes = all_solved_nodes
            print('Solved sequence found of length: %d'%(len(all_solved_nodes)+1))
            lengths.append(len(all_solved_nodes)+1)
    
    
        if settings.output_filetype == 'json':
            json.dump(best_solved_nodes,open(path.join(settings.output_folder,'solved_contig.json'),'w'),indent = 2)
            best_solved_nodes = ParseJson(json.load(open(path.join(settings.output_folder,'solved_contig.json'),'r')))
        elif settings.output_filetype == 'dat':
            pickle.dump(best_solved_nodes,open(path.join(settings.output_folder,'solved_contig.dat'),'wb'))
            best_solved_nodes = pickle.load(open(path.join(settings.output_folder,'solved_contig.dat'),'rb'))
    
        final_sequence = construct_sequence_from_nodes(best_solved_nodes)
        print('\n\nWinning sequence length: %d \nWinning sequence: %s'%(len(final_sequence),final_sequence))
        print('\n\nWriting out winning sequence into FASTA and output table.')
        WriteFasta(final_sequence,path.join(settings.output_folder,'ALLELES.FASTA'))
    
    output_table = pd.DataFrame(GenerateOutputTable(final_sequence, reads, query_seq)).sort_values(['qstart'])
    for x in ['qstart','qend']:
        output_table[x] = output_table[x].astype('Int64')
    output_table.to_csv(path.join(settings.output_folder,'Alleles.aln'),index=False, sep = '\t')

