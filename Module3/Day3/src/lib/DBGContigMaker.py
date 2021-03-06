import numpy as np
import pandas as pd
from os import path
import re
import pdb
from .FileIOUtils import ReadFasta, LoadTableFromFile
from .CreateDBG import GetReverseComplement
from .DeBruijnGraph import *


def ConstructSequenceFromNodes(node_list):
    """
    Function to reconstruct a sequence from a list of node objects. Returns a string of the nucleotide sequence.
    Inputs:
        node_list (list): A list of GraphNode objects
    Output:
        out_sequence (str): A DNA sequence reconstructed from the node list.
    """

    out_sequence = ''.join([n.seq[0] for n in node_list[0:-1]])+node_list[-1].seq
    return out_sequence

#def RemoveEdgeFromGraph(node,edge_to_remove, direction = 'next'):
#    """
#    Removes edges from a list and returns the node object.
#    Inputs:
#        node object (GraphNode): A GraphNode object.
#        edge_to_remove (str): The edge that will be removed.
#    Output:
#        node (GraphNode): node object with the edge_to_remove removed, type GraphNode.
#    """

#    remove_edge(node, edge_to_remove,direction)
#    if direction == 'next':
#        ind = np.where(node.next == edge_to_remove)[0]
#        if len(ind)==0:
#            raise ValueError('Edge %s not found in list! Please check.')
#        else:
#            node.next = np.delete(node.next,ind[-1])
#        #add to our visited list
#        node = append_to_visited(node,np.append(node.visited,edge_to_remove))
#    elif direction == 'prev':
#        ind = np.where(node.prev == edge_to_remove)[0]
#        if len(ind)==0:
#            raise ValueError('Edge %s not found in list! Please check.')
#        else:
#            node.prev = np.delete(node.prev,ind[-1])
#        #add to our visited list
#        node.visited = np.append(node.visited,'-%s'%edge_to_remove)

#    return node

def SelectNextNode(current_node,node_table, mode = 'forwards'):
    """
    Function to select the next node when generating a contig.
    Input:
        current_node (GraphNode): A GraphNode object of the current node we are on.
        node_table (dict): The dictionary of all nodes in our graph.
        mode (str):   'forwards' extends forwards in the graph using node.next
                'reverse' extends backwards in the graph using node.prev
    Output:
        next_na (str, False): The string value of the next nucleotide in the sequence. Will give value of False if no next nucleotide exists.
    """
    if mode == 'forwards':
        options = get_options(current_node,'next')
    elif mode == 'reverse':
        options = get_options(current_node,'prev')
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
    return next_na


# solve for our query sequence
def FindQueryInGraph(query_seq,kmersize,node_table):
    """
    Function to find a given query sequence in our graph. Will return an error if the query is not found.
    Inputs:
        query_seq (str): The sequence to be queried as type string.
        kmersize (int): The kmersize used in the graph.
        node_table (dict): The dictionary with the nodes and edges of the graph.
    Outputs:
        solved_nodes (list): A list of the node objects that reconstruct the query_seq input.
    """
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
                raise ValueError('No solution found! Stopping at this sequence: %s'%ConstructSequenceFromNodes(solved_nodes))
            else:
                RemoveEdgeAndReverseComplement(node_table, next_na, cn)
                
        solved_nodes.append(cn)    
        idx+=1
    return solved_nodes

def RemoveEdgeAndReverseComplement(node_table, next_na, node, direction = 'next'):
    """
    Function to remove a given edge and its reverse complement edge from the graph.
    Inputs:
        node_table (dict): A De Bruijn Graph as constructed by CreateDBGraph
        next_na (str): The edge to remove from the current node.
        node (GraphNode): The node to remove an edge from.
    Outputs:
        None, edits the node objects passed.
    """
    #remove the edge!
    if direction == 'next':
        RemoveEdgeFromGraph(node,next_na, 'next') 
        rev_node = node_table[GetReverseComplement(node.seq[1:]+next_na)]
        #only have to remove the reverse complement if they aren't identical
        if rev_node != node:
            #remove the reverse complement!
            try:
                RemoveEdgeFromGraph(rev_node,GetReverseComplement(node.seq[0]), 'next')
            except:
                pdb.set_trace()
        #remove it from the previous list of the next node as well!
        next_node = node_table[node.seq[1:]+next_na]
        #a little bit of recursion here
        RemoveEdgeAndReverseComplement(node_table, node.seq[0], next_node, direction = 'prev') 

    if direction == 'prev':
        RemoveEdgeFromGraph(node, next_na, 'prev')
        rev_node = node_table[GetReverseComplement(next_na+node.seq[:-1])]
        #only have to remove the reverse complement if they aren't identical
        if rev_node != node:
            #remove the reverse complement!
            try:
                RemoveEdgeFromGraph(rev_node,GetReverseComplement(node.seq[-1]),'prev')
            except:
                pdb.set_trace()
        #we don't need to remove from the prev_node.next because we won't be going in that direction any more
    return

#extend our solution!
def CreateContigFromGraph(solved_nodes,node_table,direction = 'forwards'):
    """
    Creates a contig from the graph. Default direction is forwards.
    Inputs:
        solved_nodes (list): A list of GraphNode objects that seed our contig. If the list is empty a random node will be chosen.
        node_table (dict): The dictionary with the nodes and edges of the graph.
        direction (str):   'forwards' extends forwards in the graph using node.next
                'reverse' extends backwards in the graph using node.prev
    Outputs:
        solved_nodes (list): An extended version of solved_nodes, a list of GraphNode objects that can be reconstructed into a contig.
    """

    if len(solved_nodes)>0:
        cn = solved_nodes[-1]
    else:#pick a random starting node TODO: make this better
        cn = node_table[np.random.choice(list(node_table.keys()))]
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
                    RemoveEdgeAndReverseComplement(node_table, next_na, cn)
                elif direction == 'reverse':
                    next_node = node_table[next_na+cn.seq[:-1]]
                    RemoveEdgeAndReverseComplement(node_table, next_na, cn,direction = 'prev')
                solved_nodes.append(next_node)
                cn = next_node
            except KeyError:
                if direction == 'forwards':
                    seq = cn.seq[1:]+next_na
                elif direction == 'reverse':
                    seq = next_na + cn.seq[:-1]
                print('Sequence %s not found in the node table!'%seq)
                keep_looping = False
            except ValueError:
                pdb.set_trace() 
    return solved_nodes


def GenerateOutputTable(final_sequence, reads, query):
    """
    Function to generate the output table as a pandas dataframe.
    Inputs:
        final_sequence (str): A string of the contig
        reads (dict): Key, value pairs of the reads as generated by ParseReads
        query_seq (str): A string of the query sequence.
    Outputs: 
        output_table (pd.DataFrame): A pandas dataframe type of the output table.
    """

    #can pass either a filename or something already read in 
    if type(reads) == str:
        if path.isfile(reads):
            reads = ReadFasta(reads)
        else:
            raise ValueError('Reads FASTA file not found! Please check argument %s'%reads)

    #can pass either a filename or something already read in 
    if type(query) == str:
        if path.isfile(query):
            query_seq = list(ReadFasta(query).values())[0]
        else:
            raise ValueError('Query FASTA file not found! Please check argument %s'%query)
    else:
        query_seq = query

    output_table = []
    for read in reads:
        for direction in ['forwards','reverse']:
            if direction == 'reverse': search_term = GetReverseComplement(reads[read])
            elif direction == 'forwards': search_term = reads[read]
            for match in re.finditer(search_term,final_sequence):
                query_match = re.search(search_term,query_seq)
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
    output_table = pd.DataFrame(output_table).sort_values(['qstart'])
    #there will be Nans in this column, so set it to Int64 to allow us to convert to int format with nans
    for x in ['qstart','qend']:
        output_table[x] = output_table[x].astype('Int64')

    return output_table

def ResetNodeTable(node_table):
    """
    Function to reset the node table pack after a contig has been found and edges are removed.
    Inputs:
        node_table (dict): A dictionary containing our nodes and edges, with removed edges in node.visited
    Outputs:
        node_table (dict): A reset node table with node.visited empty for each node.
    """

    for key in node_table:
        reset_node(node_table[key])
        # for edge in node_table[key].visited:
        #     if '-' in edge:
        #         node_table[key].prev = np.append(node_table[key].prev,edge.replace('-',''))
        #     else:
        #         node_table[key].next = np.append(node_table[key].next,edge)
        # node_table[key].visited = np.array([])
    return node_table

def SolveGraphWithQuery(query_fasta, node_table, max_iter, kmersize):
    """
    Function to find a solution within the graph based on a query sequence.
    Inputs:
        query_fasta (str): A FASTA file (if its a valid file) or string sequence to query in the graph.
        node_table (dict): A DBG created from CreateDBGraph.
        max_iter (int): Int type, iterations to perform to find the longest solution.
        kmersize (int): The kmersize of the graph.
    Ouputs:
        solved_seq (str): String type of the longest sequence found.
    """
    #can pass either a filename or something already read in 
    if type(query_fasta) == str:
        if path.isfile(query_fasta):
            query_seq = ReadFasta(query_fasta)
            query_seq = list(query_seq.values())[0]
        else:
            query_seq = query_fasta
    else:
        raise ValueError('Query FASTA file not found! Please check argument %s'%query_fasta)

    #can pass either a filename or something already read in 
    if type(node_table) == str:
        if path.isfile(node_table):
            print('Loading existing DBG from file...')
            node_table = LoadTableFromFile(node_table,parse_json=node_table.endswith('.json'))
        else:
            raise ValueError('Node table file not found! Please check argument %s'%node_table)

    best_solved_nodes = ['One']
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

    solved_seq = ConstructSequenceFromNodes(best_solved_nodes)
    return solved_seq
