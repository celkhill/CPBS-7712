import pandas as pd
from os import path, makedirs
try:
    import ujson as json
except:
    import json

from DeBruijnGraph import GraphNode

def ReadFasta(filename):
    """
    Function to read all lines of the fasta file format
    """
    with open(filename,'r') as f:
        data = f.readlines()
    return ParseReads(data)

def ParseReads(data):
    """
    Function to parse the reads from a fasta file into a dictionary with key: value as read ID: sequence
    """
    reads = {}
    for num in range(len(data)):
        if data[num].startswith('>'):
            key = data[num].rstrip() # remove newline char
            #handle multiline FASTA
            idx = 1
            val = []
            while not data[num+idx].startswith('>'):
                val += data[num+idx].rstrip() # remove newline char
                idx+=1
                #break out if we hit the end
                if num+idx == len(data):
                    break

            reads[key] = ''.join(val)
    return reads

def WriteFasta(sequence, filename):
    """
    Function to write the sequence out to a FASTA file.
    Inputs:
        sequence: The nucleotide sequence as type string.
        filename: The filename to save the sequence out to.
    Outputs:
        None.
    """
    with open(filename,'w') as f:
        f.write('>FinalContig:Length:%s'%len(sequence))
        f.write('\n')
        f.write(sequence)
        f.close()
    return

def ParseJson(jdata):
    """
    Function to parse the JSON back into a node table or a list of node objects.
    JSON files save everything as strings, so we need to read our dictionary of strings back into GraphNode namespace.
    Inputs:
        jdata: The data read from a JSON file.
    Outputs:
        The node_table dictionary or node list with GraphNode objects instead of dicts of strings.
    """

    if type(jdata) == dict:
        #delete the kmersize so we can just get the node table
        return {x: GraphNode(**jdata[x]) for x in jdata}
    elif type(jdata)==list:
        return [GraphNode(**x) for x in jdata]
    else:
        raise ValueError('Input type not recognized. Please pass either a dictionary or a list to parse!')

def SaveTableToFile(table, outname):
    if not path.exists(path.dirname(outname)):
        makedirs(path.dirname(outname))
    print('Writing out DBG to file %s'%outname)
    if outname.endswith('.json'):
        with open(outname, 'w') as fp:
            json.dump(table,fp)
    elif outname.endswith('.dat'):
        with open(outname, 'wb') as fp:
            pickle.dump(table,fp)

def LoadTableFromFile(filepath, parse_json = False):
    if not path.exists(filepath):
        raise ValueError('Cannot load file! Filepath not found '+ filepath)
    else:
        if filepath.endswith('.json'):
            with open(filepath, 'r') as fp:
                data = json.load(fp)
            if parse_json:
                data = ParseJson(data)
        elif filepath.endswith('.dat'):
            with open(filepath, 'rb') as fp:
                data = pickle.load(fp)
    return data

def SaveOutputTable(output_table,filename):
    if not type(output_table) == pd.DataFrame:
        raise TypeError('Output table is not a pandas dataframe! Please output table.')
    if not path.isdir(path.dirname(filename)):
        makedirs(path.dirname(filename))
    output_table.to_csv(filename, index=False)

