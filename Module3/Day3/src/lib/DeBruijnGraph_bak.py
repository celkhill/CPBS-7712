import numpy as np
import pdb
class GraphNode:
    """
    Class of a GraphNode in the de bruijn graph.
    Attributes:
        seq (str): String of the sequence relating to this node.
        next (list): Nodes with edges connecting to this node (forwards)
        prev (list): Nodes with edges connecting to this node (reverse)
        visited (list): A list of nodes visited during a solving iteration.
    """
    def __init__(self, seq='', next = [],  prev = [], visited = []):
        alphabet = ['G','C','T','A']
        self.seq = seq
        if next == []:
            next = np.array([(0,0,0,0)],dtype=[('A','int'),('C','int'),('G','int'),('T','int')])
        self.next = next
        # self.prev = np.array(prev)
        if prev == []:
            prev = np.array([(0,0,0,0)],dtype=[('A','int'),('C','int'),('G','int'),('T','int')])
        self.prev = prev
        if visited == []:
            visited = np.array([(0,0,0,0,0,0,0,0)],dtype=[('A','int'),('C','int'),('G','int'),('T','int'),('-A','int'),('-C','int'),('-G','int'),('-T','int')])
        self.visited = visited #np.array(visited)
    
def append_to_next(node,sequence,count = 1):
    """
    Method to append a sequence to the node.next attribute, appends = count times
    """
    # node.next = np.append(node.next,[sequence]*count)
    node.next[sequence] += count
    return node

def append_to_prev(node,sequence, count = 1):
    """
    Method to append a sequence to the node.prev attribute, appends = count times
    """
    # node.prev = np.append(node.prev,[sequence]*count)
    node.prev[sequence] += count
    return node

def remove_from_next(node, sequence, count = 1):
    if not sequence in node.next:
        raise ValueError('Edge %s not found in list! Please check.'%sequence)
    node.next[sequence] -= count
    node.visited[sequence] += count
    return node

def remove_from_prev(node, sequence, count = 1):
    if not sequence in node.prev:
        raise ValueError('Edge %s not found in list! Please check.'%sequence)
    node.prev[sequence] -= count
    node.visited['-'+sequence] += count
    return node

def reset_node(node):
    for sequence in node.visited:
        count = node.visited[sequence]
        if sequence[0] == '-':
            node.prev[sequence.replace('-','')] += count
        else:
            node.next[sequence] += count
        node.visited[sequence] == 0
    return node

def get_options(node, mode, threshold = 0):
    if mode == 'next':
        options = [x for x in node.next[1] if node.next[x] > threshold]
    elif mode == 'prev':
        options = [x for x in node.prev[1] if node.prev[x] > threshold]
    else:
        raise ValueError('Attribute not found in node! Please check input %s.'%mode)
    return options

def get_node_degree(node, mode):
    if mode == 'next':
        degree = np.sum([x for x in node.next[0]])
        # degree = np.sum(node.next.values())
    elif mode == 'prev':
        degree = np.sum([x for x in node.prev[0]])
        # degree = np.sum(node.prev.values())
    else:
        raise ValueError('Attribute not found in node! Please check input %s.'%mode)
    return degree
