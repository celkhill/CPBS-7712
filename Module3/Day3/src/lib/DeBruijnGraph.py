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
        self.seq = seq
        self.next = np.array(next)
        self.prev = np.array(prev)
        self.visited = np.array(visited) #np.array(visited)
    
def append_to_node(node, sequence, direction, count = 1):
    if direction == 'next':
        node.next = np.append(node.next,[sequence]*count)
    elif direction == 'prev':
        node.prev = np.append(node.prev,[sequence]*count)

def RemoveEdgeFromGraph(node, edge_to_remove,direction):
    if direction == 'next':
        ind = np.where(node.next == edge_to_remove)[0]
        if len(ind)==0:
            raise ValueError('Edge %s not found in list! Please check.')
        else:
            node.next = np.delete(node.next,ind[-1])
        #add to our visited list
        node.visited = np.append(node.visited,edge_to_remove)
    elif direction == 'prev':
        ind = np.where(node.prev == edge_to_remove)[0]
        if len(ind)==0:
            raise ValueError('Edge %s not found in list! Please check.')
        else:
            node.prev = np.delete(node.prev,ind[-1])
        #add to our visited list
        node.visited = np.append(node.visited,'-%s'%edge_to_remove)

# def remove_from_next(node, sequence, count = 1):
#     if not sequence in node.next:
#         raise ValueError('Edge %s not found in list! Please check.'%sequence)

#     node.next[sequence] -= count
#     node.visited[sequence] += count
#     return node

# def remove_from_prev(node, sequence, count = 1):
#     if not sequence in node.prev:
#         raise ValueError('Edge %s not found in list! Please check.'%sequence)
#     node.prev[sequence] -= count
#     node.visited['-'+sequence] += count
#     return node

def reset_node(node):
    for sequence in node.visited:
        if sequence[0] == '-':
            append_to_node(node, sequence.replace('-',''), 'prev')
        else:
            append_to_node(node, sequence, 'next')
    node.visited = []

def get_options(node, mode, threshold = 0):
    if mode == 'next':
        un, counts = np.unique(node.next,return_counts= True)
    elif mode == 'prev':
        un, counts = np.unique(node.prev,return_counts= True)
    else:
        raise ValueError('Attribute not found in node! Please check input %s.'%mode)
    locs = np.where(counts > threshold)
    options = un[locs]
    return options

def get_node_degree(node, mode):
    if mode == 'next':
        degree = len(node.next)
        # degree = np.sum(node.next.values())
    elif mode == 'prev':
        degree = len(node.prev)
        # degree = np.sum(node.prev.values())
    else:
        raise ValueError('Attribute not found in node! Please check input %s.'%mode)
    return degree
