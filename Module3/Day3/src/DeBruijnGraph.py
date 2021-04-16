import numpy as np
class GraphNode:
    """
    Class of a GraphNode in the de bruijn graph.
    Attributes are: seq, next, prev, visited
    """
    def __init__(self, seq='', next = [],  prev = [], visited = []):
        self.seq = seq
        self.next = np.array(next)
        self.prev = np.array(prev)
        self.visited = np.array(visited)
    
    def append_to_next(self,letter,count = 1):
        self.next = np.append(self.next,[letter]*count)

    def append_to_prev(self,letter, count = 1):
        self.prev = np.append(self.prev,[letter]*count)
