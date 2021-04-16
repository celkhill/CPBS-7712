
def ReadFasta(filename):
    """
    Function to read all lines of the fasta file format
    """
    with open(filename,'r') as f:
        data = f.readlines()
    return data

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
