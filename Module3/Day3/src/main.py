import argparse
from os import path
from FileIOUtils import LoadTableFromFile, SaveOutputTable, WriteFasta, ReadFasta
from CreateDBG import CreateDeBruijnGraphAndSave, GenerateHistograms
from DBGContigMaker import SolveGraphWithQuery, GenerateOutputTable

def ConstructArgs(argv):
    """
    Function to construct the arguments for the CLI of this program.
    Inputs:
        argv: Arguments passed via CLI or through a python script. Type is list.
    Outputs:
        args: The parsed arguments.

    """
    defaults = {'output_folder':'../analysis/',
            'output_filetype': 'json',
            'kmersize':27,
            'generate_histograms':False,
            'existing_graph_filename':None,
            'solving_iterations':10,
            'existing_allele': None,
            'existing_kmer_table':None}

    parser = argparse.ArgumentParser(description='Here are the arguments for the query sequence finder tool:', formatter_class = argparse.ArgumentDefaultsHelpFormatter)

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
    opt_group.add_argument('-existing_graph_filename',action = 'store',dest = 'existing_graph_filename',
            help='Argument to use an existing graph if one has already been created.',default= defaults['existing_graph_filename'])
    opt_group.add_argument('-existing_kmer_table',action = 'store',dest = 'existing_kmer_table',
            help='Argument to use an existing kmer table if one has already been created.',default= defaults['existing_kmer_table'])
    opt_group.add_argument('-solving_iterations',action = 'store',dest='solving_iterations',
            type = int, help='Number of solving iterations to try before selecting the longest contig found.', default = defaults['solving_iterations'])
    opt_group.add_argument('-existing_allele', action = 'store',dest = 'existing_allele',
            help='Argument to use an existing solved allele that has already been created.', default = defaults['existing_allele'])
    opt_group.add_argument('--save_kmer_table',action = 'store_true', dest = 'save_kmer_table',
            help = 'Option to save the kmer table generated from the reads after error correction (if performed). Simple pass --save_kmer_table to set to True.')
    opt_group.add_argument('--do_error_correction',action = 'store_true',dest = 'do_error_correction',
            help = 'Option to perform error correction when generating the De Bruijn Graph. Simply pass --do_error_correction argument to set to True.',default=True)
    opt_group.add_argument('--generate_histograms',action='store_true', dest = 'generate_histograms',
            help = 'Option to generate histograms of the read length and kmer edge degree. Simply pass --generate_histograms argument to set to True.')

    return parser.parse_args(argv)

def RunAllAnalysis(settings):
    
    de_bruijn_graph_filename = path.join(settings.output_folder,'DeBruijnGraph_kmer%d.'%settings.kmersize+ settings.output_filetype)

    #######################################

    ### Generating the Graph ###

    #if we just want to use the existing allele, do it right now!
    if settings.existing_allele is not None:
        print('Solved allele file found! Generating output table...')
        final_sequence = list(ReadFasta(settings.existing_allele).values())[0]
    else:
        #load an existing graph if we have one, otherwise make a new one!
        if settings.existing_graph_filename is None:
            node_table = CreateDeBruijnGraphAndSave(settings.reads_fasta,
                    settings.kmersize,
                    de_bruijn_graph_filename,
                    settings.existing_kmer_table,
                    settings.do_error_correction,
                    settings.generate_histograms,
                    settings.save_kmer_table)
        else:#pass the loading on to the "SolveGraphWithQuery" function
            node_table = settings.existing_graph_filename
    
        #now let's solve!
        final_sequence = SolveGraphWithQuery(settings.query_fasta, node_table, settings.solving_iterations, settings.kmersize)

        print('\n\nWinning sequence length: %d \nWinning sequence: %s'%(len(final_sequence),final_sequence))
        print('\n\nWriting out winning sequence into FASTA and output table.')
        WriteFasta(final_sequence,path.join(settings.output_folder,'ALLELES.FASTA'))
    SaveOutputTable(GenerateOutputTable(final_sequence, settings.reads_fasta, settings.query_fasta), path.join(settings.output_folder,'Alleles.aln'))
    print('Program finished!')
if __name__ == "__main__":

    ### Parsing arguments and settings ###
    #create the arg parser
    settings = ConstructArgs(sys.argv[1:])

    RunAllAnalysis(settings)
