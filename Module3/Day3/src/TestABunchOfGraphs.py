from main import ConstructArgs, RunAllAnalysis
import timeit
import pandas as pd


kmersizes = [14, 18, 22, 26,30]
# kmersizes = [22]
execution_times = []
for kmersize in kmersizes:
    start = timeit.default_timer()
    args = ['-reads_fasta', '../data/READS.fasta',
            '-query_fasta','../data/QUERY.fasta',
            '-output_folder', '../analysis/kmer_%d'%kmersize,
            # '-existing_graph_filename','../analysis/kmer_%d/DeBruijnGraph_kmer%d.json'%(kmersize,kmersize),
            '--save_kmer_table',#'../analysis/kmer_%d/kmertable.json'%kmersize,
            '-kmersize',str(kmersize),
            '-solving_iterations',str(50),
            '--generate_histograms']
    settings= ConstructArgs(args)
    RunAllAnalysis(settings)
    end = timeit.default_timer()
    execution_times.append((end-start)/60)
    print('Done!!! Execution time was %s'%execution_times[-1])

pd.DataFrame(zip(kmersizes,execution_times),columns=['Kmersize','Execution time (min)']).to_csv('./execution_times.csv')
