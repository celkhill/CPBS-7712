from DBGContigMaker import RunAllAnalysis, ConstructArgs
import pdb
import timeit
import pandas as pd

kmersizes = [13, 17, 21, 27]
execution_times = []
for kmersize in kmersizes:
    start = timeit.default_timer()
    args = ['-reads_fasta', '../data/READS.fasta',
            '-query_fasta','../data/QUERY.fasta',
            '-output_folder', '../analysis/kmer_%d'%kmersize,
            '-kmersize',str(kmersize),
            '-solving_iterations',str(50),
            '--generate_histograms']
    settings= ConstructArgs(args)
    RunAllAnalysis(settings)
    end = timeit.default_timer()
    execution_times.append(end-start)
    print('Done!!! Execution time was %s'%execution_times[-1])

pd.DataFrame(zip(kmersizes,execution_times),columns=['Kmersize','Execution time']).to_csv('./execution_times.csv')
