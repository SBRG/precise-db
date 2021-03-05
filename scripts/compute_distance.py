"""
Clusters the S vectors generated from random_restart_ica.py
The output files are S.csv and A.csv.

To execute the code:

mpiexec -n <n_cores> python cluster_components.py -i ITERATIONS [-o OUT_DIR ]

n_cores: Number of processors to use
OUT_DIR: Path to output directory
ITERATIONS: Total number of ICA runs
"""


import numpy as np
import pandas as pd
from scipy import sparse
from mpi4py import MPI
import time,sys,os,shutil,argparse

# Argument parsing
parser = argparse.ArgumentParser(description='Generates Distance Matrix')
parser.add_argument('-i',type=int,dest='iterations',required=True,
                    help='Number of ICA runs')
parser.add_argument('-o',dest='out_dir',default='',
                    help='Path to output file directory (default: current directory)')
args = parser.parse_args()


# -----------------------------------------------------------
# Split the work

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nWorkers = comm.Get_size()

n_iters = args.iterations
# -----------------------------------------------------------

#-----------------------------------------------------------
# Parse directories

if args.out_dir == '':
    OUT_DIR = os.getcwd()
else:
    OUT_DIR = args.out_dir

tmp_dir = os.path.join(OUT_DIR,'tmp')

#-----------------------------------------------------------

def timeit(start):
    end = time.time()
    t = end-start
    if t < 60:
        print('{:.2f} seconds elapsed'.format(t))
    elif t < 3600:
        print('{:.2f} minutes elapsed'.format(t/60))
    else:
        print('{:.2f} hours elapsed'.format(t/3600))
    return end

t = time.time()

# ----------------------------------------------------------

# Define tasks as blocks for outer products
tasks = [(i,j) for i in range(nWorkers) for j in range(i,nWorkers)]

# Split up tasks evenly

worker_tasks = {w:[] for w in range(nWorkers)}
w_idx = 0
for task in tasks:
    worker_tasks[w_idx].append(task)
    w_idx = (w_idx + 1) % nWorkers

n_tasks = len(worker_tasks[rank])

if rank == 0:
    print('\nComputing clusters...')

t1 = time.time()

counter = 0
for i,j in worker_tasks[rank]:
    S1 = pd.read_csv(os.path.join(tmp_dir,'proc_{}_S.csv'.format(i)),index_col=0)
    S2 = pd.read_csv(os.path.join(tmp_dir,'proc_{}_S.csv'.format(j)),index_col=0)
    dist = abs(np.dot(S1.T,S2))
    dist[dist < .5] = 0
    sparse.save_npz(os.path.join(tmp_dir,'dist_{}_{}.npz'.format(i,j)),sparse.coo_matrix(np.clip(dist,0,1)))
    counter += 1
    #print('Completed distance block {} of {} on Processor {}'.format(counter,n_tasks,rank)
    #t = timeit(t)

# Wait for processors to finish
if rank == 0:
    test = 1
else:
    test = 0
test = comm.bcast(test,root=0)
if rank == 0:
    print('\nDistance matrix completed!')
    timeit(t1)

