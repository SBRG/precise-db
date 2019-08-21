"""
Performs ICA multiple times on a dataset to identify robust independent components.
In order to run this code, you must first install mpi4py (pip install mpi4py). 
The output files are stored in a temporary directory for the next processing step.

To execute the code:

mpiexec -n <n_cores> python random_restart_ica.py -f FILENAME -i ITERATIONS [-o OUT_DIR -t TOL]

n_cores: Number of processors to use
FILENAME: Path to log TPM data file
OUT_DIR: Path to output directory
ITERATIONS: Total number of ICA runs
TOL: Tolerance for ICA (optional, default: 1e-7)
"""


from sklearn.decomposition import FastICA,PCA
import numpy as np
import pandas as pd
from mpi4py import MPI
import time,sys,os,shutil,argparse

# Argument parsing
parser = argparse.ArgumentParser(description='Performs ICA with random initialization')
parser.add_argument('-f',dest='filename',required=True,
                    help='Path to expression data file')
parser.add_argument('-i',type=int,dest='iterations',required=True,
                    help='Number of ICA runs')
parser.add_argument('-t','--tol',type=float,default=1e-7,
                    help='ICA convergence tolerance (default: 1e-7)')
parser.add_argument('-o',dest='out_dir',default='',
                    help='Path to output file directory (default: current directory)')
args = parser.parse_args()


# -----------------------------------------------------------
# Split the work

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

n_iters = args.iterations

worker_tasks = {w:[] for w in range(size)}
w_idx = 0
for i in range(n_iters):
    worker_tasks[w_idx].append(i)
    w_idx = (w_idx + 1) % size

n_tasks = len(worker_tasks[rank])

# -----------------------------------------------------------
# Define parameters

x_file = os.path.abspath(args.filename)
tol =  args.tol# Tolerance for ICA. Larger values run faster,
           # but provide less accurate components.

pca_dims = 0.99 # Percent variance to keep for ICA

# Set output files
if args.out_dir == '':
    OUT_DIR = os.getcwd()
else:
    OUT_DIR = args.out_dir
    if rank == 0:
        if not os.path.isdir(OUT_DIR):
            os.makedirs(OUT_DIR)

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
    
# -----------------------------------------------------------
# Load Data
DF_data = pd.read_csv(x_file,index_col=0)
X = DF_data
n_genes,m_samples = X.shape

# Reduce dimensionality using PCA
pca = PCA().fit(X.transpose())
pca_var = np.cumsum(pca.explained_variance_ratio_)
k_comp = np.where(pca_var > pca_dims)[0][0] + 1
if rank == 0:
    print('Data: {} genes x {} samples'.format(n_genes,m_samples))
    print('Found {} dimensions from PCA'.format(k_comp))


# Create temporary directory for files
tmp_dir = os.path.join(OUT_DIR,'tmp')

if rank == 0:
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

# -----------------------------------------------------------
# Run ICA

if rank == 0:
    t = timeit(t)
    print('\nRunning ICA...')

S = []
A = []

t1 = time.time()

for counter,i in enumerate(worker_tasks[rank]):
    ica = FastICA(whiten=True,max_iter=int(1e10),tol=tol,n_components=k_comp)
    S.append(pd.DataFrame(ica.fit_transform(X),index=X.index))
    A.append(pd.DataFrame(ica.mixing_,index=X.columns))
    if rank == 0:
        print('Completed run {} of {} on Processor {}'.format(counter+1,n_tasks,rank))
        t = timeit(t)

S_all = pd.concat(S,axis=1)
S_all.columns = range(S_all.shape[1])
S_all.to_csv(os.path.join(tmp_dir,'proc_{}_S.csv'.format(rank)))
A_all = pd.concat(A,axis=1)
A_all.columns = range(A_all.shape[1])
A_all.to_csv(os.path.join(tmp_dir,'proc_{}_A.csv'.format(rank)))

# Wait for processors to finish
if rank == 0:
    test = 1
else:
    test = 0
test = comm.bcast(test,root=0)

if rank == 0:
    print('\nAll ICA runs complete!')
    timeit(t1)
