"""
Runs DBSCAN and post-processing on randomly restarted ICAs.
The output files are S.csv and A.csv.

To execute the code:

mpiexec -n <n_cores> python run_ica.py -f FILENAME [-o OUT_DIR -i ITERATIONS -d DISTANCE -m MIN_FRAC]

n_cores: Number of processors to use
FILENAME: Name of original expression datafile
OUT_DIR: Path to output directory
ITERATIONS: Total number of ICA runs
DISTANCE: Max distance between points in a cluster for DBSCAN (optional, default 0.1)
MIN_FRAC: Minimum number of samples in a cluster (optional, default 0.5)
"""


from sklearn.cluster import DBSCAN
import numpy as np
import pandas as pd
from scipy import stats,sparse
from mpi4py import MPI
import time,sys,os,shutil,argparse

# Argument parsing
parser = argparse.ArgumentParser(description='Clusters independent components using DBSCAN')
parser.add_argument('-i',type=int,dest='iterations',required=True,
                    help='Number of ICA runs')
parser.add_argument('-d',type=float,default=0.1,dest='dist',
                    help='Maximum distance between points in a cluster for DBSCAN')
parser.add_argument('-m',type=float,default=0.5,dest='min_frac',
                    help='Minimum fraction of samples in a cluster for DBSCAN')
parser.add_argument('-o',dest='out_dir',default='',
                    help='Path to output file directory (default: current directory)')
args = parser.parse_args()


# -----------------------------------------------------------
# Split the work

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nWorkers = comm.Get_size()

# -----------------------------------------------------------

# Define parameters

n_iters = args.iterations

eps = args.dist # Clustering dissimilarity threshold for DBSCAN

# Minimum number of samples in a neighborhood to seed a cluster
min_samples = int(round(args.min_frac*n_iters))+1

# Set output files
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

# -----------------------------------------------------------
# Combine distance matrix

if rank == 0:
    print('\nCombining D matrix')

block = []
block_size = {}
for i in range(nWorkers):
    col = []
    for j in range(nWorkers):
        if i <= j:
            mat = sparse.load_npz(os.path.join(tmp_dir,'dist_{}_{}.npz'.format(i,j)))
            col.append(mat)
            block_size[i] = mat.shape[0]
            block_size[j] = mat.shape[1]
        else:
            mat = sparse.load_npz(os.path.join(tmp_dir,'dist_{}_{}.npz'.format(j,i)))
            col.append(mat.T)
    block.append(col)
D = sparse.bmat(block,'csr')

# ----------------------------------------------------------
# Convert to true distance matrix (originally dissimilarity)

data = D.data
indices = D.indices
indptr = D.indptr
D = sparse.csr_matrix((1-D.data,D.indices,D.indptr))


if rank == 0:
    t = timeit(t)
    print('\nClustering')

# Run DBSCAN

dbscan = DBSCAN(eps=eps,min_samples=min_samples,metric='precomputed')
labels = dbscan.fit_predict(D)
n_clusters = max(labels)+1

if rank == 0:
    t = timeit(t)
    print('\nIdentified',n_clusters,'clusters')
#-----------------------------------------------------------

# Place clustered components into correct bins

if rank == 0:
    print('Loading individual S and A matrices')
    start = 0
    end = 0
    S_bins = {i:[] for i in range(n_clusters)}
    A_bins = {i:[] for i in range(n_clusters)}

    for i in range(nWorkers):
        # Get labels for each partial matrix
        end += block_size[i]
        proc_labels = labels[start:end]
        start = end
    
        # Add parts of matrix to full table
        S_part = pd.read_csv(os.path.join(tmp_dir,'proc_{}_S.csv'.format(i)),index_col=0)
        A_part = pd.read_csv(os.path.join(tmp_dir,'proc_{}_A.csv'.format(i)),index_col=0)
        S_part.columns = range(S_part.shape[1])
        A_part.columns = range(A_part.shape[1])
        for i,lab in enumerate(proc_labels):
            if lab != -1:
                S_bins[lab].append(S_part[i])
                A_bins[lab].append(A_part[i])

# ----------------------------------------------------------

# Gather final S and A matrices

if rank == 0:

    t = timeit(t)
    print('\nGathering final S and A matrices')

    S_final = pd.DataFrame(columns=range(n_clusters),index=S_part.index)
    A_final = pd.DataFrame(columns=range(n_clusters),index=A_part.index)
    df_stats = pd.DataFrame(columns=['S_mean_std','A_mean_std','count'],index=range(n_clusters))

    for lab in range(n_clusters):
        S_clust = S_bins[lab]
        A_clust = A_bins[lab]
        
        # First item is base component
        Svec0 = S_clust[0]
        Avec0 = A_clust[0]

        # Make sure base component is facing positive
        if abs(min(Svec0)) > max(Svec0):
            Svec0 = -Svec0
            Avec0 = -Avec0

        S_single = [Svec0]
        A_single = [Avec0]

        # Add in rest of components
        for j in range(1,len(S_clust)):
            Svec = S_clust[j]
            Avec = A_clust[j]
            if stats.pearsonr(Svec,Svec0)[0] > 0:
                S_single.append(Svec)
                A_single.append(Avec)
            else:
                S_single.append(-Svec)
                A_single.append(-Avec)

        # Add centroid of cluster to final S matrix        
        S_final[lab] = np.array(S_single).T.mean(axis=1)
        A_final[lab] = np.array(A_single).T.mean(axis=1)

        # Get component stats
        df_stats.loc[lab,'S_mean_std'] = np.array(S_single).std(axis=1).mean()
        df_stats.loc[lab,'A_mean_std'] = np.array(A_single).std(axis=1).mean()
        df_stats.loc[lab,'count'] = len(S_single)
        
    print('\nFinal components created')
    t = timeit(t)


# Write to file
if rank == 0:
    # Remove components that exist in under 50% of runs
    good_comps = df_stats[df_stats['count'] > n_iters*.5].index
    
    S_final.columns = range(len(S_final.columns))
    A_final.columns = range(len(A_final.columns))
    df_stats = df_stats.reindex(good_comps).reset_index(drop=True)
    
    # Add L2 norm to stats
    df_stats['L2_norm'] = np.sqrt((A_final**2).sum(axis=1))

    print('Writing files to '+OUT_DIR)
    S_final.to_csv(os.path.join(OUT_DIR,'S.csv'))
    A_final.T.to_csv(os.path.join(OUT_DIR,'A.csv'))
    df_stats.to_csv(os.path.join(OUT_DIR,'component_stats.csv'))
    # Clean up tmp directories
    # shutil.rmtree(tmp_dir)
    
    t = timeit(t)
    print('Complete!')

