#!/bin/bash

function usage {
    printf "\nUsage: run_ica.sh [options] FILE\n"
    printf "\n"
    printf "OPTIONS\n"
    printf "  -i|--iter <n_iter>	Number of random restarts (default: 100)\n"
    printf "  -t|--tolerance <tol>      Tolerance (default: 1e-6)\n"
    printf "  -n|--n-cores              Number of cores to use (default: 8)\n"
    printf "  -o|--outdir <path>        Output directory for files (default: current directory)\n"
    printf "  -l|--logfile              Name of log file to use if verbose is off (default: ica.log)\n"
    printf "  -v|--verbose              Send output to stdout rather than writing to file\n"
    printf "  -h|--help                 Display help information\n"
    printf "\n"
    exit 1
}

# Handle arguments

OUTDIR=$(pwd)
TOL="1e-6"
ITER=100
CORES=8
VERBOSE=false
LOGFILE="ica.log"

POSITIONAL=()

while [[ $# -gt 0 ]]; do
    case $1 in
	-i|--iter)
            ITER=$2
	    shift; shift ;;
        -o|--out) 
            OUTDIR=$2
            shift; shift ;;
        -t|--tolerance)
            TOL=$2
            shift; shift ;;
        -n|--n-cores)
            CORES=$2
            shift; shift ;;
        -l|--logfile)
            LOGFILE=$2
            shift; shift ;;
        -v|--verbose)
            VERBOSE=true
            shift ;;
        -h|--help)
            usage ;;
        --) shift ; break ;;
        *) 
            POSITIONAL+=("$1")
            shift ;;
    esac
done

set -- "${POSITIONAL[@]}"

FILE="$1"

# Error checking

if [ "$FILE" = "" ]; then
    printf "Filename for expression data is required\n"
    usage
fi

if [ ! -f $FILE ]; then
    printf "ERROR: $FILE does not exist\n"
    exit 1
fi

# Run code

if [ "$VERBOSE" = true ]; then
    
    mpiexec -n $CORES python -u random_restart_ica.py -f $FILE -i $ITER -o $OUTDIR -t $TOL 2>&1
    mpiexec -n $CORES python -u compute_distance.py -i $ITER -o $OUTDIR 2>&1
    mpiexec -n $CORES python -u cluster_components.py -i $ITER -o $OUTDIR 2>&1
else
    echo "" > $LOGFILE
    mpiexec -n $CORES python -u random_restart_ica.py -f $FILE -i $ITER -o $OUTDIR -t $TOL >> $LOGFILE 2>&1
    mpiexec -n $CORES python -u compute_distance.py -i $ITER -o $OUTDIR >> $LOGFILE 2>&1
    mpiexec -n $CORES python -u cluster_components.py -i $ITER -o $OUTDIR >> $LOGFILE 2>&1
fi
