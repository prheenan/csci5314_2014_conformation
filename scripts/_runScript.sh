#!/bin/bash
# Author: prheenan (patrick.heenan@colorado.edu)
# -J: name of job
#SBATCH -J movie_generation
#SBATCH -A CLCSCI48300115
#--qos: specify which (janus is def)
#SBATCH --qos=janus
# -t specifies runtime hours:minutes:seconds
#SBATCH -t 12:00:00
# -n: 5 cores (one per parallel process)
#SBATCH -n 5
# memory per CPU, in MB
#SBATCH --mem-per-cpu 1000
# -N: number of cores
#SBATCH -N 5
# -o: where to put output. %j: job ID
#SBATCH -o %j.out
# also can have stuff for CPUs

set -euo pipefail
IFS=$'\n\t'
reRunData=${1:-False}
defInput=${2:-"/lustre/janus_scratch/pahe3165/Data-csci7000/"}
defOutput=${3:-"/lustre/janus_scratch/pahe3165/Output-csci7000/"}
workingDir=`pwd`
# PRE:
# load the latest anaconda for python 2.0

# POST: all needed modules are loaded
cd $workingDir
if [ "$reRunData" = True ]; then
    cd ./Stage0_Run
    python Main.py --inPath ${defInput} --outPath ${defOutput}
fi
# POST: we have should data, assuming the user didn't incorrectly specify reRun
cd ${workingDir}
cd ./StageViz
# visualization uses the *output* from the analysis
args="--inPath $defOutput --outPath ${defOutput}viz/vizOut/ --cachePath ${defOutput}viz/vizTmp/"
eval "python vizMain.py $args"


