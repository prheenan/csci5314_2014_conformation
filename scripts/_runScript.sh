#!/bin/bash
# courtesy of: http://redsymbol.net/articles/unofficial-bash-strict-mode/
# (helps with debugging)
# set -e: immediately exit if we find a non zero
# set -u: undefined references cause errors
# set -o: single error causes full pipeline failure.
set -euo pipefail
IFS=$'\n\t'
# next lines courtesy of:
#  http://redsymbol.net/articles/unofficial-bash-striyesct-mode/
# (helps with debugging)
# set -e: immediately exit if we find a non zero
# set -u: undefined references cause errors
# set -o: single error causes full pipeline failure.
set -euo pipefail
IFS=$'\n\t'
reRunData=${1:-True}
defInput=${2:-"/lustre/janus_scratch/pahe3165/Data-csci7000/"}
defOutput=${3:-"/lustre/janus_scratch/pahe3165/Output-csci7000/"}

workingDir=`pwd`
# PRE:
# load the latest anaconda for python 2.0
module load python/anaconda-2.1.0 

# POST: all needed modules are loaded
cd $workingDir
echo $reRunData
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


