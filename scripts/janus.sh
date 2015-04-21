#!/bin/bash
# 
#to be used on the JANUS cluster at CU Boulder
#
# Author: prheenan (patrick.heenan@colorado.edu)
# -J: name of job
#SBATCH -J movie_generation
#SBATCH -A CLCSCI48300115
#--qos: specify which (janus is def)
#SBATCH --qos=janus
# -t specifies runtime hours:minudtes:seconds
#SBATCH -t 02:00:00
# -n: 4 cores (one per parallel process)
#SBATCH -n 4
# -N: number of cores
# -o: where to put output. %j: job ID
#SBATCH -o %j.out
# also can have stuff for CPUs

# add "~/bin" to the path
PATH=$PATH:$HOME/bin
bash ./scripts/_runScript.sh
