import os
from collections import defaultdict
import numpy as np
# need to import the sparse matrices and stacking funciton
from scipy.sparse import csr_matrix,csc_matrix
import sys
sys.path.append('../Utilities')
import CheckpointUtilities as cPoint
import GenUtilities as gUtil
import PlotUtilities as pUtil
import matplotlib.pyplot as plt
from vizUtil import getCheckpointFileDict,getDirs
from vizPlot import saveAsSubplot,getAllStages
from multiprocessing import Process # use parallel processing 

import subprocess
import argparse

def parseCmdLine():
    # get the input/output directories
    parser = argparse.ArgumentParser(description='Protein Visualizatio args')
    parser.add_argument('--inPath', type=str, default="../output/",
                        help="Folder where formatted .dat file reside")
    parser.add_argument('--outPath', type=str, default="./vizOut/",
                        help="Default output directory (base)")
    parser.add_argument('--cachePath', type=str, default="./vizTmp/",
                        help="Default cache output directory (base)")
    args = parser.parse_args()
    dataDir = args.inPath
    mWorking = args.cachePath
    mOut = args.outPath
    return dataDir,mWorking,mOut

def generateMovie(args):
    command = 'ffmpeg'
    return os.system(command + ' ' + args)

def saveSingleTrial(mWorking,mOut,condition,condNum,trial,trialNum,fps=10,
                    vizFileFormat="t{:05d}",ffmpegFormat = "t%05d",
                    vizExt=".png"):
        # for this condition and trial, get the directories and all the 
    # data, then save
    trialDir,allStageDir =  getDirs(mOut,condition,trial,trialNum)
    X,Y,c1,c2 = getAllStages(fileDict,condition,trial,mWorking,condNum,trialNum)
    # format the ffmpeg arguments as we want them
    argStr = ('-i {:s} -c:v libx264 -r {:d} -y -pix_fmt yuv420p '+
              '{:s}1_movie_{:s}_trial_{:d}.mov').\
        format(allStageDir+ffmpegFormat+vizExt,fps,allStageDir,
               condition,trialNum)
    saveAsSubplot(X,Y,c1,c2,allStageDir,vizFileFormat)
    # POST: all videos saved for this trial. make the movie
    generateMovie(argStr)

def saveConditions(condition,condNum,conditionKeys,workDir,outDir):
    for j,trial in enumerate(conditionKeys):
        saveSingleTrial(workDir,outDir,condition,condNum,trial,j)

if __name__ == '__main__':
    inDir,workDir,outDir = parseCmdLine()
    # next two must match, for the automatic video encoding to work
    gUtil.ensureDirExists(outDir)
    gUtil.ensureDirExists(workDir)
    # get all the files. returns a dictionary of dictionaries.
    # each key in the outer (first) dictionary is a condition
    # each key in the innter (second) dictionary is a trial for that condition
    fileDict = getCheckpointFileDict(inDir)
    # loop through each condition and trial
    conditionArr = fileDict.keys()
    processes= []
    for i,condition in enumerate(conditionArr):
        print("Forking off a process for condition {:s}".format(condition))
        func = saveConditions
        funcArgs = (condition,i,fileDict[condition].keys(),workDir,outDir)
        p = (Process(target=func, args=funcArgs))
        processes.append(p)
        p.start()
    # wait until all processes are done...
    for p in processes:
        p.join()

