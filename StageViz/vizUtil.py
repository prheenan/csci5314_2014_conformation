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

def vizFormatIOFile(fileName,outDir,numpy):
    base = outDir + fileName
    # either pickle or npz
    if (numpy):
        base += ".npz"
    else:
        base += ".pkl"
    return base

def getCheckpointFileDict(dataDir):
    # file dict has a key for each surface condition (e.g. 5x_PBS_25C_) and 
    # within that condition, a key for each trial (e.g. 5x_PBS_25C_200ms01)
    fileDict = defaultdict(defaultdict)
    checkDir = "/DataCheckpoints/"
    nStages = 3
    stageNames = [ "Stage_{:d}.npz".format(i+1) for i in range(nStages)]
    for conditionLabel in os.listdir(dataDir):
        fullPath = dataDir + conditionLabel + '/'
        trialDict = defaultdict(str)
        if (os.path.isdir(fullPath)):
            # found a directory in this trial
            for trialLabel in os.listdir(fullPath):
                pathToStages = fullPath + trialLabel +  checkDir
                files = []
                for fileName in stageNames:
                    fullFilePath = pathToStages + fileName
                    if (os.path.isfile(fullFilePath)):
                        files.append(fullFilePath)
                trialDict[trialLabel] = files
        fileDict[conditionLabel] = trialDict
    return fileDict

def getSparseData(data):
    times = data["times"]
    objIdx = data["idx"]
    flattened = np.concatenate(times)
    uniqueTimes = np.unique(flattened)
    frameTime = np.min(np.diff(uniqueTimes))
    nObjs = len(times)
    nTimes = len(uniqueTimes)
    lengths = [0]
    lengths.extend([len(t) for  t in times])
    indptr = np.cumsum(lengths)
    # get the zero-indexed time (e.g., 0.1 --> 0 if 0.1 is the frame time)
    # XXX fix min time bug...
    timeIdx = map( lambda x: int((x-frameTime)/frameTime), flattened)
    X = np.concatenate(data["X"])
    Y = np.concatenate(data["Y"])
    xSparse = csr_matrix((X,timeIdx,indptr),shape=(nObjs,nTimes))
    ySparse = csr_matrix((Y,timeIdx,indptr),shape=(nObjs,nTimes))
    # each row is an object, each columni is a time.
    # XXX add in others..
    normalize = lambda x: (x-np.min(x))/(np.max(x)-np.min(x))
    flatC1 = normalize(np.concatenate(data["FRET_C1"]))
    flatC2 = normalize(np.concatenate(data["FRET_C2"]))
    # get the CSR matrix
    fretC1 = csr_matrix((flatC1,timeIdx,indptr),shape=(nObjs,nTimes))
    fretC2 = csr_matrix((flatC2,timeIdx,indptr),shape=(nObjs,nTimes))
    return xSparse,ySparse,fretC1,fretC2

def getDirs(mOut,condLabel,trial,trialNum):
    trialLabel = condLabel + '_trial{:d}/'.format(trialNum)
    trialDir = mOut + trialLabel
    allStageDir = trialDir + "allStages/"
    gUtil.ensureDirExists(allStageDir)
    return trialLabel,allStageDir
