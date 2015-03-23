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
    times = data['times']
    objIdx = data['idx']
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
    X = np.concatenate(data['X'])
    Y = np.concatenate(data['Y'])
    xSparse = csr_matrix((X,timeIdx,indptr),shape=(nObjs,nTimes))
    ySparse = csr_matrix((Y,timeIdx,indptr),shape=(nObjs,nTimes))
    # each row is an object, each column is a time.
    # XXX add in others..
    #fretC1 = data['FRET_C1']
    #fretC2 = data['FRET_C2']
    #velY = data['velY']
    #velX = data['velX']
    return xSparse,ySparse

def vizFormatIOFile(fileName,outDir,numpy):
    base = outDir + fileName
    # either pickle or npz
    if (numpy):
        base += ".npz"
    else:
        base += ".pkl"
    return base

def plotSingle(X,Y,maxX,maxY,timeIdx):
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)
    xTmp = X[:,timeIdx]
    yTmp = Y[:,timeIdx]
    if (yTmp.nnz == 0 or xTmp.nnz == 0):
        return
    xVals = xTmp[xTmp.nonzero()]
    yVals = yTmp[yTmp.nonzero()]
    plt.plot(xVals,yVals,'ro')
    plt.xlim([0,maxX])
    plt.ylim([0,maxY])

def saveImagesSingle(X,Y,outputDir,stageNum):
    nTimes = X.shape[1] # columns give times
    maxX = X.max()
    maxY = Y.max()
    for i in range(nTimes):
        fig = plt.figure()
        # Stack the two rows together
        plotSingle(X,Y,maxX,maxY,i)
        pUtil.saveFigure(fig,outputDir + "t{:05d}".format(i),
                         overrideIO=True)

def saveAsSubplot(XByStages,YByStages,outputDir):
    maxX = 0
    maxY = 0
    # number of times given by columns
    nTimes = [ X.shape[1] for X in XByStages]
    maxTimes = np.max(nTimes)
    numStages = len(XByStages)
    for i in range(numStages):
        maxX = max(XByStages[i].max(),maxX)
        maxY = max(YByStages[i].max(),maxY)
    for t in range(maxTimes):
        fig = plt.figure()
        for i in range(numStages):
            plt.subplot(1,numStages,i+1)
            timeIdx = min(t,nTimes[i])
            plotSingle(XByStages[i],YByStages[i],maxX,maxY,timeIdx)
        pUtil.saveFigure(fig,outputDir + "t{:05d}".format(timeIdx),
                             overrideIO=True)

def vizIOSparse(inputFile):
    data = np.load(inputFile)
    return getSparseData(data)

dataDir = "../output/"
mWorking = "./vizTmp/"
mOut = "./vizOut/"
gUtil.ensureDirExists(mOut)
gUtil.ensureDirExists(mWorking)
fileDict = getCheckpointFileDict(dataDir)

for i,key in enumerate(fileDict.keys()):
    for j,keyV in enumerate(fileDict[key].keys()):
        trialLabel = key + '_trial{:d}/'.format(j)
        trialDir = mOut + trialLabel
        allStageDir = trialDir + "allStages/"
        gUtil.ensureDirExists(trialDir)
        gUtil.ensureDirExists(allStageDir)
        X = []
        Y = []
        for k,stageFile in enumerate(fileDict[key][keyV]):
            stageDir = trialDir + "Stage{:d}/".format(k)
            gUtil.ensureDirExists(stageDir)
            # get the output file as a pickle (False, not numpy)
            mFile = vizFormatIOFile("IO_condition{:d}_trial{:d}_stage{:d}".
                                    format(i,j,k),mWorking,False)
            # get the x and y sparse matrices
            xSparse,ySparse = \
                cPoint.getCheckpoint(mFile,vizIOSparse,False,stageFile)
            X.append(xSparse)
            Y.append(ySparse)
            #saveImagesSingle(xSparse,ySparse,stageDir,k)
        saveAsSubplot(X,Y,allStageDir)

