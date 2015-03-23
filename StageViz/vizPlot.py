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
from vizUtil import vizFormatIOFile,getCheckpointFileDict,getSparseData,getDirs

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
    nTimes = [ X.shape[1]-1 for X in XByStages]
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
        pUtil.saveFigure(fig,outputDir + "t{:05d}".format(t),
                             overrideIO=True)

def vizIOSparse(inputFile):
    data = np.load(inputFile)
    return getSparseData(data)

def getAllStages(fileDict,conditionLabel,trialLabel,mWorking,condIdx,trialIdx):
    X = []
    Y = []
    for k,stageFile in enumerate(fileDict[conditionLabel][trialLabel]):
        # get the output file as a pickle (False, not numpy)
        mFile = vizFormatIOFile("IO_condition{:d}_trial{:d}_stage{:d}".
                                format(condIdx,trialIdx,k),mWorking,False)
        # get the x and y sparse matrices
        xSparse,ySparse = \
                    cPoint.getCheckpoint(mFile,vizIOSparse,False,stageFile)
        X.append(xSparse)
        Y.append(ySparse)
    return X,Y
