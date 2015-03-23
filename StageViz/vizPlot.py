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
import matplotlib.pyplot  as plt

def plotSingle(X,Y,c1ByStages,c2ByStages,maxX,maxY,timeIdx):
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)
    xTmp = X[:,timeIdx]
    yTmp = Y[:,timeIdx]
    c1 = c1ByStages[:,timeIdx]
    c2 = c2ByStages[:,timeIdx]
    if (yTmp.nnz == 0 or xTmp.nnz == 0):
        return
    okayIdx = xTmp.nonzero()
    nTimes = len(okayIdx[0])
    ratio = np.reshape(c1[okayIdx]/c2[okayIdx],
                       (nTimes,1))
    xVals = xTmp[okayIdx]
    yVals = yTmp[okayIdx]
    cmap = plt.cm.brg
    nColors = 64
    colors = cmap(np.linspace(0,1,nColors))
    indices=  np.array([int(min(r*nColors,nColors-1)) for r in ratio] )
    objColor = colors[indices]
    plt.scatter(xVals,yVals,c=objColor)
    plt.xlim([0,maxX])
    plt.ylim([0,maxY])

def saveAsSubplot(XByStages,YByStages,c1ByStages,c2ByStages,outputDir):
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
            plotSingle(XByStages[i],YByStages[i],c1ByStages[i],c2ByStages[i],
                       maxX,maxY,timeIdx)
        pUtil.saveFigure(fig,outputDir + "t{:05d}".format(t),
                             overrideIO=True)

def vizIOSparse(inputFile):
    data = np.load(inputFile)
    return getSparseData(data)

def getAllStages(fileDict,conditionLabel,trialLabel,mWorking,condIdx,trialIdx):
    X = []
    Y = []
    fretC1 = []
    fretC2 = []
    for k,stageFile in enumerate(fileDict[conditionLabel][trialLabel]):
        # get the output file as a pickle (False, not numpy)
        mFile = vizFormatIOFile("IO_condition{:d}_trial{:d}_stage{:d}".
                                format(condIdx,trialIdx,k),mWorking,False)
        # get the x and y sparse matrices
        xSparse,ySparse,c1,c2 = \
                    cPoint.getCheckpoint(mFile,vizIOSparse,True,stageFile)
        X.append(xSparse)
        Y.append(ySparse)
        fretC1.append(c1)
        fretC2.append(c2)
    return X,Y,fretC1,fretC2
