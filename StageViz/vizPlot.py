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
import matplotlib.cm as cmaps

import matplotlib.pyplot as plt
from vizUtil import vizFormatIOFile,getCheckpointFileDict,getSparseData,getDirs
import matplotlib

def hideAxis():
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)

def plotSingle(X,Y,c1ByStages,c2ByStages,maxX,maxY,timeIdx,cmap,nColors):
    hideAxis()
    xTmp = X[:,timeIdx]
    yTmp = Y[:,timeIdx]
    c1 = c1ByStages[:,timeIdx]
    c2 = c2ByStages[:,timeIdx]
    if (yTmp.nnz == 0 or xTmp.nnz == 0):
        return
    okayIdx = xTmp.nonzero()
    nTimes = len(okayIdx[0])
    ratio = np.reshape(c2[okayIdx]/c1[okayIdx],
                       (nTimes,1))
    xVals = xTmp[okayIdx]
    yVals = yTmp[okayIdx]
    colors = cmap(np.linspace(0,1,nColors))
    indices=  np.array([int(min(r*nColors,nColors-1)) for r in ratio] )
    objColor = colors[indices]
    rawImage = np.zeros((maxX,maxY))
    ax = plt.scatter(xVals,yVals,c=objColor,linewidth=0.0,alpha=0.7)
    plt.xlim([0,maxX])
    plt.ylim([0,maxY])
    return ax

def saveAsSubplot(XByStages,YByStages,c1ByStages,c2ByStages,outputDir,
                  fileFormat):
    maxX = 0
    maxY = 0
    # number of times given by columns
    nTimes = [ X.shape[1]-1 for X in XByStages]
    maxTimes = np.max(nTimes)
    numStages = len(XByStages)
    subPlots = numStages+1
    for i in range(numStages):
        maxX = max(XByStages[i].max(),maxX)
        maxY = max(YByStages[i].max(),maxY)
    rawImage = np.zeros((maxY,maxX))
    normV = matplotlib.colors.Normalize()
    cmap = plt.cm.brg
    nColors = 64
    for t in range(maxTimes):
        fig = plt.figure(dpi=200)
        rawImage[:,:] = 0
        xVals = XByStages[0][:,t]
        goodIdx = xVals.nonzero()
        rawX = np.floor(xVals[goodIdx]).astype(np.uint64)-1
        rawY = np.floor(YByStages[0][:,t][goodIdx]).astype(np.uint64)-1
        rawIntensity = (c1ByStages[0][:,t][goodIdx]+
                        c2ByStages[0][:,t][goodIdx])
        rawImage[rawY,rawX] = 1.
        counter =1
        ax = plt.subplot(1,subPlots,counter)
        plt.title('Raw')
        plt.imshow(rawImage,interpolation='sinc',cmap=cmaps.gray,
                   aspect='auto',filterrad=4,origin='lower')
        hideAxis()
        counter += 1
        for i in range(numStages):
            plt.subplot(1,subPlots,counter)
            timeIdx = min(t,nTimes[i])
            plt.title('Filter #{:d}'.format(i))
            subAx = plotSingle(XByStages[i],YByStages[i],c1ByStages[i],
                               c2ByStages[i],maxX,maxY,timeIdx,cmap,nColors)
            counter += 1
        ### XXX TODO: add in colorbar based on fluorescence
        axis = fig.add_axes([0.03, 0.03,0.9, 0.025])
        m = plt.cm.ScalarMappable(cmap=cmap)
        offsetTick = 4
        m.set_array([0,nColors-1])
        cbar = fig.colorbar(cax=axis,mappable=m,orientation='horizontal',
                            ticks=[offsetTick,nColors-1-offsetTick],
                            ticklocation='top')
        cbar.ax.set_xticklabels(['Unfolded', 'Folded'],fontsize=12)
        # save it
        plt.savefig(outputDir + fileFormat.format(t),dpi=300)


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
                    cPoint.getCheckpoint(mFile,vizIOSparse,False,stageFile)
        X.append(xSparse)
        Y.append(ySparse)
        fretC1.append(c1)
        fretC2.append(c2)
    return X,Y,fretC1,fretC2
