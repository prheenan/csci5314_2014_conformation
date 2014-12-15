# use numpy for array operations
import numpy as np
# import sys for file and path stuff
import sys
# add the utilities and model helper so we can fit
sys.path.append('../Utilities')
sys.path.append('../Stage4_FindModelVariables')

# import the util file.
import Utilities as util
# import the model utility file
import ModUtilities as modUtil
# import the plotting utility file
import PlotUtilities as plotUtil
import getModelVariables as modHelp

# import cycling..
from itertools import cycle
import matplotlib.pyplot as plt

def getIndicesWithChanges(counts):
    cumDiff = np.abs(np.diff(counts))
    idx = (np.where(cumDiff != 0)[0])
    return idx

def mSecondAxis(ax,numProteins):
    # get the proper lower bound from the normalized y axis
    minY = numProteins * np.min(ax.get_ylim())
    plotUtil.secondAxis(ax,'Absolute Number of Proteins',[minY,numProteins])


def GetModelMain(residenceTimes,diffCoeffs,modelNums,frameRate = 0.1,plotRTD=False):
    modParams = []
    modStdev = []
    modRSQ = []
    # get the 'overall' parameters for the time grid
    ravelTimes = np.concatenate(residenceTimes).flatten()
    numProteins = len(ravelTimes)
    sortedRavelTimes = np.sort(ravelTimes)
    diffTimes = np.gradient(np.unique(sortedRavelTimes))
    minTimeStep = max(np.min(diffTimes),frameRate)
    numTrials = len(residenceTimes)
    # get the time grid based on the parameters.
    timeGrid = np.arange(0,np.max(ravelTimes),minTimeStep)
    numTimes = len(timeGrid)
    # for each trial, we will
    countGrid = np.zeros((numTimes,numTrials))
    # for the cummulative and count, might as well wait until everything
    # is averaged to get the totals...
    cummulativeGrid = np.zeros((numTimes,1))
    countStd = np.zeros((numTimes,1))
    cummuStd = np.zeros((numTimes,1))
    # next, count the numbers in each trial...
    for trialIdx,timeTrial in enumerate(residenceTimes):
        sortedTimes = np.sort(timeTrial)
        hist,bins = np.histogram(sortedTimes, bins=timeGrid)
        countGrid[:-1,trialIdx] = hist
    # get the average counts within the windows and such...
    countAvg = np.sum(countGrid,axis=1)
    nonZeroIdx = np.where(countAvg > 0)[0]
    countAvg = countAvg[nonZeroIdx]
    timeRtd = timeGrid[nonZeroIdx]
    countErr = np.std(countGrid[nonZeroIdx,:],axis=1)
    # get the proportion unfolded, cummulatively (cummu residence time dist)
    # XXX could label the dots with only a single protein...
    RTD = countAvg/numProteins
    CRTD = 1.-np.cumsum(RTD,dtype=np.float64)
    rtdErr = countErr/numProteins
    crtdErr = np.cumsum(rtdErr,dtype=np.float64)


        # POST: we have all the final data we need. Go ahead and save.
    util.saveAll([timeGrid,timeRtd,CRTD,RTD],
                 [util.IO_time_CRTD,util.IO_time_RTD,util.IO_CRTD,util.IO_RTD],
                 [util.IO_Stage3Folder])

    for modelNumber in modelNums:
        modelLabels = []
        modelGuess = []
        # create lists to make all the tau and fs.
        freqLabels = ['f0']
        tauLabels = ['tau0']
        # start off with an initial guess
        tauGuess = [1]
        freqGuess = [1]
        #modelNumber-1 gives number of taus
        #modelNumber-1 gives number of taud (but only fit modelNumber-1)
        for i in range(modelNumber-1):
            freqLabels.append('f' + str(i))
            tauLabels.append('tau'+str(i))
            tauGuess.append(tauGuess[i]/10)
            freqGuess.append(freqGuess[i]/10)
        if (modelNumber > 1):
            modelGuess.extend(freqGuess[1:])
        modelGuess.extend(tauGuess)
        modelLabels.extend(freqLabels)
        modelLabels.extend(tauLabels)
        modelUnits = []
        for m in modelLabels:
            if ['A' in m or 'f' in m]:
                modelUnits.append('au')
        else:
            modelUnits.append('s')
        modelProb,fitUnordered,stdUnordered,rsq =modHelp.\
                                    getModelVariables(timeRtd,CRTD,modelNumber,modelGuess)
        # fit parameters may be out of order (don't enforce that in the model)
        #  go ahead and sort them. This is just a labelling thing, and doesn't change the physics
        numParams = len(fitUnordered)
        middle = numParams/2
        indicesTau = np.arange(middle,numParams,dtype=np.int32)
        indicesF = np.arange(0,middle,dtype=np.int32)
        taus = fitUnordered[indicesTau]
        sortIndices = np.argsort(taus)
        taus = taus[sortIndices]
        fVals = fitUnordered[indicesF][sortIndices]
        stdevTaus = stdUnordered[indicesTau][sortIndices]
        stdevFvals = stdUnordered[indicesF][sortIndices]
        fitParams = np.array( ((fVals,taus)) ).flatten()
        # second is stdev
        stdev = np.array( ((stdevFvals,stdevTaus)) ).flatten()
        modStr = modUtil.modelString("W1 Fit",fitParams,stdev,modelLabels,
                                     modelUnits,rsq)

        util.saveAll([fitParams,stdev,[rsq]],
                     [util.getModelFile(util.IO_model,modelNumber),
                      util.getModelFile(util.IO_model_stdevs,modelNumber),
                      util.getModelFile(util.IO_model_rsq,modelNumber)],
                     [util.IO_Stage3Folder])
        modParams.append(fitParams)
        modStdev.append(stdev)
        modRSQ.append([rsq])

        # all that remains is plotting.
        fig = plotUtil.pFigure()
        if (plotRTD):
            numPlots = 2
        else:
            numPlots = 1
        pltCounter = 1
        # plot the residence time distribution (RTD)
        #plot the cummulative res time dist (CRTD)
        axCRTD = fig.add_subplot(numPlots,1,pltCounter)
        axCRTD.set_yscale('log')
        plt.title('Cummulative Residence Time Distribution [CRTD, P(t >t_u)] vs t_u is exponential')
        plt.errorbar(timeRtd,CRTD,rtdErr,fmt='ro',label='Experimental CRTD')
        plt.plot(timeRtd,modelProb,'b-',label=modStr)
        # re-center the plot around the real data, so that the messed up model is OK.
        plt.ylim( 0.95*np.min(CRTD),np.max(CRTD)*1.05)
        plt.xlabel('Time (s)')
        plt.ylabel('P(t)')
        plt.legend()
        mSecondAxis(axCRTD,numProteins)
        nStr = plotUtil.getNStr(numProteins)
        pltCounter += 1
        if (plotRTD):
            axRTD = fig.add_subplot(numPlots,1,pltCounter)
            axRTD.set_yscale('log')
            plt.errorbar(x=timeRtd,y=RTD,yerr=rtdErr,fmt='ro-')
            mSecondAxis(axRTD,numProteins)
            plt.ylabel('Probability to fold during time window t')
            plt.title('Residence Time Distribution: P(t) vs t for' + nStr + ' is noisier than CRTD')
            pltCounter += 1
        plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0)
        # save the figure
        plotUtil.saveFigure(fig,'RTD_model' +str(modelNumber))
    return modParams,modStdev,modRSQ

    
    
