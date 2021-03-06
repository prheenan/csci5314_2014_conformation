# force floating point division. Can still use integer with //
from __future__ import division
# use numpy for array operations
import numpy as np
# import sys for file and path stuff
import sys
sys.path.append('../Utilities')
# import the util file.
import Utilities as util
# import the plotting utility file
import PlotUtilities as plotUtil
# import cycling..
from itertools import cycle
import matplotlib.pyplot as plt
import scipy.cluster.vq as cluster 
from scipy import interpolate
# use math for isnan and such
import math
sys.path.append('../Base')
import Filter

def convertTimeToIndex(times,frameRate=0.1):
    return np.round(times / frameRate)

def getMinimumTime(values,times,value,lessThanBool=False):
    minimumTimes = -1 * np.ones((len(values),1))
    # the times to interpolate the times. XXX probably 
    # want to make this a constant of some kind.
    numPoints = 1000
    # get the bounds for interpolation
    # get the interpolation grid
    for i,listV in enumerate(values):
        timeRaw = times[i]
        timeGrid = np.linspace(0,np.max(timeRaw),numPoints)
        listTmp = np.interp(timeGrid,timeRaw,listV)
        # need to make sure we start below the lower bound
        # or send above the upper bound
        lessThan =(listTmp < value)
        greaterThan = (listTmp > value)        
        if (lessThanBool):
            toTest = lessThan
        else:
            toTest= greaterThan
        # if our condition isn't true for any time,
        # ignore the protein (ie: it wasn't ever less than
        # the lower bound or greater than the upper bound)
        if not np.any(toTest):
            continue
        # look for the lowest squared difference from the desired value.
        diffArr = (listTmp-value)**2
        firstIndex = np.argmin(diffArr)
        # get the best interpolated time.
        bestTime = timeGrid[firstIndex]
        minimumTimes[i] = bestTime
    return minimumTimes

def getDifferentialTime(valuesBefore,valuesAfter):
    # returns when we have a well-defined unfolding...
    # also returns the indices of the proteins which were used.
    numVals = len(valuesBefore)
    indices=  []
    diffTime = -1. * np.ones(numVals)
    for i in range(numVals):
        first = valuesBefore[i]
        second = valuesAfter[i]
        if (first < 0 or second < 0):
            continue
        if (first >= second):
            continue
        # POST: non negative, values 'match'
        indices.append(i)
        # append the proper residence time.
        # all others will be -1...
        diffTime[i] = second-first
    return diffTime,indices

def getDistances(goodFRET,goodTimes):
    # goodFRET is defined by the API as the FRET ratio from GetTraces.
    # This essentially means it has been corrected for objects which have
    # ill-defined trajectories or difusion coeffs
    # RETURNS: the 'good' distances and time found, as well as the indexing array
    distances = []
    times = []
    indices = []
    nodeIndices = []
    count = 0
    for t,FRET in zip(goodTimes,goodFRET):
        # this is the trace, FRET, and diff for each protein
        # force some reason we have negative data? (Figure this out)
        # only look at where fret is finite and ge 0.
        goodIndices = np.where( (FRET >= 0) 
                                &
                                np.isfinite(FRET))[0]
        if ((len(goodIndices) == 0) or 
            (min(goodIndices) == max(goodIndices))):
            indices.append([])
            continue
        betterFRET = FRET[goodIndices]
        betterTime = t[goodIndices]
        betterTime -= np.min(betterTime)
        distanceValues = betterFRET**(1/6)
        distances.append(np.array(distanceValues))
        indices.append(np.array(goodIndices))
        times.append(np.array(betterTime))
        nodeIndices.append(count)
        count += 1
    return np.array(distances),np.array(times),\
        np.array(indices),np.array(nodeIndices)



def GetPhysicsMain(goodTimes,goodFRET,goodDiff):
    source = 'Step2::GetPhysics'
    util.ReportMessage("Starting",source)
    colors = ['r','g','b','k','m','c','y','0.33','0.66']
    colorCycle = cycle(colors)
    count = 0
    fig = plotUtil.pFigure()
    # XXXfill all these in! based on video size
    frameRate = 0.1
    maxNumTimes = 30*10
    distances,times,definedDistancesIdx,nodeIdx =\
                                getDistances(goodFRET,goodTimes)
    # get just the 'nodes' with valid valued.
    # flatten the distances to get a notion of the 'all time'
    # distance information. We can use this and kmeans to find a 'folded
    # and unfolded state
    
    distMaxDiff = [np.max(d) - np.min(d) for d in distances]
    flattenDistances = np.concatenate(distances)
    # use two clusters; folded and unfolded
    numClusters = 2
    # lots of iterations (this number seems to work well; the 'smooth'
    # running time / convergence (?) of kmeans is polynomial
    # http://en.wikipedia.org/wiki/K-means_clustering
    numIters = int(1e3)
    clusters,ignore = cluster.kmeans(flattenDistances,numClusters,iter=numIters)
    # the clusters give us the 'folded' and 'unfolded' groups. 
    #between those, we have a fairly undefined state.
    folded = min(clusters)
    unfolded = max(clusters)
    # go ahead and start plotting so we can use the histograms to 
    # get the average change. this will be used insted of a cluster 
    fig = plotUtil.pFigure()
    numPlots = 3
    plotCount = 1
    fretLabel = 'FRET d Distance (arb)'
    # plot the distribution of maximum FRET distances
    ax = plt.subplot(numPlots,1,plotCount)
    xlims = [np.min(flattenDistances),np.max(flattenDistances)]

    # N here is the normalized number in each bin, so this gives
    # the normaized average (good for plotting)
    medianDist = np.median(distMaxDiff)
    meanDist = np.mean(distMaxDiff)
    bestDistChangeAbs = medianDist/2
    medianStr = "Half Median Change: {:.2f}".format(bestDistChangeAbs)
    # in order to actually cluster our data, we need a real, absolute number
    # for the average. de-normalize it.
    bestDistChangeNorm  = plotUtil.normalizeTo(distMaxDiff,bestDistChangeAbs)

    ax,N,bins = plotUtil.histogramPlotFull(ax,fretLabel,'# Proteins',
                           'Maximum FRET change, ' + medianStr +',',distMaxDiff,
                                len(flattenDistances)/100,True,True)
    plt.axvline(bestDistChangeNorm,color='r',label=medianStr)

    plt.legend()

    belowIsFolded = folded
    aboveIsUnfolded = folded + bestDistChangeAbs

    # now get the folded and unfolded dynamics, based on the cutoffs from above
    # the folded times must be less than what we saw below
    foldedArr = getMinimumTime(distances,times,
                               belowIsFolded,True)
    unfoldedArr = getMinimumTime(distances,times,
                                 aboveIsUnfolded,False)

    unfoldingTimes, definedUnfoldingIdx = \
                        getDifferentialTime(foldedArr,unfoldedArr)
    # create our new diffusion coefficients by first taking the ones with
    # well-defined FRET, then taking the ones with well-defined folding.
    goodIndices = nodeIdx[definedUnfoldingIdx]
    indexArr = [nodeIdx,definedUnfoldingIdx]

    filterUnfolding = util.takeSubset(unfoldingTimes,indexArr)
    goodTime = util.takeSubset(times,indexArr)
    plotCount +=1
    # plot the FRET distances
    ax = plt.subplot(numPlots,1,plotCount)
    plotUtil.histogramPlot(ax,fretLabel,'# Frames',
                           'FRET Distance histogram',flattenDistances,
                           len(flattenDistances)/100,True,True)
    # plot guiding lines for the two clusters we found
    normalClusters = plotUtil.normalizeTo(flattenDistances,
                                          [belowIsFolded,aboveIsUnfolded])
    foldedNorm = min(normalClusters)
    unfoldedNorm = max(normalClusters)
    lowerNorm = [0,foldedNorm, unfoldedNorm]
    upperNorm = [foldedNorm,unfoldedNorm,100]

    colors = ['r','y','b']
    labels = ['folded','intermediate','unfolded']
    alphaVal = 0.25
    for low,up,col,lab in zip(lowerNorm,upperNorm,colors,labels):
        plt.axvspan(low, up, color=col, alpha=alphaVal,label=lab)
    plt.legend()
    plotCount += 1

    # plot something like the residence time
    ax = plt.subplot(numPlots,1,plotCount)
    plotUtil.histogramPlot(ax,'Unfolding time distribution','# Proteins',
                           'Unfolding time (seconds) ',filterUnfolding,
                           len(filterUnfolding)/100,True,True)


    plotUtil.saveFigure(fig,'FRET_and_ResidenceTimeHistograms')
    # return the good unfolding times and differential coefficients
    return unfoldingTimes,goodIndices

class GetPhysics(Filter.DataFilter):
    def __init__(self,data,outFile,forced):
        super(GetPhysics,self).__init__(data,outFile,forced)
    def getFilterIdx(self):
        goodTimes,goodFRET,goodDiff = self.get([Filter.prop_time,
                                                Filter.prop_ratio,
                                                Filter.prop_diff])
        residenceTime,goodIndices= \
                                GetPhysicsMain(goodTimes,goodFRET,goodDiff)
        # add in the residence times (unfiltered, we will filter later).
        self.setElement(Filter.prop_resi,residenceTime) 
        return goodIndices
        
