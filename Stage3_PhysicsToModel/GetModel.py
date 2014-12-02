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

def GetModelMain(residenceTimes,diffCoeffs,frameRate = 0.1):
    sortedTimes = np.sort(residenceTimes)
    numProteins = len(residenceTimes)
    diffTimes = np.gradient(np.unique(sortedTimes))
    # get the minimum step we care about, limited by our frame rate
    minTimeStep = max(min(diffTimes),frameRate)
    timeGrid = np.arange(0,max(sortedTimes),minTimeStep)
    # count the number unfolded
    countUnfoldedCummu = [sum(sortedTimes > t) for t in timeGrid]
    hist,bins = np.histogram(sortedTimes, bins=timeGrid)
    goodIndices = np.where(hist > 0)[0]
    hist = hist[goodIndices]
    bins = bins[goodIndices]
    # get the proportion unfolded
    proportionUnfoldedCummu = np.divide(countUnfoldedCummu,numProteins)
    # all that remains is plotting.
    fig = plotUtil.pFigure()
    numPlots = 2
    pltCounter = 1
    # plot the residence time distribution
    fig.add_subplot(numPlots,1,pltCounter)
    plt.semilogy(bins,hist/numProteins,'ro-')
    plt.ylabel('Probability to fold during time window t')
    plt.title('Residence Time Distribution: P(t) vs t')
    pltCounter += 1
    #plot the cummulative res time dist (CRTD)
    fig.add_subplot(numPlots,1,pltCounter)
    plt.semilogy(timeGrid,proportionUnfoldedCummu,'ro-')
    plt.xlabel('Time (s)')
    plt.ylabel('Probability to fold after time t')
    plt.title('Cummulative Residence Time Distribution: P(tau>t) vs t')
    pltCounter += 1

    # save the figure
    plotUtil.saveFigure(fig,'RTD')

    
    
