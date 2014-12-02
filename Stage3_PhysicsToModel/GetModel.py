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
    print(residenceTimes)
    sortedTimes = np.sort(residenceTimes)
    print(sortedTimes)
    numProteins = len(residenceTimes)
    diffTimes = np.gradient(np.unique(sortedTimes))
    # get the minimum step we care about, limited by our frame rate
    minTimeStep = max(min(diffTimes),frameRate)
    timeGrid = np.arange(0,np.max(sortedTimes),minTimeStep)
    # count the number unfolded
    countUnfoldedCummu = [sum(sortedTimes > t) for t in timeGrid]
    hist,bins = np.histogram(sortedTimes, bins=timeGrid)
    goodIndices = np.where(hist > 0)[0]
    hist = hist[goodIndices]
    bins = bins[goodIndices]
    # get the proportion unfolded, cummulatively (cummu residence time dist)
    CRTD = np.divide(countUnfoldedCummu,numProteins)
    RTD = hist/numProteins

    # POST: we have all the final data we need. Go ahead and save.
    util.saveAll([timeGrid,bins,CRTD,RTD],
                 [util.IO_time_CRTD,util.IO_time_RTD,util.IO_CRTD,util.IO_RTD],
                 [util.IO_Stage3Folder])

    # all that remains is plotting.
    fig = plotUtil.pFigure()
    numPlots = 2
    pltCounter = 1
    # plot the residence time distribution (RTD)
    fig.add_subplot(numPlots,1,pltCounter)
    plt.semilogy(bins,RTD,'ro-')
    plt.ylabel('Probability to fold during time window t')
    plt.title('Residence Time Distribution: P(t) vs t')
    pltCounter += 1
    #plot the cummulative res time dist (CRTD)
    fig.add_subplot(numPlots,1,pltCounter)
    plt.semilogy(timeGrid,CRTD,'ro-')
    plt.xlabel('Time (s)')
    plt.ylabel('Probability to fold after time t')
    plt.title('Cummulative Residence Time Distribution: P(tau>t) vs t')
    pltCounter += 1

    # save the figure
    plotUtil.saveFigure(fig,'RTD')

    
    
