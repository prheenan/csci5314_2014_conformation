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

def GetPhysicsMain(goodTimes,goodFRET,goodDiff):
    colors = ['r','g','b','k','m','c']
    colorCycle = cycle(lines)
    for times,FRET,diff in zip(goodTimes,goodFRET,goodDiff):
        # this is the trace, FRET, and diff for each protein
        fig = plotUtil.pFigure()
        distanceValues = FRET**(1/6)
        normTime = util.normalize(times)
        normDist = util.normalize(distValues)
        plt.plot(normTime,distanceValues,'-',color=next(colorCycle))
    
