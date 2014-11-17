#for matrices, use numpy
import numpy as np
# for file information, use of 
import os
# for adding path / utility functions, use sys
import sys
sys.path.append('../Utilities')
# import the util file.
import Utilities as util
# use the regular expression generator for parsing the files
import re
# use matplotlib for plotting
import matplotlib.pyplot  as plt
# use scipy for some stats
from scipy import stats

# create our pattern using the verbose flag for commenting.
    # http://stackoverflow.com/questions/9507173/python-regex-re-compile-match:
    # r flag: tells python to interpret escape chaarcters literally. 
    # triple quotes: prevents premature ending  by single or double quotes

singleValuePattern = re.compile(
     r"""
      \{          # open curly
      ([^{,]+) # data1
      ,        # a comma
      ([^},]+) # data2: anything not a curly or comma
      \}
      """,re.VERBOSE)
      
singleObjectStr ="""
    \{                          # single object: match multiple values
    (?:                         # non capturing group...
    \{                          # lowest level: match single value for choice: 
    (?:                         # non capturing group (just join together our choices)
    (?:["a-zA-Z]+,["a-zA-Z]+) # 1, match word,word
    |(?:[+-]?\d+\.\d*,        # 2a, match float,float, allowing for +/- and "0."
     [+-]?\d+\.\d*)           # 2b, match float,float, allowing for +/- and "0."
    |(?:\d+,\d+)              # 3, match int, int
    )
    \},?)+                      # end of single value
    \}                      # end of single object
"""

singleObjectPattern =re.compile(singleObjectStr,re.VERBOSE)

searchPattern = re.compile(
    r"""
    (
    \{                         # single dataset
    (?:                          # non capturing group...
    """
    + 
    singleObjectStr 
    + """
    ,?)+
    \})                          # end of single dataset
    """,re.VERBOSE)

def GetListOfObjectData(f):
    mSource = "Step1::GetListOfObjectData"
    fileStr = "file [" + f + "]"
    util.ReportMessage("Reading: " + fileStr,mSource)
    fileHandle = open(f)
    # read in the file and remove all whitespace by splitting and joining
    # by spaces
    fileBuffer = ''.join(fileHandle.read().split())
    util.ReportMessage("Matching: " + fileStr,mSource)
    groups=searchPattern.findall(fileBuffer)
    numMatches = len(groups)    
    if (groups):
        dataByGroups = []
        for i in range(numMatches):
            match = groups[i]
            objects = singleObjectPattern.findall(match)
            allObjects = []
            for o in objects:
                pattern = singleValuePattern.split(o)
                # split returns the end of the string
                # (which we don't carew about, as long as it is empty)
                # as the last argument and the start as the first
                numPatterns = len(pattern)
                step = 3
                objFirstVals = [float(pattern[i])
                                for i in range(1,numPatterns,step)]
                objSecondVals = [float(pattern[i])
                                 for i in range(2,numPatterns,step)]
                # POST: have all the first and second values for this object
                toAdd =  [objFirstVals, objSecondVals]
                allObjects.append(toAdd)
                dataByGroups.append(allObjects)
    else:
        util.ReportError(True,
                         ("Could not match file [" + f +"] to regex"),
                         mSource)

    # POST: at least one macth
    # XXX check that lengths are all the same, proper numberf of groups
    numObjects = len(dataByGroups[1])
    numDataPoints = numMatches
    dataByObjects = []
    util.ReportMessage("Processing: [" + str(numObjects) 
                       + "] points for [" + str(numMatches) + "]",mSource)
    for objNum in range(numObjects):
        thisObj = []
        for dataType in range(numDataPoints):
            toAppend = dataByGroups[dataType][objNum]
            thisObj.append(toAppend)
        dataByObjects.append(thisObj)
    return dataByObjects, dataByGroups

def getVelocity(xArr,deltaTimeArr):    
    return [ np.gradient(x) / deltaTimeArr[i] 
             if np.all(deltaTimeArr[i]) >0 else 0 for i,x in enumerate(xArr)]

def sqDeltaSum(vals):
    #gives the cummulative sum of the sqaured values of val 
    # good for calculating the MSD. returns a vector of size N-1
    delta = np.subtract(vals,vals[0])
    return np.cumsum(delta*delta)

def ProcessData(dataByObject,frameRate):
    # frameRate is in second
    mSource = "Step1::ProcessData"
    # frames appearing is the first data point
    timeIndex = 0
    meanIndex = 0
    times = [ np.multiply(obj[timeIndex][meanIndex],frameRate)
              for obj in dataByObject]
    deltaTimes = [ np.gradient(t) if len(t) > 1 else 0 for t in times ]
    numTimes = [len(t) for t in times]
    # xy is the second data point
    xyIndex = 1
    xIndex = 0
    yIndex = 1
    xVals = [ obj[xyIndex][xIndex] for obj in dataByObject]
    yVals = [ obj[xyIndex][yIndex] for obj in dataByObject]
    MSD = [ 1/(np.arange(1,1+len(x))) * (sqDeltaSum(x) + sqDeltaSum(y))
            if len(x) > 1 else np.array([0]) for x,y in zip(xVals,yVals) ]
    velX = getVelocity(xVals,deltaTimes)
    velY = getVelocity(yVals,deltaTimes)
    # channel 1 (FRET donor) is the third
    donorIndex = 2
    channel1 = [ obj[donorIndex][meanIndex] for obj in dataByObject]
    # channel 2 (FRET acceptor) is the fourth
    acceptorIndex = 3
    channel2 = [ obj[acceptorIndex][meanIndex] for obj in dataByObject]
    fretRatio = [ np.divide(c2,c1) for c1,c2 in  zip(channel1,channel2)]
    return velX,velY,times,numTimes,fretRatio,MSD

def AnalyzeTraces(velX,velY,times,numTimes,fretRatio,MSD,frameRate):
    mSource = "Step1::AnalyzeTraces"
    proteinYStr = '# Proteins'
    numProteins = len(numTimes)
    numBins = max(numTimes)-1
    fig = plt.figure()
    titleStr = "Raw distribution of protein appearances, N={:d}".format(
        numProteins)
    ax = fig.add_subplot(1,1,1)
    util.histogramPlot(ax,'Frame duration of protein',proteinYStr,
                       titleStr,numTimes,max(numTimes))
    util.saveFigure(mSource,"Protein_Distribution",fig)

    fig = plt.figure(figsize=(16, 12),dpi=500)
    # save the MSD and R^2 of the MSD
    numStats = 2
    msdMatrix = np.zeros((numProteins,numStats))
    plotCount = 1
    numPlots = 3
    plt.subplot(numPlots,1,plotCount)
    allMSDs = np.concatenate(MSD)
    allTimes= np.concatenate(times)
    for i in range(numProteins):
        # get the X and Y values to fit...
        tmpMSD = MSD[i]
        # multiple the times by four, per the diffusion formulae
        tmpTimes = times[i]*2
        tmpTimes -= tmpTimes[0]
        if (len(tmpTimes) < 3):
            # nothing valid here. set everything to 0 and flag 0 RSQ later
            # we need at least three values to be able to get a Diffusion Coefficient
            # and an uncertainty.
            msdMatrix[i,:] = 0
            continue
        # linear fit
        deg = 1
        polyVals = np.polyfit(tmpTimes,tmpMSD,deg)
        polyFunc = np.poly1d(polyVals)
        fitVals = polyFunc(tmpTimes)
        slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(fitVals,tmpMSD)
        rSquared = r_value**2
        # MSD (slope) is given by the slope, first coeff returned
        diffCoeff = polyVals[0]
        if (diffCoeff < 0):
            # ignore diffusion coefficients less than 0
            continue
        msdMatrix[i,0] = diffCoeff
        msdMatrix[i,1] = rSquared
       # s: size in points^2. It is a scalar or an array of the same
       #length as x and y.
        plt.plot(tmpTimes,tmpMSD,'r.',markersize=0.5)
        plt.plot(tmpTimes,fitVals,'b-')
    plt.ylim(1,max(allMSDs))
    plt.xlim(frameRate,max(allTimes))
    plt.yscale('log', nonposy='clip')
    plt.xscale('log', nonposy='clip')
    plt.ylabel('Mean Sq Distance, (pixels/second)')
    plt.xlabel('Time (Seconds)')
    plotCount += 1
    goodIndices = np.where(msdMatrix[:,1] > 0)
    goodMsds = np.take(msdMatrix,goodIndices,axis=0)[0]
    msdVals =  goodMsds[:,0]
    # plot the histogram of MSDs
    ax = fig.add_subplot(numPlots,1,plotCount)
    numBins = max(msdVals)
    ax = util.histogramPlot(ax,'MSD (pixels/s)',proteinYStr,
                            'Histogram of Protein MSD',msdVals,numBins)
    ax.set_xlim(1,numBins)
    ax.set_xscale('log')
    plotCount+=1
    # plot the rSquared values
    RSqVals = goodMsds[:,1]*100
    ax = fig.add_subplot(numPlots,1,plotCount)
    # use 100 bins for each of the RSq values
    ax = util.histogramPlot(ax,'R Squared Coeff',proteinYStr,
                            'Histogram of Protein RSq',RSqVals,100)
    plt.tight_layout()
    util.saveFigure(mSource,"MSD",fig)
        

def GetTracesMain(fileNameList,frameRate=0.1):
    for f in fileNameList:
        if (not os.path.isfile(f)):
            util.ReportError(True,("Could not find file [" +f +"]"),mSource)
        else:
        # POST: tmpFile is a valid filename
            dataByObjects, dataByGroups = GetListOfObjectData(f)
            velX,velY,times,numTimes,fretRatio,MSD = ProcessData(dataByObjects,frameRate)
            AnalyzeTraces(velX,velY,times,numTimes,fretRatio,MSD,frameRate)
            # XXX only works one at a time for now
            return velX,velY,times,fretRatio
