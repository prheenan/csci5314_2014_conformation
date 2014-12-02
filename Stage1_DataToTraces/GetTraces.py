#for matrices, use numpy
import numpy as np
# for file information, use of 
import os
# for adding path / utility functions, use sys
import sys
sys.path.append('../Utilities')
# import the util file.
import Utilities as util
# import the plotting utility file
import PlotUtilities as plotUtil
# use the regular expression generator for parsing the files
import re
# use matplotlib for plotting
import matplotlib.pyplot  as plt
# use scipy for some stats
from scipy import stats
# import the colormap
import matplotlib.cm as cm

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
    fileStr = "file [" + f + "]"
    util.ReportMessage("Reading: " + fileStr +  " in GetListOfObjectData")
    fileHandle = open(f)
    # read in the file and remove all whitespace by splitting and joining
    # by spaces
    fileBuffer = ''.join(fileHandle.read().split())
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
                       + "] points for [" + str(numMatches) + "]")
    for objNum in range(numObjects):
        thisObj = []
        for dataType in range(numDataPoints):
            toAppend = dataByGroups[dataType][objNum]
            thisObj.append(toAppend)
        dataByObjects.append(thisObj)
    util.ReportMessage("Processed: [" + str(numObjects) 
                       + "] points for [" + str(numMatches) + "]")
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
    # frames appearing is the first data point
    util.ReportMessage("ProcessData")
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
    #FRET Ratio is donor/acceptor (1/2)
    # XXX should probably check for a div/0 here...
    channel1Flatten=np.concatenate(channel1)
    channel2Flatten=np.concatenate(channel2)
    maxXAxis = max(len(channel1Flatten),len(channel2Flatten))
    fretRatio = [ np.divide(c1,c2) for c1,c2 in  zip(channel1,channel2)]
    toCompare = [channel1Flatten,channel2Flatten]
    plotUtil.compareHist('Intensity (au)','Protein Count','FRET Intensities'
                         ,maxXAxis,toCompare,
                         ['Donor Channel','Acceptor Channel'])
    return velX,velY,times,numTimes,fretRatio,MSD

def AnalyzeTraces(velX,velY,times,numTimes,fretRatio,MSD,frameRate):
    util.ReportMessage("AnalyzeTraces")
    proteinYStr = '# Proteins'
    numProteins = len(numTimes)
    numBins = max(numTimes)
    fig = plotUtil.pFigure()
    titleStr = "Raw distribution of protein appearances"
    ax = fig.add_subplot(1,1,1)
    plotUtil.histogramPlot(ax,'Frame duration of protein',proteinYStr,
                           titleStr,numTimes,numBins)
    plotUtil.saveFigure(fig,"Protein_Distribution")

    fig = plotUtil.pFigure()
    # save the MSD and R^2 of the MSD
    numStats = 2
    msdMatrix = np.zeros((numProteins,numStats))
    plotCount = 1
    numPlots = 2
    msdAx = plt.subplot(numPlots,1,plotCount)
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
    rawRsqVals = msdMatrix[:,1]
    goodIndices = np.where(rawRsqVals > 0)[0]
    goodMsds = np.take(msdMatrix,goodIndices,axis=0)
    msdVals =  goodMsds[:,0]
    # plot the histogram of MSDs
    ax = fig.add_subplot(numPlots,1,plotCount)
    numBins = max(msdVals)
    numProteins = len(msdVals)
    axTmp = plotUtil.histogramPlot(ax,'Diffusion Coeff (pixels^2/second)',
                                   proteinYStr,
                                   'Histogram of Protein Hiffusion Coeffs',
                                   msdVals,numBins)
    plotCount+=1
    # plot the rSquared values
    RSqVals = goodMsds[:,1]*100
    # RSQ is between 0 and 1, multiply by 100 and fit with 100 bins for 'simple' normalization
    numBinsRsq= 100 
    ax = fig.add_subplot(numPlots,1,plotCount)
    # use 100 bins for each of the RSq values
    plotUtil.histogramPlot(ax,'R Squared Coeff',proteinYStr,
                            'Histogram of Protein RSq',RSqVals,numBinsRsq)
    plotCount += 1
    # next, we plot a comparison of the 'raw', 'valid' (D with uncertainty),
    # and 'processed'
    # (D fit with RSQ > cutoffRsq)
    cutoffRsq = 0.8
    plt.tight_layout()
    plotUtil.saveFigure(fig,"MSD")
    
    titleStr = proteinYStr
    xStr = 'Frames Appearing'
    xLimit = [1,max(numTimes)]
    bestIndices = np.where(rawRsqVals > cutoffRsq)[0]
    # make labels for each of the indices ('Raw' is assumed...)
    compLabel = ["All Proteins","Valid Diffusion Coeffs",
                        ("RSQ > {:.3f}".format(cutoffRsq))]
    comparisonIndices = [goodIndices,bestIndices]
    plotUtil.comparisonPlot(xStr,proteinYStr,titleStr,xLimit,numTimes,
                            comparisonIndices,compLabel)
    # XXX TODO: compare other options (e.g. x velocity, y velcocity, etc)
    # XXX TODO: try for all files, need a better way to store
    # return all the diffusion coeffs, as well as the 'best' indices we found
    return msdMatrix[:,0],bestIndices
    
def GetTracesMain(fileNameList,frameRate=0.1):
    for f in fileNameList:
        if (not os.path.isfile(f)):
            util.ReportError(True,("Could not find file [" +f +"]"))
        else:
        # POST: tmpFile is a valid filename
            dataByObjects, dataByGroups = GetListOfObjectData(f)
            velX,velY,times,numTimes,fretRatio,MSD = ProcessData(dataByObjects,frameRate)
            diffCoeffs,bestIndices = \
                                     AnalyzeTraces(velX,velY,times,numTimes,
                                                   fretRatio,MSD,frameRate)
            
            # XXX only works one at a time for now
            atGoodIndices =\
            util.saveAllAtIndices([times,numTimes,fretRatio,MSD,diffCoeffs],
                                  ['Per_Protein_Frames_Appearing',
                                   'Per_Protein_Num_Frames',
                                   'Per_Protein_Fret_Ratio','MSD_Per_Protein',
                                   'Per_Protein_Diff_Coeff'],
                                  ['Post_Stage1_Analysis'],bestIndices)
            goodTimes , ignore, goodFRET , ignore,goodDiff = atGoodIndices
            return goodTimes,goodFRET,goodDiff
