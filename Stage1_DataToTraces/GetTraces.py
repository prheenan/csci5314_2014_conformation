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

def ProcessData(dataByObject,frameRate=0.1):
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
    velX = getVelocity(xVals,deltaTimes)
    velY = getVelocity(yVals,deltaTimes)
    # channel 1 (FRET donor) is the third
    donorIndex = 2
    channel1 = [ obj[donorIndex][meanIndex] for obj in dataByObject]
    # channel 2 (FRET acceptor) is the fourth
    acceptorIndex = 3
    channel2 = [ obj[acceptorIndex][meanIndex] for obj in dataByObject]
    fretRatio = [ np.divide(c2,c1) for c1,c2 in  zip(channel1,channel2)]
    return velX,velY,times,numTimes,fretRatio

def AnalyzeTraces(velX,velY,times,numTimes,fretRatio):
    mSource = "Step1::AnalyzeTraces"
    numProteins = len(numTimes)
    numBins = numProteins/100
    fig = plt.figure()
    numPlots = 1
    plotCount = 1
    ax = fig.add_subplot(numPlots,1,1)
    vals, xBins, p= ax.hist( numTimes, numBins, normed=1,
                             facecolor='green', alpha=0.75)
    plt.title("Raw distribution of protein appearances, N={:d}".format(
        numProteins))
    plt.ylabel('Normalized protein histogram (Fraction of total N) ')
    plt.xlabel('Frames in which protein appeared')
    plt.yscale('log', nonposy='clip')
    yLimNorm = [1/numProteins,1.]
    yLimits= plt.ylim(yLimNorm)
    plt.grid(True,which='major',color='b')
    plt.grid(True,which='minor',color='r')
    # create a twin axis for the absolute count 
    yLimAbs = np.multiply(yLimNorm,numProteins)
    util.secondAxis(ax,"Absolute Protein Histogram",yLimAbs)
    plotCount += 1
    util.saveFigure(mSource,"Protein_Distribution",fig)

def GetTracesMain(fileNameList):
    for f in fileNameList:
        if (not os.path.isfile(f)):
            util.ReportError(True,("Could not find file [" +f +"]"),mSource)
        else:
        # POST: tmpFile is a valid filename
            dataByObjects, dataByGroups = GetListOfObjectData(f)
            velX,velY,times,numTimes,fretRatio = ProcessData(dataByObjects)
            AnalyzeTraces(velX,velY,times,numTimes,fretRatio)
