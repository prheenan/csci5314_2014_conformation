# import utilities for error repoorting etc
import Utilities as util
# use matplotlib for plotting
#http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
import matplotlib
matplotlib.use('AGG') 
import matplotlib.pyplot  as plt
# import numpy for array stuff
import numpy as np

def normalizedBins(array):
    return 100

def histNormalize(histIn,minV=-1,maxV=-1):
    # histograms are only truly normal if we use bins of size 0.1
    # so, normalize the values to between 0 and 1, use size of since 0.1
    if (minV < 0):
        minV = min(histIn)
    if (maxV < 0):
        maxV = max(histIn)
    # normalize the histogram and center between 0 and 1 for normal plotting
    normalizedHist = 100 *(histIn - minV)/(maxV-minV)
    return normalizedHist,minV,maxV


def histogramPlot(ax,xlabelStr,ylabelStr,titleStr,data,numBins,
                  secondX = False, absoluteNum=-1):
    n, bins, patches =ax.hist( data, numBins, normed=True,stacked=True,
                                facecolor='green', alpha=0.75)
    if (absoluteNum < 0):
        absoluteNum = len(data)
    t = plt.title(titleStr +  getNStr(absoluteNum))
    if (secondX):
        t.set_y(1.4)
    ax.set_ylabel('Norm ' + ylabelStr)
    ax.set_xlabel(xlabelStr)
    ax.set_yscale('log', nonposy='clip')
    # limits are half of an individual to slightly more than the entire population
    yLimNorm = [0.5/absoluteNum,2]
    yLimits= plt.ylim(yLimNorm)
    ax.grid(True,which='major',color='b')
    ax.grid(True,which='minor',color='r')
    # create a twin axis for the absolute count 
    yLimAbs = np.multiply(yLimNorm,absoluteNum)
    secondAxis(ax,"Abs " + ylabelStr ,yLimAbs)
    return ax


def secondAxis(ax,label,limits,secondY =True):
    #copies the first axis (ax) and uses it to overlay an axis of limits
    if(secondY):
        ax2 = ax.twinx()
        ax2.set_yscale(ax.get_yscale(), nonposy='clip')
        ax2.set_ylim(limits)
        ax2.set_ylabel(label)
    else:
        ax2 = ax.twiny()
        ax2.set_xscale(ax.get_xscale(), nonposy='clip')
        ax2.set_xlim(limits)
        ax2.set_xlabel(label)

def saveFigure(source,fileName,figure):
    # source : where to save the output iunder the output folder
    # filename: what to save the file as. automagically saved as high res pdf
    defaultOutputDir = util.defaultOutputDir()
    util.ensureDirExists(defaultOutputDir)
    sourceOutputDir = defaultOutputDir + source +'/'
    util.ensureDirExists(sourceOutputDir)
    fullPath = sourceOutputDir + fileName
    plt.tight_layout(False)
    formatStr = "png"
    figure.savefig(fullPath + '.' + formatStr,format=formatStr, 
                   dpi=figure.get_dpi())

def pFigure(ySize=8,xSize=6,dpiArg=100):
    return  plt.figure(figsize=(ySize,xSize),dpi=dpiArg)

def getNStr(n,space = " "):
    return space + "n={:d}".format(n)

def comparisonPlot(xStr,yStr,titleStr,xLimit,rawValues,indices,
                   compLabel,mSource):
    # we will plot the 'raw' values, plus whatever is in indices
    numPlots = len(indices)+1
    # first plot: raw times
    fig = pFigure()
    plotCount = 1
    ax = fig.add_subplot(numPlots,1,plotCount)
    rawValues = np.array(rawValues)
    normStrX = 'Normalized ' + xStr
    absStrX = 'Absolute ' + xStr
    normRaw,minV,maxV = histNormalize(rawValues)
    bins = normalizedBins(normRaw)
    # plot the normalizd number of frames, and add on an x label
    # for the absolute number
    maxCounts = len(normRaw)
    axTmp = histogramPlot(ax,normStrX,yStr,
                          titleStr + compLabel[1],normRaw,
                          bins,True,maxCounts)
    plotCount += 1
    secondAxis(axTmp,absStrX,xLimit,False)
    for i in range(numPlots-1):
        # now we plot the values are their given indices
        tmpHist = rawValues[indices[i]]        
        # ignore the 'new' min an max, plot according to the 'raw' values.
        # this ensures that our x axes are all the same
        normHist,junk1,junk2 = histNormalize(tmpHist,minV,maxV)
        # second plot: 'good' times
        ax = fig.add_subplot(numPlots,1,plotCount)
        axTmp = histogramPlot(ax,normStrX,yStr,
                              titleStr + compLabel[i+1],
                              normHist,bins)
        plotCount += 1
    # XXX need to save by other file name, also ...
    fileStr = 'Histogram' + titleStr
    saveFigure(mSource,fileStr,fig)
