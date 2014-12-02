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
    # return the 'normalized' bins or bins numbers for histogramming
    return 100

def denormalize(absArr,toDenorm):
    ignore,minV,maxV = histNormalize(absArr)
    return (toDenorm/normalizedBins(absArr)*(maxV-minV)) + minV

def normalizeTo(normalArr,toNorm):
    # normalize toNorm to the array 'normalArr', according to histNormalize.
    # good for plotting!
    ignore,minV,maxV = histNormalize(normalArr)
    toRet,ignore,ignore = histNormalize(toNorm,minV,maxV)    
    return toRet

def histNormalize(histIn,minV=-1,maxV=-1):
    # histograms are only truly normal if we use bins of size 0.1
    # so, normalize the values to between 0 and 1, use size of since 0.1
    if (minV < 0):
        minV = np.min(histIn)
    if (maxV < 0):
        maxV = np.max(histIn)
    # normalize the histogram and center between 0 and 1 for normal plotting
    normalizedHist = normalizedBins(histIn) *(histIn - minV)/(maxV-minV)
    return normalizedHist,minV,maxV

def histogramPlot(ax,xlabelStr,ylabelStr,titleStr,data,numBins,
                  secondX = False,Normalize=False):
    # wrapper function, ignore the bins which we don't care about.
    ax,ignore,ignore = histogramPlotFull(ax,xlabelStr,ylabelStr,titleStr,
                                         data,numBins,secondX,Normalize)
    return ax

def histogramPlotFull(ax,xlabelStr,ylabelStr,titleStr,data,numBins,
                  secondX = False,Normalize=False):
    numPoints = len(data)
    t = plt.title(titleStr +  getNStr(numPoints))
    if (Normalize):
        # if we are normalizing, the absolute number is the max of the data...
        absoluteNum = np.max(data)
        data,ignoreMin,ignoreMax = histNormalize(data)
        numBins = normalizedBins(data)
    else:
        absoluteNum = numPoints
    # POST: data and absolute scaling number are set
    if (secondX):
        t.set_y(1.4)
        xlim = [0,absoluteNum]
        secondAxis(ax,'Absolute' + xlabelStr,xlim,False)
    n, bins, patches =ax.hist( data, numBins, normed=True,stacked=True,
                                facecolor='green', alpha=0.75)

    ax.set_xlabel(xlabelStr)
    ax.set_ylabel('Norm ' + ylabelStr)
    ax.set_yscale('log', nonposy='clip')
    # limits are half of an individual to slightly more than the entire population
    yLimNorm = [0.5/numPoints,2]
    yLimits= plt.ylim(yLimNorm)
    ax.grid(True,which='major',color='b')
    ax.grid(True,which='minor',color='r')
    # create a twin axis for the absolute count 
    yLimAbs = np.multiply(yLimNorm,numPoints)
    secondAxis(ax,"Abs " + ylabelStr ,yLimAbs)
    # return the *absolute* number in the bins
    return ax,n,bins


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

def saveFigure(figure,fileName):
    # source : where to save the output iunder the output folder
    # filename: what to save the file as. automagically saved as high res pdf
    step = util.globalIO.st
    fullPath = util.globalIO.getOutputDir([step],fileName)
    plt.tight_layout(False)
    formatStr = "png"
    figure.savefig(fullPath + '.' + formatStr,format=formatStr, 
                   dpi=figure.get_dpi())
    plt.close(figure)

def pFigure(ySize=10,xSize=8,dpiArg=100):
    return  plt.figure(figsize=(ySize,xSize),dpi=dpiArg)

def getNStr(n,space = " "):
    return space + "n={:d}".format(n)


def compareHist(xStr,yStr,titleStr,xLimit,matrices,
                   compLabel):
    # we will plot the 'raw' values, plus whatever is in indices
    numPlots = len(matrices)
    # first plot: raw times
    fig = pFigure()
    plotCount = 1
    ax = fig.add_subplot(numPlots,1,plotCount)
    rawValues = np.array(matrices[0])
    normRaw,minV,maxV = histNormalize(rawValues)
    bins = normalizedBins(normRaw)
    # plot the normalizd number of frames, and add on an x label
    # for the absolute number
    axTmp = histogramPlot(ax,xStr,yStr,
                          titleStr + compLabel[0],normRaw,
                          bins,True)
    xlimits = axTmp.get_xlim()
    plotCount += 1
    for i in range(numPlots-1):
        # now we plot the values are their given indices
        tmpHist = matrices[i+1]       
        # ignore the 'new' min an max, plot according to the 'raw' values.
        # this ensures that our x axes are all the same
        normHist,junk1,junk2 = histNormalize(tmpHist,minV,maxV)
        # get the number of bins such that each are of length 1
        numBins = round(max(normHist))
        # second plot: 'good' times
        ax = fig.add_subplot(numPlots,1,plotCount)
        axTmp = histogramPlot(ax,xStr,yStr,
                              titleStr + compLabel[i+1],
                              normHist,numBins)
        axTmp.set_xlim(xlimits)
        plotCount += 1
    # XXX need to save by other file name, also ...
    fileStr = 'Histogram' + titleStr
    saveFigure(fig,fileStr)


def comparisonPlot(xStr,yStr,titleStr,xLimit,rawValues,indices,
                   compLabel):
    # we will plot the 'raw' values, plus whatever is in indices
    rawValues = np.array(rawValues)
    toCompare = [rawValues]
    # number of plots (besides the raw)...
    numPlots = len(indices)
    for i in range(numPlots):
        # now we plot the values are their given indices
        toCompare.append(rawValues[indices[i]])
    # call compare hist with the 'subset' we want
    compareHist(xStr,yStr,titleStr,xLimit,toCompare,compLabel)
        
