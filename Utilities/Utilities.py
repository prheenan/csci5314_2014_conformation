# use OS for file IO and the like
import os
# use matplotlib for plotting
import matplotlib.pyplot  as plt
# use numpy for array operations
import numpy as np

def histogramPlot(ax,xlabelStr,ylabelStr,titleStr,data,numBins,
                  onlyPositive=False):
    n, bins, patches =\
                       ax.hist( data, numBins, normed=True,stacked=True,
                                facecolor='green', alpha=0.75)
    plt.title(titleStr)
    plt.ylabel('Norm ' + ylabelStr)
    plt.xlabel(xlabelStr)
    plt.yscale('log', nonposy='clip')
    absoluteNum = len(data)
    if (onlyPositive):
        plt.xlim(left=0)
    yLimNorm = [0.5/absoluteNum,1.05]
    yLimits= plt.ylim(yLimNorm)
    plt.grid(True,which='major',color='b')
    plt.grid(True,which='minor',color='r')
    # create a twin axis for the absolute count 
    yLimAbs = np.multiply(yLimNorm,absoluteNum)
    secondAxis(ax,"Abs " + ylabelStr ,yLimAbs)
    return ax


def secondAxis(ax,label,limits):
    #copies the first axis (ax) and uses it to overlay an axis of limits
    ax2 = ax.twinx()
    # Create a dummy plot
    ax2.set_yscale(ax.get_yscale(), nonposy='clip')
    ax2.set_ylim(limits)
    ax2.set_ylabel(label)

def ReportError(isFatal=False, description="None Given",source="None Given"):
    ReportMessage("Error [" + description +"]",source)
    if (isFatal):
        print("\tError was fatal. Exiting.")
        exit(-1)

def ReportMessage(description="None Given",source="None Given"):
    print("Message [" + str(description) + "] Received by [" 
          + str(source) + "]")


def ensureDirExists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def saveFigure(source,fileName,figure):
    # source : where to save the output iunder the output folder
    # filename: what to save the file as. automagically saved as high res pdf
    defaultOutputDir ="../output/"
    ensureDirExists(defaultOutputDir)
    sourceOutputDir = defaultOutputDir + source +'/'
    ensureDirExists(sourceOutputDir)
    fullPath = sourceOutputDir + fileName
    figure.savefig(fullPath + '.png', dpi=500)

    
