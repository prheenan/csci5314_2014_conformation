import os

import matplotlib.pyplot  as plt

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

    
