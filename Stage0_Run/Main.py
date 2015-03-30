# include os for file IO
import os
# include sys to add things to the path
import sys
# more or less a way of 'C-style' importing,  to be able to run wherever
# XXX TODO : should fix this out here
import argparse


foldersToImport = ['Stage0_Run','Stage1_DataToTraces',
                   'Stage2_TracesToPhyscs','Stage3_PhysicsToModel',
                   'Testing','Utilities','Base']
for f in foldersToImport:
    # do a relative and absoute add.
    tmpF = '/' + f + '/'
    sys.path.append('..' + tmpF)
    sys.path.append('.' + tmpF)
    
# import the files needed
# XXX will need to import your files when you need them (whatever the 
# python files are named, just import that: "GetTraces.py" --> GetTraces)
import GetTraces
import GetPhysics
import GetModel
import Utilities
import PlotUtilities as pltUtil
import matplotlib.pyplot as plt
import numpy as np
import Filter

def plotStats(numBars,numTrials,trialStats):
    f = pltUtil.pFigure()
    tickLabelX = np.arange(0,numBars)
    colorCycler = pltUtil.cycleColors()
    ax = plt.subplot(1,1,1)
    axRSQ = pltUtil.secondAxis(ax,'RSQ',[0,1.0])
    wid = 1/(numTrials*numBars)
    plt.title('Comparison of Surface Modeling Parameters')
    barDict = dict(elinewidth=4, ecolor='m')
    for t in range(numTrials):
        # XXX fix this to make it more general...
        RSQ = trialStats[t,-1]
        numMeanStd = numParams-1
        middle = int(numMeanStd/2)
        meanVals = trialStats[t,0:middle]
        stdVals = trialStats[t,middle:numMeanStd] 
        numVals = len(meanVals)
        xVals = tickLabelX+wid*t
        mColor = next(colorCycler)
        mLabel = trials[t]
        ax.bar(xVals[0:numVals],meanVals,width=wid,yerr=0.2*stdVals,
               color=mColor,error_kw = barDict,align='center',label=mLabel)
        axRSQ.bar(xVals[numVals],RSQ,width=wid,color=mColor,align='center',
                  label=mLabel)
        # only the second time, add a third axis
    ax.legend(loc='upper left')
    ax.set_ylabel('Tau value (seconds)')
    plt.xlim([-0.5,max(xVals)*1.1])
    plt.xticks(tickLabelX, statsLabel)
    pltUtil.saveFigure(f,"Comparison")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Protein Data Analysis Args')
    parser.add_argument('--inPath', type=str, default="../Data/ObjTr/",
                        help="Folder where formatted .dat file reside")
    parser.add_argument('--outPath', type=str, default="../output/",
                        help="Default output directory (base)")
    args = parser.parse_args()
    # XXX below lines test this functioon
    testDir=args.inPath
    Utilities.globalIO.globalOutput = args.outPath
    #  get a list of files in the test directory
    # sort alphabetically
    filesFound=sorted(os.listdir(testDir))
    testDataType=".dat"
    # append path to files, since os.listdir just returns file name
    trials = ["5x_PBS_25C_","AP2_minus6_AHAPS","AP2_minus6_pbs","AP2minus4_FS"]
    numTrials = len(trials)
    # mean, stdev, and RSQ for the first tau (for now
    statsLabel = ['A[seconds]','B[au]','RSQ']
    # each tau value (everything except RSQ) has an uncertainty
    numBars = len(statsLabel)
    numParams = (numBars-1)*2+1
    trialStats = np.zeros((numTrials,numParams))
    # mean will just have error bars...
    fileNamesStage3 = [Utilities.IO_model,Utilities.IO_model_stdevs,
                       Utilities.IO_model_rsq]
    # can 'force' a re-calulation of the data
    forceStage1 = False
    forceStage2 = False
    forceStage3 = False
    # the numbers in the modeling file, determining which models we will use..
    modelsToUse = [4]
    # look for 'cached' files, by previous stages, speeding up everything a lot.
    extFile = ".npy"
    # loop through each trial separately, averaging the independent runs
    filesFound=[ testDir+f for f in filesFound if str(f).endswith(testDataType)]
    if (len(filesFound) < 1):
        msg = "Couldn't find any file like [{}] in [{}]".format(testDataType,
                                                                testDir)
        Utilities.ReportError(True,msg,"Main.py")

    for tNum,t in enumerate(trials):
        # get the matching files for this trial
        trialFiles = [ fileN for fileN in filesFound if t in fileN]
        trialTimes = []
        trialDiff = []
        Utilities.globalIO.setTrial(t)
        stage3Folder = (Utilities.globalIO.getOutputDir([],"") +
                        Utilities.CheckpointDir)
        # stage 3 recquires knowing which model we will use...
        # XXX fix this to allow for different models? Pass as a param to model.
        cacheFileS3 = stage3Folder + Utilities.getModelFile(fileNamesStage3[0],
                                                            1) + extFile
        # if we have stage3, then just skip to the final ....
        if (not Utilities.dirExists(cacheFileS3)
            or forceStage3):
            fileData = []
            for f in trialFiles:
                fileNameExt = os.path.basename(f)
                fileName = os.path.splitext(fileNameExt)[0]
                Utilities.globalIO.setFile(fileName)
                Utilities.globalIO.setStep("Step1::GetTraces")

                # get the output dir with no prior / after path
                outputDir = Utilities.globalIO.getOutputDir([],"")
                checkPointDir = outputDir + Utilities.CheckpointDir
                checkPointFilePath = checkPointDir + 'Stage_'
                # ensure the saving directory exists
                Utilities.ensureDirExists(checkPointDir)
                stageCount = 1
                # Next, start stage 1 if we need it; 
                #otherwise just pull in the needed variables.
                # get the X,Y velocity, times, and ratio on a per protein basis
                # each of these is a list. Each element  corresponds to data
                # for a single protein
                # save the raw data as well
                rawFile = checkPointFilePath + str(stageCount)
                stageCount += 1
                traceFilter = GetTraces.GetTraces(f,checkPointFilePath + 
                                                  str(stageCount),rawFile,
                                                  forceStage1)
                filteredByTrace = traceFilter.filterAndCheckPoint()
                # get the output dir for the previous path. 
                # Move to stage 2. where we get the unfolding times
                stageCount += 1
                Utilities.globalIO.setStep("Step2::GetPhysics")
                physFilter = GetPhysics.GetPhysics(filteredByTrace,
                                                   checkPointFilePath +
                                                   str(stageCount),
                                                   forceStage2)
                filteredByPhysics = physFilter.filterAndCheckPoint()
                fileData.append(filteredByPhysics)
            # POST: valid unfolding time and diffusion coefficients.
            # same the unfolding time and diffusion coefficients across trial
            # Move to Stage 3, where we get the residence time distribution
            Utilities.globalIO.setStep("Step3::GetModel")
            # no longer in any particular file, so save in 'top level'
            Utilities.globalIO.setFile("")
            params,stdev,RSQ = GetModel.GetModelMain(fileData,modelsToUse)
        else:
            # stage 3 files exist. Just need to load the files to be able 
            # to compare...
            params = []
            stdev = []
            RSQ = []
            for num in modelsToUse:
                thisModelFiles = [Utilities.getModelFile(f,num) 
                                  for f in fileNamesStage3]
                paramsTmp, stdevTmp,RSQTmp = \
                                    Utilities.loadAll(stage3Folder,
                                                      thisModelFiles)
                params.append(paramsTmp)
                stdev.append(stdevTmp)
                RSQ.append(RSQTmp)
            Utilities.ReportMessage("Skipping Stage3 since [" + cacheFileS3 + 
                                    "] already exists, plotting comparison.")
        # for each trial, for now just look at a single Model
        # XXX fix this!...    
        modelToCheck = 0
        fitStdRsq = [params[modelToCheck][:],stdev[modelToCheck][:],
                     RSQ[modelToCheck]]
        stats = np.concatenate(fitStdRsq)
        trialStats[tNum,:] = stats
    # POST: all trials completed.
    Utilities.globalIO.reset()
    outputDir = Utilities.globalIO.getOutputDir([],"") 
    plotStats(numBars,numTrials,trialStats)







