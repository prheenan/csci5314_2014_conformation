# include os for file IO
import os
# include sys to add things to the path
import sys
# more or less a way of 'C-style' importing,  to be able to run wherever
foldersToImport = ['Stage0_Run','Stage1_DataToTraces','Stage2_TracesToPhyscs',
                   'Stage3_PhysicsToModel','Testing','Utilities']
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

# XXX below lines test this functioon
testDir="../Data/ObjTr/"
#  get a list of files in the test directory
# sort alphabetically
filesFound=sorted(os.listdir(testDir))
testDataType=".dat"
# append path to files, since os.listdir just returns file name

trials = ["5x_PBS_25C_","AP2_minus6_AHAPS","AP2_minus6_pbs","AP2minus4_FS"]
numTrials = len(trials)
# mean, stdev, and RSQ for the first tau (for now
statsLabel = ['Tau0','Tau1','1-RSQ']
# each tau value (everything except RSQ) has an uncertainty
numBars = len(statsLabel)
numParams = (numBars-1)*2+1

trialStats = np.zeros((numTrials,numParams))
# mean will just have error bars...

# file names for stage 1, defined in the common utility class.
fileNamesStage1 = [Utilities.IO_frames,Utilities.IO_diff,Utilities.IO_fret]
fileNamesStage2 = [Utilities.IO_frames,Utilities.IO_diff,Utilities.IO_indices]
fileNamesStage3 = [Utilities.IO_model,Utilities.IO_model_stdevs,Utilities.IO_model_rsq]

# can 'force' a re-calulation of the data
forceStage1 = False
forceStage2 = False
forceStage3 = False
# the numbers in the modeling file, determining which models we will use..
modelsToUse = [1,2,3]
# look for 'cached' files, which allow us to save the previous stages, speeding up everything a lot.
extFile = ".npy"

# loop through each trial separately, averaging the independent
# runs
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
    stage3Folder = Utilities.globalIO.getOutputDir([],"") + Utilities.IO_Stage3Folder
    # stage 3 recquires knowing which model we will use...
    # XXX fix this to allow for different models? Pass as a param to model.
    cacheFileS3 = stage3Folder + Utilities.getModelFile(fileNamesStage3[0],1) \
                  + extFile
    # if we have stage3, then just skip to the final ....
    if (not Utilities.dirExists(cacheFileS3)
        or forceStage3):
        for f in trialFiles:
            fileNameExt = os.path.basename(f)
            fileName = os.path.splitext(fileNameExt)[0]
            Utilities.globalIO.setFile(fileName)
            Utilities.globalIO.setStep("Step1::GetTraces")

            # get the output dir with no prior / after path
            outputDir = Utilities.globalIO.getOutputDir([],"")
            stage1Folder = outputDir + Utilities.IO_Stage1Folder
            stage2Folder = outputDir + Utilities.IO_Stage2Folder

            cacheFileS1 = stage1Folder + fileNamesStage1[0] + extFile
            cacheFileS2 = stage2Folder + fileNamesStage2[0] + extFile

            # Next, start stage 1 if we need it; otherwise just pull in the needed 
            # variables.
            if (not Utilities.dirExists(cacheFileS1)
                or forceStage1):
                # get the X,Y velocity, times, and ratio on a per protein basis
                # each of these is a list. Each element in a list corresponds to data
                # for a single protein
                trackedTimes,trackedFRET,trackedDiffusion = \
                                                            GetTraces.GetTracesMain(f)
                # get the output dir for the previous path. 
            else:
                # the output directory exists; no sense in re-running the analysis.
                # just use the text output.
                Utilities.ReportMessage("Skipping Stage1 since [" + cacheFileS1 + 
                                        "] already exists")
                trackedTimes, trackedDiffusion,trackedFRET = \
                                            Utilities.loadAll(stage1Folder,fileNamesStage1)
                # POST: now have valid, trackable diffusion, times, and FRET
                # Move to stage 2. where we get the unfolding times
                Utilities.globalIO.setStep("Step2::GetPhysics")
                if (not Utilities.dirExists(cacheFileS2)
                    or forceStage2):
                    unfoldingTimes, diffusionCoeffs,trackedIndices = \
                                        GetPhysics.GetPhysicsMain(trackedTimes
                                                        ,trackedFRET,trackedDiffusion)
                else:
                    Utilities.ReportMessage("Skipping Stage2 since [" + cacheFileS2 + 
                                            "] already exists")
                    unfoldingTimes, diffusionCoeffs,trackedIndices = \
                                                Utilities.loadAll(stage2Folder,fileNamesStage2)
                    # XXX fix ths flattening
            trialTimes.append(unfoldingTimes)
            trialDiff.append(diffusionCoeffs)
        # POST: valid unfolding time and diffusion coefficients.
        # same the unfolding time and diffusion coefficients across all the trial
        # Move to Stage 3, where we get the residence time distribution
        Utilities.globalIO.setStep("Step3::GetModel")
        # no longer in any particular file, so save everything in the 'top level'
        Utilities.globalIO.setFile("")

        params,stdev,RSQ = GetModel.GetModelMain(trialTimes,trialDiff,modelsToUse)
    else:
        # stage 3 files exist. Just need to load the files to be able to compare...
        params = []
        stdev = []
        RSQ = []
        for num in modelsToUse:
            thisModelFiles = [Utilities.getModelFile(f,num) for f in fileNamesStage3]
            paramsTmp, stdevTmp,RSQTmp = \
                                Utilities.loadAll(stage3Folder,thisModelFiles)
            params.append(paramsTmp)
            stdev.append(stdevTmp)
            RSQ.append(RSQTmp)
        Utilities.ReportMessage("Skipping Stage3 since [" + cacheFileS3 + 
                                "] already exists, plotting comparison.")
    # for each trial, for now just look at a single Model
    modelToCheck = 1
    fitStdRsq = [params[modelToCheck][-2:],stdev[modelToCheck][-2:],RSQ[modelToCheck]]
    stats = np.concatenate(fitStdRsq)
    trialStats[tNum,:] = stats

# POST: all trials completed.
Utilities.globalIO.reset()
outputDir = Utilities.globalIO.getOutputDir([],"") 
f = pltUtil.pFigure()
tickLabelX = np.arange(0,numBars)
colorCycler = pltUtil.cycleColors()
ax = plt.subplot(1,1,1)
axRSQ = pltUtil.secondAxis(ax,'1-RSQ',[0,0.002])
wid = 1/(numTrials*numBars)
plt.title('Comparison of Surface Modeling Parameters')
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
           color=mColor,ecolor='m',align='center',label=mLabel)
    axRSQ.bar(xVals[numVals],1-RSQ,width=wid,color=mColor,align='center',label=mLabel)
    # only the second time, add a third axis
#axRSQ.set_yscale('log')
ax.legend(loc='upper left')
ax.set_ylabel('Tau value (seconds)')
plt.xlim([-0.5,max(xVals)*1.1])
plt.xticks(tickLabelX, statsLabel)
pltUtil.saveFigure(f,"Comparison")







