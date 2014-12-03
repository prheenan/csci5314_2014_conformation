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
import numpy as np

# XXX below lines test this functioon
testDir="../Data/ObjTr/"
#  get a list of files in the test directory
# sort alphabetically
filesFound=sorted(os.listdir(testDir))
testDataType=".dat"
# append path to files, since os.listdir just returns file name

trials = ["5x_PBS_25C_","AP2_minus6_AHAPS","AP2_minus6_pbs","AP2minus4_FS",
          "AP2minus4_OEG"]

# file names for stage 1, defined in the common utility class.
fileNamesStage1 = [Utilities.IO_frames,Utilities.IO_diff,Utilities.IO_fret]
fileNamesStage2 = [Utilities.IO_frames,Utilities.IO_diff,Utilities.IO_indices]
# can 'force' a re-calulation of the data
forceStage1 = False
forceStage2 = False

# loop through each trial separately, averaging the independent
# runs
filesFound=[ testDir+f for f in filesFound if str(f).endswith(testDataType)]
if (len(filesFound) < 1):
    msg = "Couldn't find any file like [{}] in [{}]".format(testDataType,
                                                            testDir)
    Utilities.ReportError(True,msg,"Main.py")

for t in trials:
    # get the matching files for this trial
    trialFiles = [ fileN for fileN in filesFound if t in fileN]
    trialTimes = []
    trialDiff = []
    for f in trialFiles:
        fileNameExt = os.path.basename(f)
        fileName = os.path.splitext(fileNameExt)[0]
        Utilities.globalIO.setTrial(t)
        Utilities.globalIO.setFile(fileName)
        Utilities.globalIO.setStep("Step1::GetTraces")
        outputDir = Utilities.globalIO.getOutputDir([],"")
        stage1Folder = outputDir + Utilities.IO_Stage1Folder
        stage2Folder = outputDir + Utilities.IO_Stage2Folder
        extFile = ".npy"
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
        unfoldingTimes = np.ravel(unfoldingTimes)
        diffusionCoeffs = np.ravel(diffusionCoeffs)
        trialTimes.append(unfoldingTimes)
        trialDiff.append(diffusionCoeffs)
        # POST: valid unfolding time and diffusion coefficients.
        # same the unfolding time and diffusion coefficients across all the trial
    # Move to Stage 3, where we get the residence time distribution
    Utilities.globalIO.setStep("Step3::GetModel")
    # no longer in any particular file, so save everything in the 'top level'
    Utilities.globalIO.setFile("")

    GetModel.GetModelMain(np.concatenate(trialTimes),np.concatenate(trialDiff))
    exit(1)
