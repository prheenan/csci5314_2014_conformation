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
forceStage2 = True

# loop through each trial separately, averaging the independent
# runs
for t in trials:
    # get the matching files
    filesFound=[ testDir+f for f in filesFound if str(f).endswith(testDataType)]
    if (len(filesFound) < 1):
        msg = "Couldn't find any file like [{}] in [{}]".format(testDataType,
                                                                testDir)
        Utilities.ReportError(True,msg,"Main.py")
    # XXX
    # XXX need to check for the given files...
    # XXX 
    # just use a single file for this testing bit.
    filesFound = filesFound[0]
    fileNameExt = os.path.basename(filesFound)
    fileName = os.path.splitext(fileNameExt)[0]
    Utilities.globalIO.setTrial(t)
    Utilities.globalIO.setFile(fileName)
    Utilities.globalIO.setStep("Step1::GetTraces")
    outputDir = Utilities.globalIO.getOutputDir([],"")
    stage1Folder = outputDir + Utilities.IO_Stage1Folder
    stage2Folder = outputDir + Utilities.IO_Stage2Folder
    # Next, start stage 1 if we need it; otherwise just pull in the needed 
    # variables.
    if (not Utilities.dirExists(stage1Folder) or forceStage1):
        # get the X,Y velocity, times, and ratio on a per protein basis
        # each of these is a list. Each element in a list corresponds to data
        # for a single protein
        trackedTimes,trackedFRET,trackedDiffusion = \
                GetTraces.GetTracesMain(filesFound)
        # get the output dir for the previous path. 
    else:
        # the output directory exists; no sense in re-running the analysis.
        # just use the text output.
        Utilities.ReportMessage("Skipping Stage1 since [" + stage1Folder + 
                                "] already exists")
        trackedTimes, trackedDiffusion,trackedFRET = \
                        Utilities.loadAll(stage1Folder,fileNamesStage1)
    # POST: now have valid, trackable diffusion, times, and FRET
    # Move to stage 2. where we get the unfolding times
    Utilities.globalIO.setStep("Step2::GetPhysics")
    if (not Utilities.dirExists(stage2Folder) or forceStage2):
        unfoldingTimes, diffusionCoeffs,trackedIndices = \
                    GetPhysics.GetPhysicsMain(trackedTimes
                                              ,trackedFRET,trackedDiffusion)
    else:
        Utilities.ReportMessage("Skipping Stage2 since [" + stage2Folder + 
                                "] already exists")
        unfoldingTimes, diffusionCoeffs,trackedIndices = \
                        Utilities.loadAll(stage2Folder,fileNamesStage2)
    # POST: valid unfolding time and diffusion coefficientrs.
    # Move to Stage 3, where we get the residence time distribution
    Utilities.globalIO.setStep("Step3::GetModel")
    GetModel.GetModelMain(unfoldingTimes,diffusionCoeffs)
    exit(1)
