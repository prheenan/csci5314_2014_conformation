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
import Utilities
import numpy as np

# XXX below lines test this functioon
testDir="../Data/ObjTr/"
#  get a list of files in the test directory
# sort alphabetically
filesFound=sorted(os.listdir(testDir))
testDataType=".dat"
# append path to files, since os.listdir just returns file name

trials = ["5x_PBS_25C_","AP2_minus6_AHAPS","AP2_minus6_pbs","AP2minus4_FS","AP2minus4_OEG"]

# file names for stage 1
# XXX probably want to functionalize these?
fileNames = ["Per_Protein_Frames_Appearing.npy",  
             "Per_Protein_Diff_Coeff.npy",
             "Per_Protein_Fret_Ratio.npy"]

for t in trials:
    filesFound=[ testDir+f for f in filesFound if str(f).endswith(testDataType)]
    if (len(filesFound) < 1):
        msg = "Couldn't find any file like [{}] in [{}]".format(testDataType,
                                                                testDir)
        Utilities.ReportError(True,msg,"Main.py")
    # just use a single file for this testing bit.
    filesFound = [filesFound[0]]
    fileNameExt = os.path.basename(filesFound[0])
    fileName = os.path.splitext(fileNameExt)[0]
    Utilities.globalIO.setTrial(t)
    Utilities.globalIO.setFile(fileName)
    Utilities.globalIO.setStep("Step1::GetTraces")
    outputDir = Utilities.globalIO.getOutputDir([],"")
    stage1Folder = outputDir + "Post_Stage1_Analysis/"
    if (not Utilities.dirExists(stage1Folder)):
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
                        Utilities.loadAll(stage1Folder,fileNames)
    # POST: now have valid, trackable diffusion, times, and FRET
    # Move to stage 2. where we get the unfolding times
    Utilities.globalIO.setStep("Step2::GetPhysics")
    unfoldingTimes, diffusionCoeffs = \
                    GetPhysics.GetPhysicsMain(trackedTimes
                                              ,trackedFRET,trackedDiffusion)
    exit(1)
