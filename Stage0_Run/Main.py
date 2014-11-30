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
import numpy as np

# XXX below lines test this functioon
testDir="../Data/ObjTr/"
#  get a list of files in the test directory
# sort alphabetically
filesFound=sorted(os.listdir(testDir))
testDataType=".dat"
# append path to files, since os.listdir just returns file name
filesFound=[ testDir+f for f in filesFound if str(f).endswith(testDataType)]
# just use a single file for this testing bit.
filesFound = [filesFound[1]]
# get the X,Y velocity, times, and ratio on a per protein basis
# each of these is a list. Each element in a list corresponds to data
# for a single protein
trackedTimes,trackedFRET,trackedDiffusion = GetTraces.GetTracesMain(filesFound)
# XXX functionalize below.
#stage1Folder = "../output/Post_Stage1_Analysis/"
#timeFile = stage1Folder + "Per_Protein_Frames_Appearing.npy"
#diffFile = stage1Folder + "Per_Protein_Diff_Coeff.npy"
#fretFile = stage1Folder + "Per_Protein_Fret_Ratio.npy"
#trackedTimes= np.load(timeFile)
#trackedFRET= np.load(fretFile)
#trackedDiffusion = np.load(diffFile)
unfoldingTimes, diffusionCoeffs = \
                                  GetPhysics.GetPhysicsMain(trackedTimes
                                                            ,trackedFRET,
                                                            trackedDiffusion)
