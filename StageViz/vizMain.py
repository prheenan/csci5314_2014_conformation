import os
from collections import defaultdict
import numpy as np
# need to import the sparse matrices and stacking funciton
from scipy.sparse import csr_matrix,csc_matrix
import sys
sys.path.append('../Utilities')
import CheckpointUtilities as cPoint
import GenUtilities as gUtil
import PlotUtilities as pUtil
import matplotlib.pyplot as plt
from vizUtil import getCheckpointFileDict,getDirs
from vizPlot import saveAsSubplot,getAllStages

import subprocess
dataDir = "../output/"
mWorking = "./vizTmp/"
mOut = "./vizOut/"
# next two must match, for the automatic video encoding to work
vizFileFormat =  "t{:05d}"
ffmpegFormat = "t%05d"
fps = 10 # frames per second, should match video
vizExt = ".png"
gUtil.ensureDirExists(mOut)
gUtil.ensureDirExists(mWorking)
# get all the files. returns a dictionary of dictionaries.
# each key in the outer (first) dictionary is a condition
# each key in the innter (second) dictionary is a trial for that condition
fileDict = getCheckpointFileDict(dataDir)

# loop through each condition and trial
for i,key in enumerate(fileDict.keys()):
    for j,keyV in enumerate(fileDict[key].keys()):
        # for this condition and trial, get the directories and all the 
        # data, then save
        trialDir,allStageDir =  getDirs(mOut,key,keyV,j)
        #X,Y,c1,c2 = getAllStages(fileDict,key,keyV,mWorking,i,j)
        command = 'ffmpeg'
        # format the ffmpeg arguments as we want them
        args = ('-i {:s} -c:v libx264 -r {:d} -y -pix_fmt yuv420p '+
                '{:s}1_movie_{:s}.mov').\
            format(allStageDir+ffmpegFormat+vizExt,fps,allStageDir,key)
        print(args)
        #saveAsSubplot(X,Y,c1,c2,allStageDir,vizFileFormat)
        # POST: all videos saved for this trial. make the movie
        ret = os.system(command + ' ' + args)
        
