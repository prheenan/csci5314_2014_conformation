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

dataDir = "../output/"
mWorking = "./vizTmp/"
mOut = "./vizOut/"
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
        X,Y = getAllStages(fileDict,key,keyV,mWorking,i,j)
        saveAsSubplot(X,Y,allStageDir)

