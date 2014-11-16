# include os for file IO
import os
# include sys to add things to the path
import sys
sys.path.append('../Stage1_DataToTraces/')
# import the traces file
import GetTraces

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
GetTraces.GetTracesMain(filesFound)
