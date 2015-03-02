import abc
import copy
import numpy as np
import sys
sys.path.append('../Utilities')
import Utilities as Util
from os.path import expanduser
home = expanduser("~")
# get the utilties directory (assume it lives in ~/utilities/python)
# but simple to change
path= home +"/utilities/python"
import sys
sys.path.append(path)
import GenUtilities as pGenUtil

prop_time = 'times'
prop_c1 = 'FRET_C1'
prop_c2 = 'FRET_C2'
prop_ratio = 'FRET_Ratio'
prop_X = 'X'
prop_Y = 'Y'
prop_vx = 'velX'
prop_vy = 'velY'
prop_MSD= 'MSD'
prop_diff = 'diffusion'
prop_resi = 'resTimes'
prop_model = 'model'
prop_prob = 'probabilities'
prob_idx = 'idx'

class CheckpointData:
    def _sanit(self,listV):
        if (listV) is None:
            return None
        else:
            return np.array(listV,dtype=np.object)
    def __init__(self,times,FRET_C1,FRET_C2,FRET_Ratio,X,Y,velX,velY,MSD,
                 diffusion=None,resTimes = None,model=None,probabilities=None,idx=None):
        # save everything
        self._t = self._sanit(times)
        self._c1 =self._sanit( FRET_C1)
        self._c2 = self._sanit(FRET_C2)
        self._ratio = self._sanit(FRET_Ratio)
        self._x = self._sanit(X)
        self._y = self._sanit(Y)
        self._velX = self._sanit(velX)
        self._velY = self._sanit(velY)
        self._MSD = self._sanit(MSD)
        self._diff = self._sanit(diffusion)
        self._resTime = self._sanit(resTimes)
        self._model = self._sanit(model)
        self._prob = self._sanit(probabilities)
        if (idx is None):
            idx = self._sanit(range(len(times)))
        self._indices = idx
        # XXX use names a
        self._propList = \
            {prop_time: self._t,
            prop_c1:self._c1,
            prop_c2:self._c2,
            prop_ratio:self._ratio,
            prop_X:self._x,
            prop_Y:self._y,
            prop_vx:self._velX,
            prop_vy:self._velY,
            prop_MSD:self._MSD,
            prop_diff:self._diff,
            prop_resi:self._resTime,
            prop_model:self._model,
            prop_prob:self._prob,
            prob_idx: self._indices
            }
    def filterBy(self,idx):
    # filter all properties by the index
    # lists are pass be reference in python, so we can just index into the
    # property list. Assumes the index is OK
        numProps = len(self._propList)
        keys = self._propList.keys()
        for i in range(numProps):
            thisKey = keys[i]
            if (self._propList[thisKey] is not None):
                self._propList[thisKey] = self._propList[thisKey][idx]
            else:
                # nothing to do! nothing to index...
                pass
        # return a deep copy of the new list, so we can pass it to the next 
        # and avoid corrupting this checkpoint
        return copy.deepcopy(self)
        # POST: filtered everything.
    def getData(self):
        return self._propList
    def getSingleVal(self,name):
        return self._propList[name]
    def getValsAt(self,nameList):
        return [self.getSingleVal(name) for name in nameList]
    def __getitem__(self,name):
        if (name not in self._propList):
            Util.ReportError(True,
            description=(("Variable {:s} not found in keys {:s}")
                         .format(str(name),str(self._propList.keys()))),
                        source='CheckPointData :: __getitem__')
        return self._propList[name]

class DataFilter:
    __metaclass__ = abc.ABCMeta
    def __init__(self,data,mFile,force=False,frameRate=0.1,ext='.npz'):
        ''' passed in a data object ''' 
        self._data = data
        self._filePath = pGenUtil.ensureEnds(mFile,ext)
        self._force = force
        self._frameRate = frameRate
    @abc.abstractmethod
    def getFilterIdx(self):
        ''' child class must implement this ''' 
        # given the current data, returns the indices of the 'good' information
        # if necessary, also sets properties in the data object
        # for example, the GetModel function will not have the probabilities
        # set. 
        return 
    @classmethod
    def callIfNoFile(cls,toCall,fileN):
        if (pGenUtil.isfile(fileN) ):
            return np.load(fileN)
        else:
            return toCall(fileN)
    def filter(self):
        force = self._force
        if (pGenUtil.isfile(self._filePath) and not force):
            # return a new checkpoint object from the data
            data = np.load(self._filePath)
            return CheckpointData(**data)
        else:
            idx = self.getFilterIdx()
            return DataFilter.filterDataStatic(self._data,idx,
                                               self._filePath,force)
    def filterAndCheckPoint(self):
        # force is true/false if we want to manually reload.
        newData = self.filter()
        # save all the data we have
        DataFilter.checkPoint(newData,self._filePath)
        return newData
    @classmethod 
    def filterDataStatic(cls,data,idx,output,force=False):
        # POST: need to filter all the input data by the index.
        return data.filterBy(idx)
        # done!
    @classmethod 
    def filterAndCheckStatic(cls,data,idx,output,force=False):
        Filtered = filterDataStatic(cls,data,idx,output,force=False)
        DataFilter.checkPoint(filtered,output)
        # done!
    @classmethod 
    def checkPoint(cls,data,output):
        toSave = dict()
        for key,val in data.getData().iteritems():
            if (val is not None):
                toSave[key] = val
            else:
                # dont save none variables
                pass
        np.savez_compressed(output,**toSave)
    @classmethod
    def concatenate(cls,dataObjects):
        # dataobjects is a list of objects with the same not-non propteries
        # name of the properties
        properties = dataObjects[0].getData().keys()
        # get the last index, then take the max accross all videos
        # allows us to maintain unique indices later...
        maxIdx = np.max([ dataObj[prob_idx][-1] 
                                for dataObj in dataObjects])
        concatProps = dict()
        # we assume all properties are either set or not -- nothing in between
        for prop in properties:
            # need to treat ID a little special, since we want to preserve
            # unique ID's accross all the data objects. Idx is just an internal
            # reference..
            if (dataObjects[0][prop] is None):
                concatProps[prop] = None
            else:
                # must all be defined
                if (prop == prob_idx):
                    arr= [dataObj[prop]+maxIdx*i 
                          for i,dataObj in enumerate(dataObjects) ]
                else:
                    arr = [dataObj[prop] for dataObj in dataObjects]
                    # POST: R is set with the properties we want, now array-ify
                    # and flatten. lose the video-specific information
                concatProps[prop] = np.array(arr,dtype=np.object).flatten()
        return CheckpointData(**concatProps)
    def get(self,nameList):
        return self._data.getValsAt(nameList)
    def _isNoneOrError(self,toCheck):
        if (toCheck is not None):
            Util.ReportError(True,description="Variable to check was none",
                             source='Filter::_isNoneOrError')
        return
    def setData(self,newData):
        self._isNoneOrError(self._data)
        # POST: data was none
        self._data = newData
    def setElement(self,name,value):
        self._isNoneOrError(self._data._propList[name])
        # POST: wasnt set
        self._data._propList[name] = value


            
        
