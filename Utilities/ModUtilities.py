# use numpy for mathmatical functions
import numpy as np

def modelString(equation,values,stdevs,labels,units,RSQ):
    mStr = "Equation: {:s}, RSQ: {:.3f}".format(equation,RSQ)
    for val,std,label,unit in zip(values,stdevs,labels,units):
        mStr += "\n{:}: {:3.2g} +/- {:2.1g} [{:}]".format(label,val,std,unit)
    return mStr


def modelStrGen(params,paramsStd,predicted,actual,labels,modelString,
                returnParams=False):
    # generic string generation. Useful for labelled plots
    model = '{0:s}, RSQ: {1:.3f}\n'.format(modelString,
                                         RSQ(predicted,actual))
    paramStr = []
    count =0
    for l,p,std in zip(labels,params,paramsStd):
        thisModel = "{0:s}={1:.2f}±{2:.2f}".format(l,p,std)
        paramStr.append(thisModel)
        if (count % 2 == 0 and count > 0 ):
            model += '\n' + thisModel
        elif (count > 0):
            model += ', ' + thisModel
        else:
            model += thisModel
        count += 1
    if (returnParams):
        return model,paramStr
    else:
        return model



def RSQ(yPred,yActual):
    mean = np.mean(yActual)
    SStot = np.sum((yActual-mean)**2)
    SSres = np.sum((yPred-yActual)**2)
    return 1-SSres/SStot

def exponentialGetFinalParams(fitParams,stdev):
    if (len(fitParams) > 1):
        # parameters go like [f1,f2,...,tau1,tau2]
        numParams = len(fitParams)
        numFits = numParams/2.
        numfVals = np.floor(numFits)
        fIndices = np.arange(0,numfVals,dtype=np.int32)
        fVals = fitParams[fIndices]
        # get the last f based on the models...
        fStdevs = stdev[fIndices]
        lastF = 1 - np.sum(fVals)
        lastFStdev = np.sqrt(np.sum(fStdevs**2))
        fitParams = np.insert(fitParams,[numfVals],[lastF])
        stdev = np.insert(stdev,[numfVals],[lastFStdev])
    else:
        fitParams = np.insert(fitParams,[0],[1])
        stdev = np.insert(stdev,[0],[0])
    return fitParams,stdev

def hyperfitLogStr(params,paramsStd,predicted,actual,returnParams=False):
    modelStr = 'b * t/(t+a)'
    labels=['a','b']
    return modelStrGen(params,paramsStd,predicted,actual,labels,modelStr,
                       returnParams)
def hyperfitLog(t,a,b):
    return  (-b*t/(t+a))

def modEqGen(t,fList,tauList):
    # general 'wrapper' for all tau...
    numTimes = np.shape(t)
    toRet = np.zeros(numTimes)
    for tauI,fI in zip(tauList,fList):
        if (tauI < 0 or fI < 0):
            return np.zeros(numTimes)
        toRet += fI*np.exp(-t/tauI)
    if (not np.isfinite(toRet).any()):
        return np.zeros(numTimes)
    return toRet

def modEq1(t, tau):
    #Equation 2 from Walder 2014
    return modEqGen(t,[tau],[1])

def modEqx(t, f, r, D, tau):
    #Equation 6 from Walder 2014
    return np.sum(f * np.exp((-r**2)/(4 * D * t)))


def modEq2(t,f1,tau1,tau2):
    toRet = 0
    tau  = [tau1,tau2]
    f = [f1,1-f1]
    return modEqGen(t,tau,f)

def modEq3(t,f1,f2,tau1,tau2,tau3):
    toRet = 0
    tau  = [tau1,tau2,tau3]
    f = [f1,f2,1-(f1+f2)]
    return modEqGen(t,f,tau)


def modEquations(eqNumber):
    if (eqNumber == 1):
        #Equation 2 from Walder 2014
        model = modEq1
        
    if (eqNumber == 2):
        model = modEq2
    if (eqNumber == 3):
        model = modEq3
    if (eqNumber == 4):
        #Poisson eq
        model = hyperfitLog
        
    #Return the equation wanted as 'model'
    return model
        
