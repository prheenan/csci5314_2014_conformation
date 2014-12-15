import sys
sys.path.append('../Utilities')
# import model equations from modUtilities
import ModUtilities as modU
# use numpy for any additional calculations
import numpy as np
# import curve_fit for calculations
from scipy.optimize import curve_fit

def getModelVariables(residenceTimes, probability, eqNumber,
                      paramGuess = [1,0.5,1.5]):
    #Use ModUtilities to get model equation 'eqNumber' to fit
    # x values 'residence times' against y values 'probability'
    #Initial guess for variables A, f, tau in that order is given by paramGuess 
    model = modU.modEquations(eqNumber)
    xdata = residenceTimes
    ydata = probability
    #Use curve_fit with given times and probability on the given eq. Number
    fitParams, fitCov = curve_fit(model, xdata,ydata, p0=paramGuess)
    stdev = np.sqrt(np.diag(fitCov))
    predicted = model(xdata,*fitParams)
    rsq = modU.RSQ(predicted,probability)
    if (len(fitParams) > 1):
        # parameters go like [f1,f2,...,tau1,tau2]
        numParams = len(fitParams)
        numFits = numParams/2.
        numfVals = np.floor(numFits)
        fIndices = np.arange(0,numfVals,dtype=np.int32)
        fVals = fitParams[fIndices]
        # get the last f based on the models...
        fStdevs = fitParams[fIndices]
        lastF = 1 - np.sum(fVals)
        lastFStdev = np.sqrt(np.sum(fStdevs**2))
        #
        fitParams = np.insert(fitParams,[numfVals],[lastF])
        stdev = np.insert(stdev,[numfVals],[lastFStdev])
    else:
        fitParams = np.insert(fitParams,[0],[1])
        stdev = np.insert(stdev,[0],[0])
    return predicted,fitParams,(stdev),rsq
    
    
