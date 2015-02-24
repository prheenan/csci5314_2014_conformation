import sys
sys.path.append('../Utilities')
# import model equations from modUtilities
import ModUtilities as modU
# use numpy for any additional calculations
import numpy as np
# import curve_fit for calculations
from scipy.optimize import curve_fit

def getModelVariables(residenceTimes, probability, eqNumber,
                      paramGuess = [1,0.5,1.5],finalParameters=None):
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
    if (finalParameters is not None):
        predicted,stdev = finalParameters(predicted,stdev)
    return predicted,fitParams,(stdev),rsq
    
    
