# import model equations from modUtilities
import ModUtilities as modU
# use numpy for any additional calculations
import numpy as np
# import curve_fit for calculations
from scipy.optimize import curve_fit

def getModelVariables(residenceTimes, probability, eqNumber):
    #Use ModUtilities to get model equation 'eqNumber'
    model = modU.modEquations(eqNumber)
    xdata = residenceTimes
    #Initial guess for variables A, f, tau in that order 
    #tempData = model(xdata, 1, 0.5, 1.5)
    ydata = probability
    #Use curve_fit with given times and probability on the given eq. Number
    fitParams, fitCovariances = curve_fit(model, xdata, ydata)
    return fitParams
    
    
