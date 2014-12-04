# use numpy for mathmatical functions
import numpy as np

def modEq1(t, A, f, tau):
    #Equation 2 from Walder 2014
    return A * np.sum(f*np.exp(-t/tau))
    
def modEq2(t, f, r, D, tau):
    #Equation 6 from Walder 2014
    return np.sum(f * np.exp((-r**2)/(4 * D * t)))

def modEquations(eqNumber):
    if (eqNumber == 1):
        #Equation 2 from Walder 2014
        model = modEq1
        
    if (eqNumber == 2):
        #Equation 6 from Walder 2014
        model = modEq2
        
    #Return the equation wanted as 'model'
    return model
        
