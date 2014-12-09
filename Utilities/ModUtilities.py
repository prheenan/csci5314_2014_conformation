# use numpy for mathmatical functions
import numpy as np

def modelString(equation,values,stdevs,labels,units,RSQ):
    mStr = "Equation: {:s}, RSQ: {:.3f}".format(equation,RSQ)
    for val,std,label,unit in zip(values,stdevs,labels,units):
        mStr += "\n{:}: {:3.2g} +/- {:2.1g} [{:}]".format(label,val,std,unit)
    return mStr

def RSQ(yPred,yActual):
    mean = np.mean(yActual)
    SStot = np.sum((yActual-mean)**2)
    SSres = np.sum((yPred-yActual)**2)
    return 1-SSres/SStot

def modEq1(t, A, f, tau):
    #Equation 2 from Walder 2014
    toRet = A * f*np.exp(-t/tau)
    return toRet

def modEqx(t, f, r, D, tau):
    #Equation 6 from Walder 2014
    return np.sum(f * np.exp((-r**2)/(4 * D * t)))

def modEq2(t,A,f1,f2,tau1,tau2):
    toRet = 0
    tau  = [tau1,tau2]
    f = [f1,f2]
    for tauI,fI in zip(tau,f):
        toRet += modEq1(t,A,tauI,fI)
    return toRet

def modEq3(t,A,f1,f2,f3,tau1,tau2,tau3):
    toRet = 0
    tau  = [tau1,tau2,tau3]
    f = [f1,f2,f3]
    for tauI,fI in zip(tau,f):
        toRet += modEq1(t,A,tauI,fI)
    return toRet


def modEquations(eqNumber):
    if (eqNumber == 1):
        #Equation 2 from Walder 2014
        model = modEq1
        
    if (eqNumber == 2):
        model = modEq2
    if (eqNumber == 3):
        model = modEq3
    if (eqNumber == 4):
        #Equation 6 from Walder 2014
        model = modEq4
        
    #Return the equation wanted as 'model'
    return model
        
