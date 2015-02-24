import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../../Utilities')
from scipy.optimize import curve_fit
import Utilities as util
import ModUtilities as Mod
import PlotUtilities as pltUtil

def hyperfit(t,a,b):
    return np.exp(hyperfitLog(t,a,b))

def fitStr(params,paramsStd,predicted,actual,returnParams=False):
    modelStr = 'e^(-b*(t/t+a))'
    labels=['a','b']
    return modelStrGen(params,paramsStd,predict,actual,labels,modelStr,
                       returnParams)
    

def hyperfitLogStr(params,paramsStd,predicted,actual,returnParams=False):
    modelStr = 'b * t/(t+a)'
    labels=['a','b']
    return modelStrGen(params,paramsStd,predicted,actual,labels,modelStr,
                       returnParams)

def hyperfitLog(t,a,b):
    return  (-b*t/(t+a))

CRTD = [1.0 , 0.48386793340231393, 0.23356222368544821, 0.14250776032358203, 0.093876399209858019, 0.066973944125670259, 0.05041858715078551, 0.038848650174019395, 0.03235819772363846, 0.028689681121249255, 0.0242686482927289, 0.021917035086069125, 0.018718841125011876, 0.016649421503151296, 0.015050324522622671, 0.012980904900762091, 0.012040259618098181, 0.011381807920233467, 0.010535227165835992, 0.010158969052770472, 0.0092183237701065623, 0.0079014203743771327, 0.0075251622613116131, 0.0071489041482460935, 0.0066785815069141385, 0.0061141943373157481, 0.0056438716959837931, 0.0053616781111845979, 0.0051735490546518381, 0.0050794845263854027, 0.0049854199981189673, 0.0047972909415862075, 0.0044210328285206879, 0.0040447747154551683, 0.0039507101871887329, 0.0038566456589222975, 0.0035744520741231023, 0.0033863230175903425, 0.0032922584893239071, 0.0031981939610574717, 0.0031041294327910363, 0.0030100649045246008, 0.0029160003762581654, 0.00282193584799173, 0.0027278713197252946, 0.0026338067914588592, 0.0025397422631924238, 0.0024456777349259884, 0.0022575486783932286, 0.0021634841501267932, 0.0020694196218603578, 0.001881290565327598, 0.0016931615087948382, 0.0015990969805284028, 0.0015050324522619674, 0.0014109679239955319, 0.0013169033957290965, 0.0011287743391963367, 0.00094064528266357694, 0.00084658075439714153, 0.00075251622613070612, 0.00065845169786427071, 0.0005643871695978353, 0.00047032264133139989, 0.00037625811306496448, 0.00028219358479852907, 0.00018812905653209366]
times  = [0.0, 0.10000000000000001, 0.20000000000000001, 0.30000000000000004, 0.40000000000000002, 0.5, 0.60000000000000009, 0.70000000000000007, 0.80000000000000004, 0.90000000000000002, 1.0, 1.1000000000000001, 1.2000000000000002, 1.3, 1.4000000000000001, 1.5, 1.6000000000000001, 1.7000000000000002, 1.8, 1.9000000000000001, 2.0, 2.1000000000000001, 2.2000000000000002, 2.3000000000000003, 2.4000000000000004, 2.5, 2.6000000000000001, 2.7000000000000002, 2.8000000000000003, 2.9000000000000004, 3.0, 3.1000000000000001, 3.2000000000000002, 3.4000000000000004, 3.5, 3.6000000000000001, 4.0, 4.1000000000000005, 4.2000000000000002, 4.2999999999999998, 4.4000000000000004, 4.5, 4.7000000000000002, 4.8000000000000007, 4.9000000000000004, 5.0, 5.3000000000000007, 5.4000000000000004, 5.5, 5.7000000000000002, 6.1000000000000005, 6.2000000000000002, 6.7000000000000002, 6.9000000000000004, 7.8000000000000007, 9.0, 9.4000000000000004, 9.6000000000000014, 9.9000000000000004, 11.0, 11.4, 11.700000000000001, 12.0, 12.100000000000001, 12.700000000000001, 16.400000000000002, 21.700000000000003]

times = np.array(times)
CRTD = np.array(CRTD)
logP = np.log(CRTD)
# fit to the normal CRTD
popt, pcov = curve_fit(hyperfit,times,CRTD)
paramsStd = np.sqrt(np.diag(pcov))
# get the predicted values
predict = hyperfit(times,*popt)
# get the string for the model label
modelStr,modelParams = fitStr(popt,paramsStd,predict,CRTD,True)
fig = pltUtil.pFigure()
numPlots = 3
plotCount = 1
xLabelStr = 'Time to unfold, t (seconds)'
# plot it all
timeConstant = popt[0]
timeLabel = modelParams[0]
ax = plt.subplot(numPlots,1,plotCount)
ax.axvline(timeConstant,color='k',label=timeLabel)
ax.plot(times,CRTD,'ro',label="Data")
ax.plot(times,predict,'b-',label=modelStr)
ax.set_ylim([min(CRTD),max(CRTD)])
ax.set_yscale('log')
plt.grid(b=True, which='major', color='b', linestyle='-')
plt.grid(b=True, which='minor', color='r', linestyle='--')
plt.legend(loc='best')
plt.title('CRTD fit by exponential fit')
plt.ylabel('Unfolding probability (>=t), CRTD')
plt.xlabel(xLabelStr)
# next, plot log(1/CRTD) = -log(CRTD) and the model...
plotCount += 1
plt.subplot(numPlots,1,plotCount)
# get the logatirhmic fit...
yVals = logP
popt, pcov = curve_fit(hyperfitLog,times,yVals)
paramsStd = np.sqrt(np.diag(pcov))
# get the predicted values
predictLog = hyperfitLog(times,*popt)
# get the string for the model label
modelStr = hyperfitLogStr(popt,paramsStd,predictLog,yVals)
plt.plot(times,yVals,'kx-',label='data')
plt.plot(times,predictLog,'b-',label=modelStr)
plt.title('Logarithmic fit results in slightly different parameters.')
plt.ylabel('CRTD')
plt.xlabel(xLabelStr)
plt.legend()
plotCount += 1
plt.subplot(numPlots,1,plotCount)
expErrRel=np.abs(predict-CRTD)/CRTD
plt.semilogy(times,expErrRel,'ro-',
             label='Error, Exponential Fit')

logErrRel = np.abs((predictLog-yVals)/yVals)
# divide by zero will give us an nan, impossibe to plot
# so, if we find NAN, (which isn't equal to itself), set to
# zero. kinda hacky, but I wrote this on a plane..
plt.title('Logarithmic fit has lower average error than exponential fit')
logErrRel[np.where(logErrRel != logErrRel)[0]] = 0
plt.semilogy(times,logErrRel,'kx-',
             label='Error, Linear fit')
# get the average (mean?) errors. first order for fit goodness.
averageExp = np.mean(expErrRel)
averageLog = np.mean(logErrRel)
plt.axhline(averageExp,color='r',linestyle='--',
            label='Average Error, Exp: {:4.2f}'.format(averageExp))
plt.axhline(averageLog,color='k',linestyle='--',
            label='Average Error, Log: {:4.2f}'.format(averageLog))
plt.ylabel('|Residuals|/CRTD)')
plt.xlabel(xLabelStr)

plt.legend(loc='best')
plotCount += 1

plt.tight_layout()
pltUtil.saveFigure(fig,'../../output/bestFit',overrideIO=False,close=False)

