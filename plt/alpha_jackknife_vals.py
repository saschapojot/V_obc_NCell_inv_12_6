import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats

#This script loads csv data and plot alpha, with condidence interval

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()
rowNum=0#int(sys.argv[1])
unitCellNum=int(sys.argv[1])


csvDataFolderRoot="../dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/csvOutAll/"

inCsvFile="../V_inv_12_6Params.csv"

TVals=[]
TFileNames=[]

for TFile in glob.glob(csvDataFolderRoot+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]


def alpha_estimator(LVec, UVec,T):
    """

    :param LVec:
    :param UVec:
    :param T:
    :return:
    """
    LVec=np.array(LVec)
    UVec=np.array(UVec)
    LU=LVec*UVec

    LUMean=np.mean(LU)

    LMean=np.mean(LVec)

    UMean=np.mean(UVec)

    alpha=1/T**2*1/LMean*(LUMean-LMean*UMean)/unitCellNum

    return alpha


def alpha_jackknife(LVec, UVec,T):
    """

    :param LVec:
    :param UVec:
    :param T:
    :return:
    """
    n=len(UVec)
    jackknife_samples = np.zeros(n)
    for i in range(0,n):
        sampleL=np.delete(LVec,i)
        sampleU=np.delete(UVec,i)
        jackknife_samples[i]=alpha_estimator(sampleL,sampleU,T)

    # Jackknife estimate of the statistic
    jackknife_estimate = np.mean(jackknife_samples)
    variance_estimate = (n - 1) / n * np.sum((jackknife_samples - jackknife_estimate) ** 2)

    return jackknife_estimate, variance_estimate

def alpha_confidence_interval(LVec, UVec,T,confidence_level=0.95):
    """

    :param LVec:
    :param UVec:
    :param T:
    :param confidence_level:
    :return:
    """
    jackknife_estimate, jackknife_variance=alpha_jackknife(LVec, UVec,T)
    n=len(UVec)
    alpha = 1 - confidence_level
    t_critical = stats.t.ppf(1 - alpha / 2, df=n-1)
    # Calculate the standard error
    standard_error = np.sqrt(jackknife_variance)

    # Calculate the confidence interval
    ci_lower = jackknife_estimate - t_critical * standard_error
    ci_upper = jackknife_estimate + t_critical * standard_error

    return jackknife_estimate,ci_lower, ci_upper


def generate_one_alpha_point(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return:
    """
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    U_distPath=oneTFile+"/U_dist/U_distData.csv"

    df=pd.read_csv(U_distPath)

    UVec=np.array(df.iloc[:,0])

    print("T="+str(TVal)+", data num="+str(len(UVec)))
    distArray=np.array(df.iloc[:,1:])
    nRow,nCol=distArray.shape
    xA_array=distArray[:,:unitCellNum]
    xB_array=distArray[:,unitCellNum:]
    LVec=xB_array[:,-1]-xA_array[:,0]

    jackknife_estimate,ci_lower, ci_upper=alpha_confidence_interval(LVec,UVec,TVal)

    return [jackknife_estimate,ci_lower, ci_upper]



alphaValsAll=[]
interval_lowerValsAll=[]
interval_upperValsAll=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    jackknife_estimate,ci_lower, ci_upper=generate_one_alpha_point(oneTFile)
    alphaValsAll.append(jackknife_estimate)
    interval_lowerValsAll.append(ci_lower)
    interval_upperValsAll.append(ci_upper)

sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals<100)
TToPlt=sortedTVals[TInds]
interval_lowerValsAll=np.array(interval_lowerValsAll)
interval_upperValsAll=np.array(interval_upperValsAll)
alphaValsAll=np.array(alphaValsAll)

alpha_err_var=alphaValsAll-interval_lowerValsAll

#plt alpha
fig,ax=plt.subplots()
ax.errorbar(TToPlt,alphaValsAll[TInds],yerr=alpha_err_var[TInds],fmt='o',color="blue", ecolor='r', capsize=5,label='mc')

ax.set_xlabel('$T$')
ax.set_ylabel("$\\alpha$")
ax.set_title("$\\alpha$ per unit cell, unit cell number="+str(unitCellNum))
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/alpha.png")
plt.close()