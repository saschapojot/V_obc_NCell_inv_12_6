import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats

#This script loads csv data and plot C, with condidence interval

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()
rowNum=int(sys.argv[1])
unitCellNum=5


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

def C_estimator(UVec,U2Vec,T):
    """

    :param UVec: U values
    :param U2Vec: U^{2} values
    :param T: temperature
    :return: C
    """
    EU=np.mean(UVec)
    EU2=np.mean(U2Vec)

    CTmp=(EU2-EU**2)/T**2/unitCellNum
    return CTmp


def C_jackknife(UVec,U2Vec,T):
    """

    :param UVec: U values
    :param U2Vec: U^{2} values
    :param T: temperature
    :return: jackknife_estimate, variance_estimate
    """

    n=len(UVec)
    jackknife_samples = np.zeros(n)
    for i in range(0,n):
        sampleU=np.delete(UVec,i)
        sampleU2=np.delete(U2Vec,i)
        jackknife_samples[i]=C_estimator(sampleU,sampleU2,T)

    # Jackknife estimate of the statistic
    jackknife_estimate = np.mean(jackknife_samples)
    variance_estimate = (n - 1) / n * np.sum((jackknife_samples - jackknife_estimate) ** 2)

    return jackknife_estimate, variance_estimate



def C_confidence_interval(UVec,U2Vec,T,confidence_level=0.95):
    """

    :param UVec:
    :param U2Vec:
    :param T:
    :param confidence_level:
    :return:
    """
    jackknife_estimate, jackknife_variance=C_jackknife(UVec,U2Vec,T)
    n=len(UVec)
    alpha = 1 - confidence_level
    t_critical = stats.t.ppf(1 - alpha / 2, df=n-1)
    # Calculate the standard error
    standard_error = np.sqrt(jackknife_variance)

    # Calculate the confidence interval
    ci_lower = jackknife_estimate - t_critical * standard_error
    ci_upper = jackknife_estimate + t_critical * standard_error

    return jackknife_estimate,ci_lower, ci_upper





def generate_one_point(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return:
    """
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    U_distPath=oneTFile+"/U_dist/U_distData.csv"

    df=pd.read_csv(U_distPath)

    UVec=np.array(df.iloc[:,0])
    U2Vec=UVec**2
    print("T="+str(TVal)+", data num="+str(len(UVec)))

    jackknife_estimate,ci_lower, ci_upper=C_confidence_interval(UVec,U2Vec,TVal)
    return [jackknife_estimate,ci_lower, ci_upper]


CValsAll=[]
interval_lowerValsAll=[]
interval_upperValsAll=[]
for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    jackknife_estimate,ci_lower, ci_upper=generate_one_point(oneTFile)
    CValsAll.append(jackknife_estimate)
    interval_lowerValsAll.append(ci_lower)
    interval_upperValsAll.append(ci_upper)
sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals<100)
TToPlt=sortedTVals[TInds]
interval_lowerValsAll=np.array(interval_lowerValsAll)
interval_upperValsAll=np.array(interval_upperValsAll)
CValsAll=np.array(CValsAll)

C_err_bar=CValsAll-interval_lowerValsAll
#plt C
fig,ax=plt.subplots()
ax.errorbar(TToPlt,CValsAll[TInds],yerr=C_err_bar[TInds],fmt='o',color="black", ecolor='r', capsize=5,label='mc')

ax.set_xlabel('$T$')
ax.set_ylabel("C")
ax.set_title("C per unit cell, unit cell number="+str(unitCellNum))
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/C.png")
plt.close()