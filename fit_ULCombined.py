import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats
import os
from decimal import Decimal
import statsmodels.api as sm


#this script combines all fit data for E(U)/N
def format_using_decimal(value):
    # Convert the float to a Decimal using string conversion to avoid precision issues
    decimal_value = Decimal(str(value))
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)


#search U json files

dataRoot="./statsAll/"
fitU_jsonFileVec=[]
for file in glob.glob(dataRoot+"/fitU_T*.json"):
    fitU_jsonFileVec.append(file)

def extractT(fileName):
    matchTRegex=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",fileName)
    if matchTRegex:
        return float(matchTRegex.group(1))

json_TVec=[extractT(file ) for file in fitU_jsonFileVec]

sorted_by_TInds=np.argsort(json_TVec)
sorted_TVec=[json_TVec[ind] for ind in sorted_by_TInds]

sorted_TStrVec=[format_using_decimal(T) for T in sorted_TVec]

sorted_fitU_jsonFileVec=[fitU_jsonFileVec[ind] for ind in sorted_by_TInds]


def createDataFromOneFile(fileName):
    with open(fileName,"r") as fptr:
        data=json.load(fptr)
    N_inv0=data["N_inv0"]
    U0Pred=float(data["U0Pred"])
    N_invUnique=data["N_invUnique"]
    U_pred=data["U_pred"]
    L0Pred=float(data["L0Pred"])
    L_pred=data["L_pred"]
    return [N_inv0,U0Pred,N_invUnique,U_pred,L0Pred,L_pred]

# U
for j in range(0,len(sorted_fitU_jsonFileVec)):
    TStr=sorted_TStrVec[j]
    fileName=sorted_fitU_jsonFileVec[j]
    N_inv0,U0Pred,N_invUnique,U_pred,_,_=createDataFromOneFile(fileName)
    intercept=float(U0Pred)
    # plt.text(0, intercept, f'{intercept:.2f}', color='blue', verticalalignment='center', horizontalalignment='right')
    plt.plot(N_invUnique,U_pred,label="T="+TStr+", intercept="+str(np.round(intercept,2)),marker=".",markersize=3,linewidth=0.8,linestyle="--")

plt.xlabel(r"$\frac{1}{N}$")
plt.legend(loc="best")
plt.ylabel(r"E(U)/N")
plt.title("E(U) per unit cell")
plt.savefig(dataRoot+"/fit_UAll"+".png")
plt.close()

#L
for j in range(0,len(sorted_fitU_jsonFileVec)):
    TStr=sorted_TStrVec[j]
    fileName=sorted_fitU_jsonFileVec[j]
    N_inv0,_,N_invUnique,_,L0Pred,L_pred=createDataFromOneFile(fileName)
    intercept=L0Pred
    # plt.text(0, intercept, f'{intercept:.2f}', color='blue', verticalalignment='center', horizontalalignment='right')
    plt.plot(N_invUnique,L_pred,label="T="+TStr+", intercept="+str(np.round(intercept,2)),marker=".",markersize=3,linewidth=0.8,linestyle="-.")


plt.xlabel(r"$\frac{1}{N}$")
plt.legend(loc="best")
plt.ylabel(r"E(L)/N")
plt.title("E(L) per unit cell")
plt.savefig(dataRoot+"/fit_LAll"+".png")
plt.close()