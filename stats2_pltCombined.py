import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
import glob

import pandas as pd
#this script reads from statsAllTxxx.txt and plot U, L, d1, d2





def removeCommentsAndEmptyLines(file):
    """

    :param file: conf file
    :return: contents in file, with empty lines and comments removed
    """
    with open(file,"r") as fptr:
        lines= fptr.readlines()

    linesToReturn=[]
    for oneLine in lines:
        oneLine = re.sub(r'#.*$', '', oneLine).strip()
        if not oneLine:
            continue
        else:
            linesToReturn.append(oneLine)
    return linesToReturn


def format_using_decimal(value):
    # Convert the float to a Decimal using string conversion to avoid precision issues
    decimal_value = Decimal(str(value))
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)



# if (len(sys.argv)!=2):
#     print("wrong number of arguments")
#     exit()
# T=float(sys.argv[1])
#
# TStr=format_using_decimal(T)

TRegex=r"T\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)"
def matchOneLine(line):
    matchTRegex=re.search(TRegex,line)
    if matchTRegex:
        TTmp=float(matchTRegex.group(1))
        matchU=re.search(r"mean_U\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        matchN=re.search(r"N\s*=\s*(\d+)",line)
        matchL=re.search(r"mean_L\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        match_hf_L=re.search(r"hf_L\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        match_hf_U=re.search(r"hf_U\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        if matchU and matchN and matchL and match_hf_L and match_hf_U:
            UMeanTmp=float(matchU.group(1))
            NTmp= int(matchN.group(1))
            LTmp=float(matchL.group(1))
            hf_LTmp=float(match_hf_L.group(1))
            hf_UTmp=float(match_hf_U.group(1))
            return [TTmp,NTmp,UMeanTmp,LTmp,hf_LTmp,hf_UTmp]
        else:
            print("mismatch in "+str(line))
statsFolder="./statsAll/"

statsFileVec=[]
for file in glob.glob(statsFolder+"/statsAllT*.txt"):
    statsFileVec.append(file)


def extractT(fileName):
    matchTRegex=re.search(r"statsAllT([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",fileName)
    if matchTRegex:
        return float(matchTRegex.group(1))



statsFile_TVec=[extractT(fileName) for fileName in statsFileVec]


sorted_by_TInds=np.argsort(statsFile_TVec)

sorted_TVec=[statsFile_TVec[ind] for ind in sorted_by_TInds]

sorted_statsFileVec=[statsFileVec[ind] for ind in sorted_by_TInds]

sorted_TStrVec=[format_using_decimal(T) for T in sorted_TVec]




def createVecFromOneFile(fileName):
    """

    :param fileName:
    :return: sorted_NVec, sorted_UMeanVec,sorted_LMeanVec,sorted_hf_LVec, sorted_hf_UVec
    """
    contentsInOneFile=removeCommentsAndEmptyLines(fileName)
    NVec=[]
    UMeanVec=[]
    LMeanVec=[]
    hf_LVec=[]
    hf_UVec=[]
    for line in contentsInOneFile:
        _,NTmp,UMeanTmp,LMeanTmp,hf_LTmp,hf_UTmp=matchOneLine(line)
        NVec.append(NTmp)
        UMeanVec.append(UMeanTmp)
        LMeanVec.append(LMeanTmp)
        hf_LVec.append(hf_LTmp)
        hf_UVec.append(hf_UTmp)
    sorted_NInds=np.argsort(NVec)
    sorted_NVec=[NVec[ind] for ind in sorted_NInds]
    sorted_UMeanVec=[UMeanVec[ind] for ind in sorted_NInds]
    sorted_LMeanVec=[LMeanVec[ind] for ind in sorted_NInds]
    sorted_hf_LVec=[hf_LVec[ind] for ind in sorted_NInds]
    sorted_hf_UVec=[hf_UVec[ind] for ind in sorted_NInds]

    sorted_NVec=np.array(sorted_NVec)
    sorted_UMeanVec=np.array(sorted_UMeanVec)
    sorted_LMeanVec=np.array(sorted_LMeanVec)
    sorted_hf_LVec=np.array(sorted_hf_LVec)
    sorted_hf_UVec=np.array(sorted_hf_UVec)

    sorted_UPerUnitCell=sorted_UMeanVec/sorted_NVec
    sorted_LPerUnitCell=sorted_LMeanVec/sorted_NVec

    return [sorted_NVec,sorted_UMeanVec,sorted_LMeanVec,sorted_hf_LVec,sorted_hf_UVec,sorted_UPerUnitCell,sorted_LPerUnitCell]

def createVecIndexed_byN(fileName):
    """
    same T, different N
    :param fileName:
    :return:
    """
    contentsInOneFile=removeCommentsAndEmptyLines(fileName)
    retArray=[]#[N,T,L,hf_L]


    for line in contentsInOneFile:
        TTmp,NTmp,_,LMeanTmp,hf_LTmp,_=matchOneLine(line)
        item=[int(NTmp),TTmp,LMeanTmp,hf_LTmp]
        retArray.append(item)
    return np.array(retArray)



#first index represents T


## E(U)/N
fig,ax=plt.subplots()
first_label = True  # Flag to ensure "mc" label is only added once
for j in range(0,len(sorted_statsFileVec)):
    TStrTmp=sorted_TStrVec[j]
    fileNameTmp=sorted_statsFileVec[j]
    sorted_NVecTmp,sorted_UMeanVecTmp,_,_,sorted_hf_UVecTmp,sorted_UPerUnitCellTmp,_=createVecFromOneFile(fileNameTmp)
    label_mc = 'mc' if first_label else None
    # ax.errorbar(sorted_NVecTmp,sorted_UPerUnitCellTmp,yerr=sorted_hf_UVecTmp/sorted_NVecTmp,fmt='.', color="black", ecolor='r', capsize=5, label=label_mc)
    ax.plot(sorted_NVecTmp,sorted_UPerUnitCellTmp,label="T="+TStrTmp,linestyle="--",linewidth=0.8,marker=".",markersize=4)
    first_label = False  # Turn off the label after the first iteration

ax.set_xlabel("$N$")
ax.set_ylabel("$E(U)/N$")
ax.set_title("E(U) per unit cell")
plt.legend(loc="best")
plt.savefig(statsFolder+"UPerAll"+".png")
plt.close()


## E(L)/N

fig,ax=plt.subplots()
first_label = True  # Flag to ensure "mc" label is only added once
L_scale=1
for j in range(0,len(sorted_statsFileVec)):
    TStrTmp=sorted_TStrVec[j]
    fileNameTmp=sorted_statsFileVec[j]
    sorted_NVecTmp,_,sorted_LMeanVecTmp,sorted_hf_LVecTmp,_,_,sorted_LPerUnitCellTmp=createVecFromOneFile(fileNameTmp)
    label_mc = 'mc' if first_label else None

    # ax.errorbar(sorted_NVecTmp,sorted_LPerUnitCellTmp,yerr=sorted_hf_LVecTmp/sorted_NVecTmp,fmt='.', color="black", ecolor='r', capsize=5, label=label_mc)
    ax.plot(sorted_NVecTmp,sorted_LPerUnitCellTmp*L_scale,label="T="+TStrTmp,linestyle="--",linewidth=1,marker=".",markersize=4)
    first_label = False  # Turn off the label after the first iteration

ax.set_xlabel("$N$")
ax.set_ylabel("$E(L)/N$")
ax.set_ylim(1.7,1.8)
ax.set_title("E(L) per unit cell")
plt.legend(loc="best")

plt.savefig(statsFolder+"LPerAll"+".png")
plt.close()


#E(L)
fig,ax=plt.subplots()
first_label = True  # Flag to ensure "mc" label is only added once
for j in range(0,len(sorted_statsFileVec)):
    TStrTmp=sorted_TStrVec[j]
    fileNameTmp=sorted_statsFileVec[j]
    sorted_NVec,_,sorted_LMeanVec,sorted_hf_LVec,_,_,sorted_LPerUnitCell=createVecFromOneFile(fileNameTmp)
    label_mc = 'mc' if first_label else None
    ax.plot(sorted_NVec,sorted_LMeanVec,label="T="+TStrTmp,linestyle="--",linewidth=1,marker=".",markersize=4)
    first_label = False  # Turn off the label after the first iteration

ax.set_xlabel("$N$")
ax.set_ylabel("$E(L)/N$")
# ax.set_ylim(40,42)
ax.set_title("E(L)")
plt.legend(loc="best")
# ax.set_yscale("log")
plt.savefig(statsFolder+"LAll"+".png")
plt.close()

#all L's table

table_L_all=createVecIndexed_byN(sorted_statsFileVec[0])
for j in range(1,len(sorted_statsFileVec)):
    fileNameTmp=sorted_statsFileVec[j]
    retArray=createVecIndexed_byN(fileNameTmp)
    table_L_all=np.r_[table_L_all,retArray]

# Get unique values of N
# Get unique values of N
unique_N_values = np.unique(table_L_all[:, 0])
fig,ax=plt.subplots()
first_label = True  # Flag to ensure "mc" label is only added once

for N in unique_N_values:
    if N<20:
        continue
    rows_for_N = table_L_all[table_L_all[:, 0] == N]
    TVec=rows_for_N[:,1]
    LVec=rows_for_N[:,2]
    LVec=np.array(LVec)
    hf_LVec=rows_for_N[:,3]
    label_mc = 'mc' if first_label else None

    ax.errorbar(TVec,LVec/N,yerr=hf_LVec/N,fmt='.', color="black", ecolor='r', capsize=5, label=label_mc)
    ax.plot(TVec,LVec/N,label="N="+str(int(N)),linestyle="--",linewidth=1,marker=".",markersize=4)
    first_label = False  # Turn off the label after the first iteration


ax.set_xlabel("$T$")
ax.set_ylabel("$E(L)$")
# ax.set_ylim(1.7,1.8)
ax.set_title("E(L) per unit cell")
plt.legend(loc="best")

plt.savefig(statsFolder+"LForN"+".png")
plt.close()