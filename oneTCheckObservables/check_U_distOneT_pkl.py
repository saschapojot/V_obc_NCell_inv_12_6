import pickle
import numpy as np
from datetime import datetime
from multiprocessing import Pool
import pandas as pd
import statsmodels.api as sm
import sys
import re
import warnings
from scipy.stats import ks_2samp
import glob
from pathlib import Path
import os
import json
import pickle

#This script checks if U, xA, xB values reach equilibrium and writes summary file of dist
#This file checks pkl files

argErrCode=2
sameErrCode=3
missingErrCode=4
if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit(argErrCode)

# print("entering")
jsonFromSummaryLast=json.loads(sys.argv[1])
jsonDataFromConf=json.loads(sys.argv[2])
# print(jsonFromSummaryLast)

TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dist_dataDir=jsonFromSummaryLast["U_dist_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])
N=int(jsonDataFromConf["unitCellNum"])

summary_U_distFile=TDirRoot+"/summary_U_dist.txt"
# print(summary_U_distFile)
lastFileNum=10
def sort_data_files_by_loopEnd(oneDir):
    dataFilesAll=[]
    loopEndAll=[]
    # print("entering sort")
    for oneDataFile in glob.glob(oneDir+"/*.pkl"):
        # print(oneDataFile)
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"loopEnd(\d+)",oneDataFile)
        if matchEnd:
            indTmp=int(matchEnd.group(1))
            loopEndAll.append(indTmp)
    endInds=np.argsort(loopEndAll)
    sortedDataFiles=[dataFilesAll[i] for i in endInds]
    return sortedDataFiles


def parseSummaryU_Dist():
    startingFileInd=-1
    startingVecPosition=-1

    summaryFileExists=os.path.isfile(summary_U_distFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition

    with open(summary_U_distFile,"r") as fptr:
        lines=fptr.readlines()
    for oneLine in lines:
        #match startingFileInd
        matchStartingFileInd=re.search(r"startingFileInd=(\d+)",oneLine)
        if matchStartingFileInd:
            startingFileInd=int(matchStartingFileInd.group(1))

        #match startingVecPosition
        matchStartingVecPosition=re.search(r"startingVecPosition=(\d+)",oneLine)
        if matchStartingVecPosition:
            startingVecPosition=int(matchStartingVecPosition.group(1))

    return startingFileInd, startingVecPosition


def auto_corrForOneColumn(colVec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    eps=5e-2
    NLags=int(len(colVec)*1/4)
    # print("NLags="+str(NLags))
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(colVec,nlags=NLags)
    except Warning as w:
        same=True
    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)

    lagVal=-1
    if minAutc<=eps:
        lagVal=np.where(acfOfVecAbs<=eps)[0][0]
    # np.savetxt("autc.txt",acfOfVecAbs[lagVal:],delimiter=',')
    return same,lagVal

def compute_acf_for_one_lag(x, lag):
    mean_x = np.mean(x)
    var_x = np.var(x)
    n = len(x)

    if lag == 0:
        return 1.0  # ACF for lag 0 is always 1
    else:
        cov = np.sum((x[:n-lag] - mean_x) * (x[lag:] - mean_x)) / n
        return cov / var_x

def ksTestOneColumn(colVec,lag):
    """

    :param colVec: a vector of data
    :param lag: auto-correlation length
    :return:
    """
    colVecSelected=colVec[::lag]

    lengthTmp=len(colVecSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    lenPart=int(lengthTmp/2)

    colVecToCompute=colVecSelected[-lengthTmp:]

    #ks test
    selectedVecPart0=colVecToCompute[:lenPart]
    selectedVecPart1=colVecToCompute[lenPart:]
    result=ks_2samp(selectedVecPart0,selectedVecPart1)
    return result.pvalue,result.statistic, lenPart*2

def checkUDataFilesForOneT(UData_dir,startingFileFraction, startingRowFraction):
    """

    :param UData_dir:
    :return:
    """
    U_sortedDataFilesToRead=sort_data_files_by_loopEnd(UData_dir)
    if len(U_sortedDataFilesToRead)==0:
        print("no data for U.")
        exit(0)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        # startingFileInd=int(len(U_sortedDataFilesToRead)*startingFileFraction)
        startingFileInd=len(U_sortedDataFilesToRead)-lastFileNum
    startingFileName=U_sortedDataFilesToRead[startingFileInd]
    # print(startingFileName)
    #read the starting U pkl file

    # print(startingFileName)

    with open(startingFileName,"rb") as fptr:
        inArrStart=pickle.load(fptr)



    in_nRowStart=len(inArrStart)
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart*startingRowFraction)
    arr=inArrStart[startingVecPosition:]


    #read the rest of the pkl files
    for pkl_file in U_sortedDataFilesToRead[(startingFileInd+1):]:
        # print("reading: "+str(pkl_file))
        with open(pkl_file,"rb") as fptr:
            inArr=pickle.load(fptr)
        arr=np.append(arr,inArr)


    sameUTmp,lagUTmp=auto_corrForOneColumn(arr)


    #if one lag==-1, then the auto-correlation is too large

    if sameUTmp==True or lagUTmp==-1:
        return [sameUTmp,lagUTmp,-1,-1,-1,-1,-1]

    pUTmp,statUTmp,lengthUTmp=ksTestOneColumn(arr,lagUTmp)
    numDataPoints=lengthUTmp

    return [sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPoints,startingFileInd,startingVecPosition]


def check_oneDistDataFilesForOneT(x1Dir, x2Dir,startingFileFraction, startingRowFraction):
    """

    :param x1Dir:
    :param x2Dir:
    :param startingFileFraction:
    :param startingRowFraction:
    :return:
    """
    x1sortedDataFilesToRead=sort_data_files_by_loopEnd(x1Dir)
    x2sortedDataFilesToRead=sort_data_files_by_loopEnd(x2Dir)
    if len(x1sortedDataFilesToRead)!= len(x2sortedDataFilesToRead):
        print("data missing.")
        exit(missingErrCode)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        # startingFileInd=int(len(x1sortedDataFilesToRead)*startingFileFraction)
        startingFileInd=len(x1sortedDataFilesToRead)-lastFileNum

    x1StartingFileName=x1sortedDataFilesToRead[startingFileInd]
    # print(x1StartingFileName)
    with open(x1StartingFileName,"rb") as fptr:
        x1_inArrStart=pickle.load(fptr)

    in_nRowStart=len(x1_inArrStart)
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart*startingRowFraction)

    x1Arr=x1_inArrStart[startingVecPosition:]

    #read the rest of the x1 pkl files
    for pkl_file in x1sortedDataFilesToRead[(startingFileInd+1):]:
        # print("reading: "+str(pkl_file))
        with open(pkl_file,"rb") as fptr:
            x1_inArr=pickle.load(fptr)
            x1Arr=np.append(x1Arr,x1_inArr)


    x2StartingFileName=x2sortedDataFilesToRead[startingFileInd]
    # print(x2StartingFileName)
    with open(x2StartingFileName,"rb") as fptr:
        x2_inArrStart=pickle.load(fptr)

    x2Arr=x2_inArrStart[startingVecPosition:]
    #read the rest of the x2 pkl files
    for pkl_file in x2sortedDataFilesToRead[(startingFileInd+1):]:
        # print("reading: "+str(pkl_file))
        with open(pkl_file,"rb") as fptr:
            x2_inArr=pickle.load(fptr)

            x2Arr=np.append(x2Arr,x2_inArr)
    x1Arr=np.array(x1Arr)
    x2Arr=np.array(x2Arr)
    distArr=x2Arr-x1Arr
    # print(distArr[-2])
    sameTmp,lagTmp=auto_corrForOneColumn(distArr)
    if sameTmp==True or lagTmp==-1:
        return [sameTmp,lagTmp,-1,-1,-1,-1,-1]

    pTmp,statTmp,lengthTmp=ksTestOneColumn(distArr,lagTmp)
    numDataPoints=lengthTmp

    return [sameTmp,lagTmp,pTmp,statTmp,numDataPoints,startingFileInd,startingVecPosition]


def wrapper_checkUDataFilesForOneT(input):
    UDataDir,startingFileFraction,startingRowFraction=input
    return checkUDataFilesForOneT(UDataDir,startingFileFraction,startingRowFraction)

def wrapper_check_oneDistDataFilesForOneT(input):
    x1Dir, x2Dir,startingFileFraction, startingRowFraction=input
    return check_oneDistDataFilesForOneT(x1Dir, x2Dir,startingFileFraction, startingRowFraction)


startingFileFraction=5/7
startingRowFraction=1/2
UDataDir=U_dist_dataDir+"/U/"


sameVec=[]
lagVec=[]
pVec=[]
statVec=[]
numDataVec=[]
# print("before ")
print("checking U")
sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPointsU,startingFileInd,startingVecPosition=checkUDataFilesForOneT(UDataDir,startingFileFraction,startingRowFraction)
print("lag="+str(lagUTmp))
sameVec.append(sameUTmp)
lagVec.append(lagUTmp)
pVec.append(pUTmp)
statVec.append(statUTmp)
numDataVec.append(numDataPointsU)

def generate_xAj_path(j):
    xAjPath=U_dist_dataDir+"/xA"+str(j)+"/"
    return xAjPath

def generate_xBj_path(j):
    xBjPath=U_dist_dataDir+"/xB"+str(j)+"/"
    return xBjPath

xAPathAll=[generate_xAj_path(j) for j in range(0,N)]
xBPathAll=[generate_xBj_path(j) for j in range(0,N)]

#checking L
# xBPathLast=xBPathAll[-1]
# xAPath0=xAPathAll[0]
# sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xAPath0,xBPathLast,startingFileFraction,startingRowFraction)
# print("lagTmp="+str(lagTmp))
# sameVec.append(sameTmp)
# lagVec.append(lagTmp)
# pVec.append(pTmp)
# statVec.append(statTmp)
# numDataVec.append(numDataPoints)
##################################################serial
#check xBj - xAj, serial
# tDist1Start=datetime.now()
# for j in range(0,N):
#     xBjPathTmp=xBPathAll[j]
#     xAjPathTmp=xAPathAll[j]
#     print("checking xB"+str(j)+"-"+"xA"+str(j))
#     sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xAjPathTmp,xBjPathTmp,startingFileFraction,startingRowFraction)
#     print("lagTmp="+str(lagTmp))
#     sameVec.append(sameTmp)
#     lagVec.append(lagTmp)
#     pVec.append(pTmp)
#     statVec.append(statTmp)
#     numDataVec.append(numDataPoints)
#
# tDist1End=datetime.now()
# print("check dist1 time: ",tDist1End-tDist1Start)

#check xAj+1 - xBj, serial
# for j in range(0,N-1):
#     xAjPlus1PathTmp=xAPathAll[j+1]
#     xBjPathTmp=xBPathAll[j]
#     print("checking xA"+str(j+1)+"-"+"xB"+str(j))
#     sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xBjPathTmp,xAjPlus1PathTmp,startingFileFraction,startingRowFraction)
#     print("lagTmp="+str(lagTmp))
#     sameVec.append(sameTmp)
#     lagVec.append(lagTmp)
#     pVec.append(pTmp)
#     statVec.append(statTmp)
#     numDataVec.append(numDataPoints)

##################################################end serial

################################################## parallel

#dist 1
# tDist1Start=datetime.now()
# procNum=12
# inputVec1=[[xAPathAll[j],xBPathAll[j],startingFileFraction,startingRowFraction] for j in range(0,N)]
#
# pool1=Pool(procNum)
# ret1=pool1.map(wrapper_check_oneDistDataFilesForOneT,inputVec1)
# tDist1End=datetime.now()
# print("check dist1 time: ",tDist1End-tDist1Start)
#
# for item in ret1:
#     sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=item
#     sameVec.append(sameTmp)
#     lagVec.append(lagTmp)
#     pVec.append(pTmp)
#     statVec.append(statTmp)
#     numDataVec.append(numDataPoints)

# dist 2
# tDist2Start=datetime.now()
# inputVec2=[[xBPathAll[j],xAPathAll[j+1],startingFileFraction,startingRowFraction] for j in range(0,N-1)]
# pool2=Pool(procNum)
# ret2=pool2.map(wrapper_check_oneDistDataFilesForOneT,inputVec2)
# tDist2End=datetime.now()
# print("check dist2 time: ",tDist2End-tDist2Start)
# for item in ret2:
#     sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=item
#     sameVec.append(sameTmp)
#     lagVec.append(lagTmp)
#     pVec.append(pTmp)
#     statVec.append(statTmp)
#     numDataVec.append(numDataPoints)

################################################## end parallel


############################################
#manual checking
values1=list(range(0,N))
j1=np.random.choice(values1, size=1, replace=True)[0]
# print(j1)
# j1=10
xBjPathTmp=xBPathAll[j1]
xAjPathTmp=xAPathAll[j1]
print("checking xB"+str(j1)+"-"+"xA"+str(j1))
sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xAjPathTmp,xBjPathTmp,startingFileFraction,startingRowFraction)
print("lagTmp="+str(lagTmp))
sameVec.append(sameTmp)
lagVec.append(lagTmp)
pVec.append(pTmp)
statVec.append(statTmp)
numDataVec.append(numDataPoints)

values2=list(range(0,N-1))
j2=np.random.choice(values2, size=1, replace=True)[0]

xAjPlus1PathTmp=xAPathAll[j2+1]
xBjPathTmp=xBPathAll[j2]
print("checking xA"+str(j2+1)+"-"+"xB"+str(j2))
sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xBjPathTmp,xAjPlus1PathTmp,startingFileFraction,startingRowFraction)
sameVec.append(sameTmp)
lagVec.append(lagTmp)
pVec.append(pTmp)
statVec.append(statTmp)
numDataVec.append(numDataPoints)
print("lagTmp="+str(lagTmp))

############################################
summary_U_distFile=TDirRoot+"/summary_U_dist.txt"

print(lagVec)
if np.min(lagVec)<0:
    msg="high correlation"
    with open(summary_U_distFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)

same_exist=any(sameVec)
if same_exist==True:
    with open(summary_U_distFile,"w+") as fptr:
        msg="error: same\n"
        fptr.writelines(msg)
        exit(sameErrCode)

#equilibrium
# pThreshHold=0.01
lagMax=np.max(lagVec)
statThreshhold=0.1
print("statVec="+str(statVec))
print("pVec="+str(pVec))

numDataPoints=np.min(numDataVec)
if np.max(statVec)<=statThreshhold and numDataPoints>=200:
    if numDataPoints>=effective_data_num_required:
        newDataPointNum=0
    else:
        newDataPointNum=effective_data_num_required-numDataPoints
    msg="equilibrium\n" \
        +"lag="+str(lagMax)+"\n" \
        +"numDataPoints="+str(numDataPoints)+"\n" \
        +"startingFileInd="+str(startingFileInd)+"\n" \
        +"startingVecPosition="+str(startingVecPosition)+"\n" \
        +"newDataPointNum="+str(newDataPointNum)+"\n"

    with open(summary_U_distFile,"w+") as fptr:
        fptr.writelines(msg)
    exit(0)

#continue
continueMsg="continue\n"
if np.max(statVec)>statThreshhold:
    #not the same distribution
    continueMsg+="stat value: "+str(np.max(statVec))+"\n"

if numDataPoints<200:
    #not enough data number

    continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
    continueMsg+="lag="+str(lagMax)+"\n"
with open(summary_U_distFile,"w+") as fptr:
    fptr.writelines(continueMsg)
exit(0)

