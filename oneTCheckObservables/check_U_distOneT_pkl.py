import pickle
import numpy as np
from datetime import datetime

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
# jsonFromSummaryLast={"startingFileInd": "74", "startingVecPosition": "500000", "newMcStepNum": "78940680", "newDataPointNum": "10920", "newFlushNum": "79", "TDirRoot": "./dataAllUnitCell10/row0/T0.5/", "U_dist_dataDir": "./dataAllUnitCell10/row0/T0.5//U_dist_dataFiles/"}

# jsonDataFromConf={"T": "0.5", "erase_data_if_exist": "False", "search_and_read_summary_file": "True", "observable_name": "U_dist", "potential_function_name": "V_inv_12_6", "effective_data_num_required": "15000", "loop_to_write": "1000000", "default_flush_num": "10", "coefs": "25,80,15,67", "confFileName": "./dataAllUnitCell10/row0/T0.5/run_T0.5.mc.conf", "unitCellNum": "10"}


TDirRoot=jsonFromSummaryLast["TDirRoot"]
U_dist_dataDir=jsonFromSummaryLast["U_dist_dataDir"]
effective_data_num_required=int(jsonDataFromConf["effective_data_num_required"])
N=int(jsonDataFromConf["unitCellNum"])

summary_U_distFile=TDirRoot+"/summary_U_dist.txt"
# print(summary_U_distFile)

def sort_data_files_by_loopEnd(oneDir):
    dataFilesAll=[]
    loopEndAll=[]
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
    eps=1e-2
    NLags=int(len(colVec)*1/10)
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
        startingFileInd=int(len(U_sortedDataFilesToRead)*startingFileFraction)
    startingFileName=U_sortedDataFilesToRead[startingFileInd]
    # print(startingFileName)
    #read the starting U pkl file

    print(startingFileName)

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
        return [sameUTmp,lagUTmp,-1,-1,-1]

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
        startingFileInd=int(len(x1sortedDataFilesToRead)*startingFileFraction)

    x1StartingFileName=x1sortedDataFilesToRead[startingFileInd]

    with open(x1StartingFileName,"rb") as fptr:
        x1_inArrStart=pickle.load(fptr)

    in_nRowStart=len(x1_inArrStart)
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart*startingRowFraction)

    x1Arr=x1_inArrStart[startingVecPosition:]

    #read the rest of the x1 pkl files
    for pkl_file in x1sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            x1_inArr=pickle.load(fptr)
        x1Arr=np.append(x1Arr,x1_inArr)


    x2StartingFileName=x2sortedDataFilesToRead[startingFileInd]

    with open(x2StartingFileName,"rb") as fptr:
        x2_inArrStart=pickle.load(fptr)

    x2Arr=x2_inArrStart[startingVecPosition:]
    #read the rest of the x2 pkl files
    for pkl_file in x2sortedDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            x2_inArr=pickle.load(fptr)

        x2Arr=np.append(x2Arr,x2_inArr)

    distArr=x2Arr-x1Arr
    # print(distArr[-2])
    sameTmp,lagTmp=auto_corrForOneColumn(distArr)
    if sameTmp==True or lagTmp==-1:
        return [sameTmp,lagTmp,-1,-1,-1]

    pTmp,statTmp,lengthTmp=ksTestOneColumn(distArr,lagTmp)
    numDataPoints=lengthTmp

    return [sameTmp,lagTmp,pTmp,statTmp,numDataPoints,startingFileInd,startingVecPosition]

startingFileFraction=5/7
startingRowFraction=1/2
UDataDir=U_dist_dataDir+"/U/"


sameVec=[]
lagVec=[]
pVec=[]
statVec=[]
numDataVec=[]

sameUTmp,lagUTmp,pUTmp,statUTmp,numDataPointsU,startingFileInd,startingVecPosition=checkUDataFilesForOneT(UDataDir,startingFileFraction,startingRowFraction)
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

#check xBj - xAj
for j in range(0,N):
    xBjPathTmp=xBPathAll[j]
    xAjPathTmp=xAPathAll[j]
    sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xAjPathTmp,xBjPathTmp,startingFileFraction,startingRowFraction)
    sameVec.append(sameTmp)
    lagVec.append(lagTmp)
    pVec.append(pTmp)
    statVec.append(statTmp)
    numDataVec.append(numDataPoints)



#check xAj+1 - xBj
for j in range(0,N-1):
    xAjPlus1PathTmp=xAPathAll[j+1]
    xBjPathTmp=xBPathAll[j]
    sameTmp,lagTmp,pTmp,statTmp,numDataPoints,_,_=check_oneDistDataFilesForOneT(xBjPathTmp,xAjPlus1PathTmp,startingFileFraction,startingRowFraction)
    sameVec.append(sameTmp)
    lagVec.append(lagTmp)
    pVec.append(pTmp)
    statVec.append(statTmp)
    numDataVec.append(numDataPoints)



summary_U_distFile=TDirRoot+"/summary_U_dist.txt"


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
pThreshHold=0.01
lagMax=np.max(lagVec)
print("statVec="+str(statVec))
print("pVec="+str(pVec))

numDataPoints=np.min(numDataVec)
if np.min(pVec)>=pThreshHold and numDataPoints>=200:
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
if np.min(pVec)<pThreshHold:
    #not the same distribution
    continueMsg+="p value: "+str(np.min(pVec))+"\n"
if numDataPoints<200:
    #not enough data number

    continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
    continueMsg+="lag="+str(lagMax)+"\n"
with open(summary_U_distFile,"w+") as fptr:
    fptr.writelines(continueMsg)
exit(0)

