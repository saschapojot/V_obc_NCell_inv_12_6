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
    print("NLags="+str(NLags))
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

def checkU_distDataFilesForOneT(U_dist_data_dir):
    """

    :param U_dist_data_dir:
    :return:
    """
    U_dist_sortedDataFilesToRead=sort_data_files_by_loopEnd(U_dist_data_dir)
    if len(U_dist_sortedDataFilesToRead)==0:
        print("no data for U_dist.")
        exit(0)

    startingFileInd,startingVecPosition=parseSummaryU_Dist()

    if startingFileInd<0:
        #we guess that the equilibrium starts at this file
        startingFileInd=int(len(U_dist_sortedDataFilesToRead)*5/7)
    startingFileName=U_dist_sortedDataFilesToRead[startingFileInd]
    # print(startingFileName)
    #read the starting U_dist pkl file
    # in_dfStart=pd.read_csv(startingFileName)
    print(startingFileName)
    with open(startingFileName,"rb") as fptr:
        inArrStart=pickle.load(fptr)
    inArrStart=np.reshape(inArrStart,(-1,2*N+1))
    in_nRowStart,in_nColStart=inArrStart.shape
    if startingVecPosition<0:
        #we guess equilibrium starts at this position
        startingVecPosition=int(in_nRowStart/2)

    UVec=inArrStart[startingVecPosition:,0]
    xAxB_array=inArrStart[startingVecPosition:,1:]
    # print(UVec[:10])

    #read the rest of the pkl files
    for pkl_file in U_dist_sortedDataFilesToRead[(startingFileInd+1):]:
        print("reading: "+str(pkl_file))
        with open(pkl_file,"rb") as fptr:
            inArr=np.reshape(pickle.load(fptr),(-1,2*N+1))


        in_U=inArr[:,0]
        in_xAxB_array=inArr[:,1:]


        UVec=np.r_[UVec,in_U]
        xAxB_array=np.r_[xAxB_array,in_xAxB_array]


    nRow,nCol=xAxB_array.shape

    unitCellNum=N
    # print(UVec[-10:])
    xA_array=xAxB_array[:,:unitCellNum]
    xB_array=xAxB_array[:,unitCellNum:]

    d1_array=np.zeros((nRow,unitCellNum),dtype=float)

    d2_array=np.zeros((nRow,unitCellNum-1),dtype=float)

    for j in range(0,unitCellNum):
        d1_array[:,j]=xB_array[:,j]-xA_array[:,j]
    # print("unitCellNum-2="+str(unitCellNum-2))
    for j in range(0,unitCellNum-1):
        d2_array[:,j]=xA_array[:,j+1]-xB_array[:,j]


    dist_same=[]
    dist_lags=[]

    sameUTmp,lagUTmp=auto_corrForOneColumn(UVec)
    dist_same.append(sameUTmp)
    dist_lags.append(lagUTmp)

    for j in range(0,unitCellNum):
        sameTmp,lagTmp=auto_corrForOneColumn(d1_array[:,j])
        dist_same.append(sameTmp)
        dist_lags.append(lagTmp)

    for j in range(0,unitCellNum-1):
        sameTmp,lagTmp=auto_corrForOneColumn(d2_array[:,j])
        dist_same.append(sameTmp)
        dist_lags.append(lagTmp)

    same_exist=any(dist_same)

    #all values are the same, exit with err code
    if same_exist==True:
        with open(summary_U_distFile,"w+") as fptr:
            msg="error: same\n"
            fptr.writelines(msg)
            exit(sameErrCode)

    pThreshHold=0.01
    #if one lag==-1, then the auto-correlation is too large

    if np.min(dist_lags)>0:
        lagMax=np.max(dist_lags)
        pValsAll=[]
        lengthValAll=[]
        statsValsAll=[]

        pUTmp,statUTmp,lengthUTmp=ksTestOneColumn(UVec,lagMax)
        pValsAll.append(pUTmp)
        lengthValAll.append(lengthUTmp)
        statsValsAll.append(statUTmp)

        for j in range(0,unitCellNum):
            pTmp,statTmp,lengthTmp=ksTestOneColumn(d1_array[:,j],lagMax)
            pValsAll.append(pTmp)
            lengthValAll.append(lengthTmp)
            statsValsAll.append(statTmp)

        for j in range(0,unitCellNum-1):
            pTmp,statTmp,lengthTmp=ksTestOneColumn(d2_array[:,j],lagMax)
            pValsAll.append(pTmp)
            lengthValAll.append(lengthTmp)
            statsValsAll.append(statTmp)
        # print(pValsAll)
        print("stats="+str(statsValsAll))
        numDataPoints=np.min(lengthValAll)
        print("pValsAll="+str(pValsAll))
        if np.min(pValsAll)>=pThreshHold and numDataPoints>=200:
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

        continueMsg="continue\n"

        if np.min(pValsAll)<pThreshHold:
            #not the same distribution

            continueMsg+="p value: "+str(np.min(pValsAll))+"\n"
        if numDataPoints<200:
            #not enough data number

            continueMsg+="numDataPoints="+str(numDataPoints)+" too low\n"
            # continueMsg+="lag="+str(lagMax)+"\n"
        continueMsg+="lag="+str(lagMax)+"\n"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(continueMsg)
        exit(0)

    else:
        msg="high correlation"
        with open(summary_U_distFile,"w+") as fptr:
            fptr.writelines(msg)
        exit(0)









checkU_distDataFilesForOneT(U_dist_dataDir)