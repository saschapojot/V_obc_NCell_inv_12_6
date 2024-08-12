import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import pandas as pd
import pickle
#this script extracts effective data from pkl files

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()


rowNum=int(sys.argv[1])

unitCellNum=int(sys.argv[2])
rowDirRoot="../dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/"
obs_U_dist="U_dist"

#search directory
TVals=[]
TFileNames=[]
TStrings=[]
for TFile in glob.glob(rowDirRoot+"/T*"):
    # print(TFile)
    matchT=re.search(r"T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",TFile)
    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))
        TStrings.append("T"+matchT.group(1))


#sort T values
sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]
sortedTStrings=[TStrings[ind] for ind in sortedInds]


def parseSummary(oneTFolder,obs_name):

    startingFileInd=-1
    startingVecPosition=-1
    lag=-1
    smrFile=oneTFolder+"/summary_"+obs_name+".txt"
    summaryFileExists=os.path.isfile(smrFile)
    if summaryFileExists==False:
        return startingFileInd,startingVecPosition,-1

    with open(smrFile,"r") as fptr:
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

        #match lag
        matchLag=re.search(r"lag=(\d+)",oneLine)
        if matchLag:
            lag=int(matchLag.group(1))
    return startingFileInd, startingVecPosition,lag



def sort_data_files_by_lpEnd(oneTFolder,obs_name,varName):
    """

    :param oneTFolder: Txxx
    :param obs_name: data files sorted by loopEnd
    :return:
    """

    dataFolderName=oneTFolder+"/"+obs_name+"_dataFiles/"+varName+"/"
    dataFilesAll=[]
    loopEndAll=[]

    for oneDataFile in glob.glob(dataFolderName+"/*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"loopEnd(\d+)",oneDataFile)
        if matchEnd:
            loopEndAll.append(int(matchEnd.group(1)))


    endInds=np.argsort(loopEndAll)
    # loopStartSorted=[loopStartAll[i] for i in startInds]
    sortedDataFiles=[dataFilesAll[i] for i in endInds]

    return sortedDataFiles

def U_dist_data2csvForOneT(oneTFolder,oneTStr,startingFileInd,startingVecPosition,lag):
    TRoot=oneTFolder
    sortedUDataFilesToRead=sort_data_files_by_lpEnd(TRoot,obs_U_dist,"U")
    startingUFileName=sortedUDataFilesToRead[startingFileInd]

    with open(startingUFileName,"rb") as fptr:
        inUStart=pickle.load(fptr)
    # print(len(sortedUDataFilesToRead))
    # print(startingUFileName)
    UVec=inUStart[startingVecPosition:]
    for pkl_file in sortedUDataFilesToRead[(startingFileInd+1):]:
        with open(pkl_file,"rb") as fptr:
            in_UArr=pickle.load(fptr)
            UVec=np.append(UVec,in_UArr)

    UVecSelected=UVec[::lag]
    # print("len(UVecSelected)="+str(len(UVecSelected)))
    dataArraySelected=UVecSelected
    #extract xA
    for j in range(0,unitCellNum):
        varName="xA"+str(j)
        sorted_xADataFilesToRead=sort_data_files_by_lpEnd(TRoot,obs_U_dist,varName)
        starting_xAFileName=sorted_xADataFilesToRead[startingFileInd]
        # print(len(sorted_xADataFilesToRead))
        # print(starting_xAFileName)
        with open(starting_xAFileName,"rb") as fptr:
            in_xAjStart=pickle.load(fptr)

        xAjVec=in_xAjStart[startingVecPosition:]
        for pkl_file in sorted_xADataFilesToRead[(startingFileInd+1):]:
            with open(pkl_file,"rb") as fptr:
                in_xAjArr=pickle.load(fptr)
                xAjVec=np.append(xAjVec,in_xAjArr)
        xAjVecSelected=xAjVec[::lag]
        # print("len(xAjVecSelected)="+str(len(xAjVecSelected)))
        dataArraySelected=np.vstack((dataArraySelected,xAjVecSelected))

    #extract xB
    for j in range(0,unitCellNum):
        varName="xB"+str(j)
        sorted_xBDataFilesToRead=sort_data_files_by_lpEnd(TRoot,obs_U_dist,varName)
        starting_xBFileName=sorted_xBDataFilesToRead[startingFileInd]
        with open(starting_xBFileName,"rb") as fptr:
            in_xBjStart=pickle.load(fptr)
        xBjVec=in_xBjStart[startingVecPosition:]
        for pkl_file in sorted_xBDataFilesToRead[(startingFileInd+1):]:
            with open(pkl_file,"rb") as fptr:
                in_xBjArr=pickle.load(fptr)
                xBjVec=np.append(xBjVec,in_xBjArr)

        xBjVecSelected=xBjVec[::lag]
        # print("len(xBjVecSelected)="+str(len(xBjVecSelected)))
        dataArraySelected=np.vstack((dataArraySelected,xBjVecSelected))

    colNamesAll=["U"]
    for j in range(0,unitCellNum):
        colNamesAll.append("xA"+str(j))
    for j in range(0,unitCellNum):
        colNamesAll.append("xB"+str(j))
    # col_names=",".join(colNamesAll)

    outCsvDataRoot=rowDirRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"+obs_U_dist+"/"
    Path(outCsvFolder).mkdir(parents=True, exist_ok=True)
    outCsvFile=outCsvFolder+"/"+obs_U_dist+"Data.csv"
    dfToSave=pd.DataFrame(dataArraySelected.T,columns=colNamesAll)
    dfToSave.to_csv(outCsvFile,index=False)





for k in range(0,len(sortedTFiles)):
    tStart=datetime.now()
    oneTFolder=sortedTFiles[k]
    oneTStr=sortedTStrings[k]

    startingfileIndTmp,startingVecIndTmp,lagTmp=parseSummary(oneTFolder,obs_U_dist)
    if startingfileIndTmp<0:
        print("summary file does not exist for "+oneTStr+" "+obs_U_dist)
        continue

    U_dist_data2csvForOneT(oneTFolder,oneTStr,startingfileIndTmp,startingVecIndTmp,lagTmp)
    tEnd=datetime.now()
    print("processed T="+str(sortedTVals[k])+": ",tEnd-tStart)


