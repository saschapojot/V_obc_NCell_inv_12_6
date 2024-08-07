import numpy as np
from datetime import datetime
import sys
import re
import glob
import os
import json
from pathlib import Path
import pandas as pd

#this script extracts effective data

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()


rowNum=int(sys.argv[1])

unitCellNum=5
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



def sort_data_files_by_lpEnd(oneTFolder,obs_name):
    """

    :param oneTFolder: Txxx
    :param obs_name: data files sorted by loopEnd
    :return:
    """

    dataFolderName=oneTFolder+"/"+obs_name+"_dataFiles/"
    dataFilesAll=[]
    loopEndAll=[]

    for oneDataFile in glob.glob(dataFolderName+"/*.csv"):
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
    sortedDataFilesToRead=sort_data_files_by_lpEnd(TRoot,obs_U_dist)
    startingFileName=sortedDataFilesToRead[startingFileInd]
    # print("startingFileInd="+str(startingFileInd))
    # print("startingVecPosition="+str(startingVecPosition))
    #read the starting U_dist csv file
    in_dfStart=pd.read_csv(startingFileName)

    dataArray=np.array(in_dfStart,dtype=float)

    #read the rest of the csv files
    for csv_file in sortedDataFilesToRead[(startingFileInd+1):]:
        in_df=pd.read_csv(csv_file)
        dataArray=np.r_[dataArray,in_df]

    dataArraySelected=dataArray[::lag,:]

    outCsvDataRoot=rowDirRoot+"/csvOutAll/"
    outCsvFolder=outCsvDataRoot+"/"+oneTStr+"/"+obs_U_dist+"/"
    Path(outCsvFolder).mkdir(parents=True, exist_ok=True)
    outCsvFile=outCsvFolder+"/"+obs_U_dist+"Data.csv"

    col_names = in_dfStart.columns

    dfToSave=pd.DataFrame(dataArraySelected,columns=col_names)
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