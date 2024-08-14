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

#This script fits the values of U for different N, same T

def format_using_decimal(value):
    # Convert the float to a Decimal using string conversion to avoid precision issues
    decimal_value = Decimal(str(value))
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)



if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()
T=float(sys.argv[1])


TStr=format_using_decimal(T)

rowNum=0
#search N files

dataFileRoot="./"
numUnitCellVec=[]
dataFolderVec=[]
for file in glob.glob(dataFileRoot+"/*"):
    matchDataAll=re.search(r"dataAllUnitCell(\d+)",file)
    if matchDataAll:
        dataFolderVec.append(file)
        numUnitCellVec.append(int(matchDataAll.group(1)))

numUnitCellInds=np.argsort(numUnitCellVec)

sortedNumUnitCellVec=[numUnitCellVec[ind] for ind in numUnitCellInds]
sortedDataFolderVec=[dataFolderVec[ind] for ind in numUnitCellInds]

sorted_csvFilesVec=[]#to stats
sorted_numUnitCellOfCsvVec=[]#to stats
sorted_statsFolder=[]#to stats
for j in range(0,len(sortedDataFolderVec)):
    folder=sortedDataFolderVec[j]
    statsFolderTmp=folder+"/row"+str(rowNum)+"/csvOutAll/T"+TStr+"/"
    csvFileTmp=statsFolderTmp+"/U_dist/U_distData.csv"
    if os.path.exists(csvFileTmp):
        sorted_csvFilesVec.append(csvFileTmp)
        sorted_numUnitCellOfCsvVec.append(sortedNumUnitCellVec[j])
        sorted_statsFolder.append(statsFolderTmp)

    else:
        print(csvFileTmp+" does not exist.")

# print(sorted_csvFilesVec)
# print(sorted_numUnitCellOfCsvVec)


def read_1_csv(csvFileName,unitCellNum,d1Ind, d2Ind):
    """

    :param csvFileName:
    :return: statistics of U and L, d1, d2
    """
    df=pd.read_csv(csvFileName)

    #U
    UVec=np.array(df["U"])
    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)
    meanU=np.mean(UVec)
    hf_U=np.sqrt(varU/len(UVec))

    #L
    distArray=np.array(df.iloc[:,1:])
    nRow,nCol=distArray.shape
    xA_array=distArray[:,:unitCellNum]
    xB_array=distArray[:,unitCellNum:]
    LVec=xB_array[:,-1]-xA_array[:,0]
    meanL=np.mean(LVec)
    varL=np.var(LVec,ddof=1)
    sigmaL=np.sqrt(varL)
    hf_L=np.sqrt(varL/len(LVec))


    #d1
    d1Vec=xB_array[:,d1Ind]-xA_array[:,d1Ind]
    mean_d1=np.mean(d1Vec)
    var_d1=np.var(d1Vec,ddof=1)
    sigma_d1=np.sqrt(var_d1)
    hf_d1=np.sqrt(var_d1/len(d1Vec))

    #d2
    d2Vec=xA_array[:,d2Ind+1]-xB_array[:,d2Ind]
    mean_d2=np.mean(d2Vec)
    var_d2=np.var(d2Vec,ddof=1)
    sigma_d2=np.sqrt(var_d2)
    hf_d2=np.sqrt(var_d2/len(d2Vec))

    return [meanU,sigmaU,hf_U,
            meanL,sigmaL,hf_L,
            mean_d1,sigma_d1,hf_d1,
            mean_d2,sigma_d2,hf_d2]


def contentsFile_U_L(csvFileName,unitCellNum,d1Ind, d2Ind,statsFolder):
    """

    :param retData: data from read_1_csv()
    :return:
    """
    retData=read_1_csv(csvFileName,unitCellNum,d1Ind, d2Ind)

    meanU,sigmaU,hf_U,\
    meanL,sigmaL,hf_L,\
    mean_d1,sigma_d1,hf_d1,\
    mean_d2,sigma_d2,hf_d2=retData

    outStatsFile=statsFolder+"/stats.txt"

    contents=[
        "#This is a file of statistics of U,L d1_"+str(d1Ind)+", d2_"+str(d2Ind)+"\n"
        "\n",
        "T="+TStr+", N="+str(unitCellNum)+", mean_U="+str(meanU)+", "\
        +"sigma_U="+str(sigmaU)+", hf_U="+str(hf_U)+", "\
        +"mean_L="+str(meanL)+", sigma_L="+str(sigmaL)+", hf_L="+str(hf_L)+", "\
        +"mean_d1_"+str(d1Ind)+"="+str(mean_d1)+", sigma_d1_"+str(d1Ind)+"="+str(sigma_d1)+", "\
        +"hf_d1_"+str(d1Ind)+"="+str(hf_d1)+", "\
        +"mean_d2_"+str(d2Ind)+"="+str(mean_d2)+", sigma_d2_"+str(d2Ind)+"="+str(sigma_d2)+", "\
        +"hf_d2_"+str(hf_d2)+"\n"
    ]

    with open(outStatsFile,"w+") as fptr:
        fptr.writelines(contents)

d1Ind=1
d2Ind=0

tContentStart=datetime.now()
for j in range(0,len(sorted_csvFilesVec)):

    csvFileTmp=sorted_csvFilesVec[j]
    print("reading file "+str(csvFileTmp))
    NTmp=sorted_numUnitCellOfCsvVec[j]
    statsFolderTmp=sorted_statsFolder[j]
    contentsFile_U_L(csvFileTmp,NTmp,d1Ind,d2Ind,statsFolderTmp)

tContentEnd=datetime.now()

print("stats time: ",tContentEnd-tContentStart)