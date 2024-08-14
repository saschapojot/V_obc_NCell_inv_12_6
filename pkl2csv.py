import pickle
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
import matplotlib.pyplot as plt
import glob
import re
from decimal import Decimal
#this script prints part of an array



def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
def sort_data_files_by_lpEnd(oneDataFolder):
    """

    :param oneDataFolder:
    :return:
    """


    dataFolderName=oneDataFolder
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


cellInd=26
unitCellNum=80
T=0.5
TStr=format_using_decimal(T)
dataRoot="./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/U_dist_dataFiles/"

print(dataRoot)
in_xAPath=dataRoot+"/xA"+str(cellInd)+"/"

in_xBPath=dataRoot+"/xB"+str(cellInd)+"/"

inUPath=dataRoot+"/U/"

sorted_inUFiles=sort_data_files_by_lpEnd(inUPath)
sorted_in_xAFiles=sort_data_files_by_lpEnd(in_xAPath)
sorted_in_xBFiles=sort_data_files_by_lpEnd(in_xBPath)

fileInd=-2
inFile1=sorted_in_xAFiles[fileInd]
inFile2=sorted_in_xBFiles[fileInd]
inFileU=sorted_inUFiles[fileInd]

with open(inFile1,"rb") as fptr:
    arr1=pickle.load(fptr)
arr1=np.array(arr1)
with open(inFile2,"rb") as fptr:
    arr2=pickle.load(fptr)
arr2=np.array(arr2)
# arr=np.reshape(arr,(-1,2*N+1))

inParamFileName="./V_inv_12_6Params.csv"
rowNum=0
inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])

def V1(r):
    return a1*r**(-12)-b1*r**(-6)
def auto_corrForOneColumn(colVec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    eps=1e-2
    NLags=int(len(colVec)*3/4)
    print("NLags="+str(NLags))
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
    try:
        acfOfVec=sm.tsa.acf(colVec,nlags=NLags)
    except Warning as w:
        same=True
    acfOfVecAbs=np.abs(acfOfVec)
    minAutc=np.min(acfOfVecAbs)
    print("minAutc="+str(minAutc))

    lagVal=-1
    if minAutc<=eps:
        lagVal=np.where(acfOfVecAbs<=eps)[0][0]
    # np.savetxt("autc.txt",acfOfVecAbs[lagVal:],delimiter=',')

    return same,lagVal

arr=arr2-arr1
same,lagVal=auto_corrForOneColumn(arr)
# print("lag="+str(lagVal))
outCsvName="./show.csv"

part=arr
V1Part=[V1(r) for r in part]
V1Part=np.array(V1Part)
V1Increment=V1Part[1:]-V1Part[:-1]
sortedV1Incre=np.sort(V1Increment)
# print(len(sortedV1Incre))
df=pd.DataFrame(part)

df.to_csv(outCsvName,index=False,header=None)
plt.figure()
plt.scatter(range(0,len(part)),part,s=1)
plt.title("dist")
plt.savefig("dist.png")
plt.close()

plt.figure()
plt.scatter(range(0,len(V1Part)),V1Part,s=1)
plt.title("U")
plt.savefig("U.png")
plt.close()

plt.figure()
plt.scatter(range(0,len(V1Increment)),V1Increment,s=1)
plt.title("V1 Increment")
plt.savefig("V1Inc.png")
plt.close()


print(inFileU)
with open(inFileU,"rb") as fptr:
    UVec=np.array(pickle.load(fptr))

ULength=int(1e6)
UPart=UVec[-ULength:]
sameU,lagU=auto_corrForOneColumn(UPart)
print("lagU="+str(lagU))
UDiff=UPart[1:]-UPart[:-1]

plt.figure()
plt.scatter(range(0,len(UPart)),UPart,s=1)
plt.title("UAll")

plt.savefig("UAll.png")
plt.close()
print(UPart)