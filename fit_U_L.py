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

#this script performs linear regression on U

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


def one_csv2_array(csvFileName,unitCellNum):

    df=pd.read_csv(csvFileName)

    #U
    UVec=np.array(df["U"])

    #L
    distArray=np.array(df.iloc[:,1:])
    nRow,nCol=distArray.shape
    xA_array=distArray[:,:unitCellNum]
    xB_array=distArray[:,unitCellNum:]
    LVec=xB_array[:,-1]-xA_array[:,0]
    LVec=np.array(LVec)

    return [UVec,LVec]

N_invToFit=[]
UToFit=[]
LToFit=[]

for j in range(0,len(sorted_csvFilesVec)):
    csvFileTmp=sorted_csvFilesVec[j]
    print("reading file "+str(csvFileTmp))
    NTmp=sorted_numUnitCellOfCsvVec[j]
    UVecTmp,LVecTmp=one_csv2_array(csvFileTmp,NTmp)
    N_invToFit.extend([1/NTmp]*len(UVecTmp))
    UToFit.extend(UVecTmp/NTmp)
    LToFit.extend(LVecTmp/NTmp)

N_invToFitWithAConst=sm.add_constant(N_invToFit)

modelU=sm.OLS(UToFit,N_invToFitWithAConst)
resultsU = modelU.fit()

# Print the detailed summary
print(resultsU.summary())

N_invUnique=np.unique(N_invToFit)
N_invUniqueWithConst=sm.add_constant(N_invUnique)

U_pred=resultsU.predict(N_invUniqueWithConst)


N_inv0=[0]
N_inv0WithConst=sm.add_constant(N_inv0)
U0Pred=resultsU.predict(N_inv0WithConst)
# print(U0Pred)

statsFolder="./statsAll/"
plt.figure()
plt.xlim(left=0,right=0.55)

plt.plot(N_invUnique,U_pred,color="red",label="linear fit",marker='o',markersize=3)
plt.scatter(N_inv0,U0Pred,color="blue",marker="D",s=30,label="Intercept")
intercept=U0Pred[0]
plt.text(0, intercept, f'{intercept:.2f}', color='blue', verticalalignment='center', horizontalalignment='right')
plt.xlabel(r"$\frac{1}{N}$")
plt.legend(loc="best")
plt.ylabel(r"E(U)/N")
plt.title("E(U) per unit cell, T="+TStr)
plt.savefig(statsFolder+"/fit_U.png")
plt.close()