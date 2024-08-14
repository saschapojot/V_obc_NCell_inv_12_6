import re
import sys
import json
import os
from decimal import Decimal
import numpy as np
import glob
from pathlib import Path
#this script combines all stats files
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


T=0.5

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

sorted_statsFilesVec=[]
for j in range(0,len(sortedDataFolderVec)):
    folder=sortedDataFolderVec[j]
    statsFileTmp=folder+"/row"+str(rowNum)+"/csvOutAll/T"+TStr+"/stats.txt"
    if os.path.exists(statsFileTmp):
        sorted_statsFilesVec.append(statsFileTmp)
    else:
        print(statsFileTmp+" does not exist.")


outStatsFolder="./statsAll/"
Path(outStatsFolder).mkdir(parents=True, exist_ok=True)
outStatsFile=outStatsFolder+"/statsAll.txt"
contentsAll=[]
for file in sorted_statsFilesVec:
    contentsTmp=removeCommentsAndEmptyLines(file)[0]+"\n"
    # print(contentsTmp)
    contentsAll.append(contentsTmp)

with open(outStatsFile,"w+") as fptr:
    fptr.writelines(contentsAll)

