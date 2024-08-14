import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal



#this script reads from statsAll.txt and plot U, L, d1, d2





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

TRegex=r"T\s*=\s*"+TStr
def matchOneLine(line):
    matchTRegex=re.search(TRegex,line)
    if matchTRegex:
        matchU=re.search(r"mean_U\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        matchN=re.search(r"N\s*=\s*(\d+)",line)
        matchL=re.search(r"mean_L\s*=\s*([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)",line)
        if matchU and matchN and matchL:
            UMeanTmp=float(matchU.group(1))
            NTmp= int(matchN.group(1))
            LTmp=float(matchL.group(1))
            return [NTmp,UMeanTmp,LTmp]
        else:
            print("mismatch in "+str(line))
statsFolder="./statsAll/"
statsAllFile=statsFolder+ "/statsAll.txt"
contents=removeCommentsAndEmptyLines(statsAllFile)

NVec=[]
UMeanVec=[]
LMeanVec=[]
for line in contents:
    NTmp,UMeanTmp,LMeanTmp=matchOneLine(line)
    NVec.append(NTmp)
    UMeanVec.append(UMeanTmp)
    LMeanVec.append(LMeanTmp)


sorted_NInds=np.argsort(NVec)
sorted_NVec=[NVec[ind] for ind in sorted_NInds]
sorted_UMeanVec=[UMeanVec[ind] for ind in sorted_NInds]
sorted_LMeanVec=[LMeanVec[ind] for ind in sorted_NInds]

sorted_NVec=np.array(sorted_NVec)
sorted_UMeanVec=np.array(sorted_UMeanVec)
sorted_LMeanVec=np.array(sorted_LMeanVec)

sorted_UPerUnitCell=sorted_UMeanVec/sorted_NVec
sorted_LPerUnitCell=sorted_LMeanVec/sorted_NVec

plt.figure()
plt.plot(sorted_NVec,sorted_UPerUnitCell,label="T="+TStr,color="blue")
plt.title("E(U) per unit cell")
plt.xlabel("$T$")
plt.ylabel("$E(U)/N$")
plt.legend(loc="best")
plt.savefig(statsFolder+"UPerN.png")


plt.figure()
plt.plot(sorted_NVec,sorted_LPerUnitCell,label="T="+TStr,color="red")
plt.title("E(L) per unit cell")
plt.xlabel("$T$")
plt.ylabel("$E(L)/N$")
plt.legend(loc="best")
plt.savefig(statsFolder+"LPerN.png")