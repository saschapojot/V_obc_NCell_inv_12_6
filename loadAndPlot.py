import numpy as np

import glob
from decimal import Decimal
import pickle
import re
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.signal import find_peaks
#This script loads  and plots the data in the last few files
def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

T=float(sys.argv[1])
unitCellNum=int(sys.argv[2])

TStr=format_using_decimal(T)
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

dataPath="./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/U_dist_dataFiles/"


plt_nameU="U"

inUPath=dataPath+"/"+plt_nameU+"/"

sorted_inUFiles=sort_data_files_by_lpEnd(inUPath)

lastFilesNum=30

files2Plot=sorted_inUFiles[-lastFilesNum:]

arrU=np.array([])
for pkl_file in files2Plot:
    with open(pkl_file,"rb") as fptr:
        arrIn=pickle.load(fptr)
        arrU=np.append(arrU,arrIn)


plt.figure()
plt.scatter(range(0,len(arrU)),arrU,s=1)
plt.title("U in last "+str(lastFilesNum)+"files")
plt.savefig("lastFilesU.png")

# peaksU, _ = find_peaks(arrU)
# periods = np.diff(peaksU)
# average_period = np.mean(periods)
#
# # Output the results
# print(f"Periods between peaks: {periods}")
# print(f"Average period: {average_period}")
#
# print(np.sqrt(np.var(periods,ddof=1)))