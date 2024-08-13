import sys
import glob
import re
import json
from decimal import Decimal
import pandas as pd
import numpy as np
import subprocess


#this script loads previous data
numArgErr=4
valErr=5
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

potential_function_name=jsonDataFromConf["potential_function_name"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]
startingVecPosition=jsonFromSummary["startingVecPosition"]
N=int(jsonDataFromConf["unitCellNum"])

if N<=0:
    print("N="+str(N)+"<=0")
    exit(valErr)
#search and read U_dist files
#give arbitrary values to L, d0Vec, d1Vec without reading data
UInit=6

# y0Init=1
# z0Init=1
# y1Init=1
coordsAll=np.array(list(range(1,2*N+1)))*0.77
xA=coordsAll[0::2]
xB=coordsAll[1::2]

loopLastFile=-1

#search csv files
csvFileList=[]
loopEndAll=[]
for file in glob.glob(U_dist_dataDir+"/*.csv"):
    csvFileList.append(file)
    matchEnd=re.search(r"loopEnd(\d+)",file)
    if matchEnd:
        loopEndAll.append(int(matchEnd.group(1)))


def create_loadedJsonData(UVal,xAVec,xBVec,loopLastFileVal):
    """

    :param UVal:
    :param xAVec:
    :param xBVec:
    :param loopLastFileVal:
    :return:
    """
    initDataDict={
        "U":str(UVal),
        "xAVec":list(xAVec),
        "xBVec":list(xBVec),
        "loopLastFile":str(loopLastFileVal)
    }
    # print(initDataDict)
    return json.dumps(initDataDict)

#if no data found, return the arbitrary values

if len(csvFileList)==0:
    loadedJsonDataStr=create_loadedJsonData(UInit,xA,xB,loopLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found csv data
sortedEndInds=np.argsort(loopEndAll)
sortedLoopEnd=[loopEndAll[ind] for ind in sortedEndInds]
sortedCsvFileNames=[csvFileList[ind] for ind in sortedEndInds]
loopLastFile=sortedLoopEnd[-1]

lastFileName=sortedCsvFileNames[-1]

def get_last_row(csv_file):
    result = subprocess.run(['tail', '-n', '1', csv_file], stdout=subprocess.PIPE)
    last_row = result.stdout.decode('utf-8').strip()
    return last_row


csvLastRowStr=get_last_row(lastFileName)
valsInLastRow = [float(value) for value in csvLastRowStr.split(',')]

UInit=valsInLastRow[0]

xA=valsInLastRow[1:1+N]
xB=valsInLastRow[1+N:2*N+1]

loadedJsonDataStr=create_loadedJsonData(UInit,xA,xB,loopLastFile)
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)