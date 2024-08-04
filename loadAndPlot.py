import numpy as np

import glob

import re
import matplotlib.pyplot as plt
import pandas as pd


T=5
unitCellNum=2


dataPath="./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+str(0.5)+"/U_dist_dataFiles/"

csvFileNamesAll=[]
loopEndAll=[]
for file in glob.glob(dataPath+"/*.csv"):
    matchEnd=re.search(r"loopEnd(\d+)",file)
    if matchEnd:
        loopEndAll.append(int(matchEnd.group(1)))
        csvFileNamesAll.append(file)

sortedLpEndInds=np.argsort(loopEndAll)

sortedCsvFilesAll=[csvFileNamesAll[ind] for ind in sortedLpEndInds]
# print(sortedCsvFilesAll)

plt_name0="x1A"
plt_name1="x1B"

vec0=np.array([])
vec1=np.array([])
for file in sortedCsvFilesAll:
    dfTmp=pd.read_csv(file)
    colValsTmp0=np.array(dfTmp[plt_name0])
    colValsTmp1=np.array(dfTmp[plt_name1])
    vec0=np.r_[vec0,colValsTmp0]
    vec1=np.r_[vec1,colValsTmp1]



startingInd0=int(len(vec0)*1/4)
endInd0=int(len(vec0)*0.1/4)
plt_vec0=vec0[startingInd0:]
plt.figure()
plt.plot(plt_vec0,color="black")
plt.title(plt_name0)
plt.savefig("tmp0.png")
plt.close()


diff_plt_vec=vec1[startingInd0:]-vec0[startingInd0:]

plt.figure()
plt.plot(diff_plt_vec,color="blue")
plt.title(plt_name1+"-"+plt_name0)
plt.savefig("tmp1.png")
plt.close()