import numpy as np

import glob

import re
import matplotlib.pyplot as plt
import pandas as pd


T=7
unitCellNum=5


dataPath="./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+str(T)+"/U_dist_dataFiles/"

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
plt_name1="x0B"
plt_name2="x4B"

vec0=np.array([])
vec1=np.array([])
vec2=np.array([])
for file in sortedCsvFilesAll:
    dfTmp=pd.read_csv(file)
    colValsTmp0=np.array(dfTmp[plt_name0])
    colValsTmp1=np.array(dfTmp[plt_name1])
    colValsTmp2=np.array(dfTmp[plt_name2])
    vec0=np.r_[vec0,colValsTmp0]
    vec1=np.r_[vec1,colValsTmp1]
    vec2=np.r_[vec2,colValsTmp2]



startingInd0=int(len(vec0)*3/4)
endInd0=int(len(vec0)*4/4)
plt_vec0=vec0[startingInd0:]
plt.figure()
plt.scatter(range(0,len(plt_vec0)),plt_vec0,color="black",s=0.05)
plt.title(plt_name0+", T="+str(T))
plt.savefig("tmp0.png")
plt.close()


diff_plt_vec=vec1[startingInd0:endInd0]-vec0[startingInd0:endInd0]

plt.figure()
plt.scatter(range(0,len(diff_plt_vec)),diff_plt_vec,color="blue",s=0.05)
plt.title(plt_name1+"-"+plt_name0+", T="+str(T))
plt.savefig("tmp1.png")
plt.close()

# indM=np.where(diff_plt_vec>=0)[0][0]
# print(indM)

plt_vec2=vec2[startingInd0:]
plt.figure()
plt.scatter(range(0,len(plt_vec2)),plt_vec2,color="green",s=0.05)
plt.title(plt_name2+", T="+str(T))
plt.savefig("tmp2.png")
plt.close()
print(np.max(plt_vec2))