import pickle
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
import matplotlib.pyplot as plt
#this script prints part of an array

inFile1="./dataAllUnitCell14/row0/T0.5/U_dist_dataFiles/xB6/loopStart9000000loopEnd9999999.xB6.pkl"
inFile2="./dataAllUnitCell14/row0/T0.5/U_dist_dataFiles/xA6/loopStart9000000loopEnd9999999.xA6.pkl"

with open(inFile1,"rb") as fptr:
    arr1=pickle.load(fptr)
arr1=np.array(arr1)
with open(inFile2,"rb") as fptr:
    arr2=pickle.load(fptr)
arr2=np.array(arr2)
# arr=np.reshape(arr,(-1,2*N+1))
def auto_corrForOneColumn(colVec):
    """

    :param colVec: a vector of data
    :return:
    """
    same=False
    eps=1e-2
    NLags=int(len(colVec)*1/10)
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
auto_corrForOneColumn(arr)

outCsvName="./show.csv"

part=arr

df=pd.DataFrame(part)

df.to_csv(outCsvName,index=False,header=None)
plt.scatter(range(0,len(arr)),arr,s=0.1)
plt.savefig("dist.png")
