import pickle
import numpy as np
import pandas as pd
import statsmodels.api as sm
import warnings
#this script prints part of an array

inFile="./dataAllUnitCell20/row0/T1/U_dist_dataFiles/U/loopStart9000000loopEnd9999999.U.pkl"


with open(inFile,"rb") as fptr:
    arr=pickle.load(fptr)

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


auto_corrForOneColumn(arr)

outCsvName="./show.csv"

part=arr

df=pd.DataFrame(part)

df.to_csv(outCsvName,index=False,header=None)
