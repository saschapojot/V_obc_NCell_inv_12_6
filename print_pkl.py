import pickle
import numpy as np
import pandas as pd

#this script prints part of an array

inFile="./dataAllUnitCell20/row0/T1/U_dist_dataFiles/U/loopStart9000000loopEnd9999999.U.pkl"
N=10

with open(inFile,"rb") as fptr:
    arr=pickle.load(fptr)

# arr=np.reshape(arr,(-1,2*N+1))

outCsvName="./show.csv"

part=arr

df=pd.DataFrame(part)

df.to_csv(outCsvName,index=False,header=None)
