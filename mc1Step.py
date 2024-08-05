import numpy as np
import pandas as pd

import scipy.integrate as integrate

inCsvFile="./dataAllUnitCell5/row0/T0.5/U_dist_dataFiles/loopStart4700000loopEnd4799999.U_dist.csv"
inParamFileName="./V_inv_12_6Params.csv"



inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[0,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])

N=5


df=pd.read_csv(inCsvFile)

row1=df.iloc[84555,:]

row2=df.iloc[84556,:]

xAVec1=np.array(row1[1:(N+1)])
xBVec1=np.array(row1[(N+1):])


xAVec2=np.array(row2[1:(N+1)])
xBVec2=np.array(row2[(N+1):])

def V1(r):
    val=a1*r**(-12)-b1*r**(-6)
    return val


def V2(r):
    val=a2*r**(-12)-b2*r**(-6)
    return val



def V(xA,xB):
    unitCellNum=len(xA)

    val=0
    for j in range(0,unitCellNum):
        val+=V1(xB[j]-xA[j])
    for j in range(0,unitCellNum-1):
        val+=V2(xA[j+1]-xB[j])
    return val




d1Vec1=[]
d1Vec2=[]
for j in range(0,N):
    d1Vec1.append(xBVec1[j]-xAVec1[j])
for j in range(0,N):
    d1Vec2.append(xBVec2[j]-xAVec2[j])
V1Vec1=[V1(r) for r in d1Vec1]

V1Vec2=[V1(r) for r in d1Vec2]




h=0.013
x1=0.3105046680000000
x2=0.3257084530000000


lm=13.4949

def f1(y):
    val=np.exp(-1/(2*h**2)*(y-x1)**2)
    return val


def f2(y):
    val=np.exp(-1/(2*h**2)*(y-x2)**2)
    return val


result1, error1 = integrate.quad(f1, 0, lm)
print(f"Result1: {result1}")
print(f"Estimated error1: {error1}")