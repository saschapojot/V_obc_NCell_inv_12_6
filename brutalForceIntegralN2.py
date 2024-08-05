import numpy as np
import scipy.integrate as integrate
import pandas as pd
from datetime import datetime
import json
#this script benchmarks the mc computation by brutal force integration
#This script computes Z
rowNum=0
inParamFileName="./V_inv_12_6Params.csv"

inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])
T=2
shift=199

def V(y0, z0, y1):
    # L, y0, z0, y1 = x
    val=a1*y0**(-12)-b1*y0**(-6)\
        +a2*z0**(-12)-b2*z0**(-6)\
        + a1*y1**(-12)-b1*y1**(-6)\
        +shift


    return val

#minimizers of V
r1=(2*a1/b1)**(1/6)
r2=(2*a2/b2)**(1/6)

eps=(r1+r1)/20

y0Range=[r1-eps,r1+eps]

z0Range=[r2-eps,r2+eps]

y1Range=[r1-eps,r1+eps]



def Z(y0, z0, y1, beta):
    #
    return np.exp(-beta * V(y0, z0, y1))

tZStart=datetime.now()

beta = 1/T

result, error = integrate.nquad(lambda y0, z0, y1: Z(y0, z0, y1, beta), [y0Range, z0Range, y1Range])

print("Z Integral result:", result)
print("Z Estimated error:", error)
tZEnd=datetime.now()
print("Z time: ",tZEnd-tZStart)
#####################################
##E(V)

tVStart=datetime.now()
rstV,errV=integrate.nquad(lambda y0, z0, y1: V(y0, z0, y1)*Z(y0, z0, y1, beta), [y0Range, z0Range, y1Range])

print("V Integral result:", rstV)
print("V Estimated error:", errV)
tVEnd=datetime.now()
print("V time: ",tVEnd-tVStart)
EV=rstV/result-shift
print("E(V)="+str(EV))
#####################################
#####################################
## E(L)
tLStart=datetime.now()
rstL,errL=integrate.nquad(lambda  y0, z0, y1: (y0+z0+y1)*Z(y0, z0, y1, beta), [y0Range, z0Range, y1Range])
print("L Integral result:", rstL)
print("L Estimated error:", errL)
tLEnd=datetime.now()
print("L time: ",tLEnd-tLStart)
EL=rstL/result
print("E(L)="+str(EL))

#####################################

