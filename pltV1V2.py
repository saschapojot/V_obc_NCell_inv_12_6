import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle


inParamFileName="./V_inv_12_6Params.csv"

unitCellNum=20
rowNum=0
inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])


def V1(r):
    return a1*r**(-12)-b1*r**(-6)
r1Min=(2*a1/b1)**(1/6)

V1Min=V1(r1Min)

r=2.61
V1Val=V1(r)
print(V1Min)
print(V1Val)



