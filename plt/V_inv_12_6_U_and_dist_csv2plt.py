import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd



#This script loads csv data and plot

if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()
rowNum=int(sys.argv[1])
unitCellNum=2

csvDataFolderRoot="../dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/csvOutAll/"

inCsvFile="../V2Params.csv"

TVals=[]
TFileNames=[]

for TFile in glob.glob(csvDataFolderRoot+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]
def pltU_dist(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: U plots, U mean, U var, dist plots, dist mean, dist var
    """
    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))

    U_distPath=oneTFile+"/U_dist/U_distData.csv"

    df=pd.read_csv(U_distPath)

    UVec=np.array(df.iloc[:,0])


    print("T="+str(TVal)+", data num="+str(len(UVec)))

    #U part
    meanU=np.mean(UVec)

    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)
    UConfHalfLength=np.sqrt(varU/len(UVec))
    nbins=500
    fig=plt.figure()
    axU=fig.add_subplot()
    (n0,_,_)=axU.hist(UVec,bins=nbins)
    meanU=np.round(meanU,4)
    print("T="+str(TVal)+", E(U)="+str(meanU))
    sigmaU=np.round(sigmaU,4)

    axU.set_title("T="+str(TVal))
    axU.set_xlabel("$U$")
    axU.set_ylabel("#")
    xPosUText=(np.max(UVec)-np.min(UVec))*1/2+np.min(UVec)
    yPosUText=np.max(n0)*2/3
    axU.text(xPosUText,yPosUText,"mean="+str(meanU)+"\nsd="+str(sigmaU))
    plt.axvline(x=meanU,color="red",label="mean")
    axU.text(meanU*1.1,0.5*np.max(n0),str(meanU)+"$\pm$"+str(sigmaU),color="red")
    axU.hlines(y=0,xmin=meanU-sigmaU,xmax=meanU+sigmaU,color="green",linewidth=15)

    plt.legend(loc="best")

    EHistOut="T"+str(TVal)+"UHist.png"
    plt.savefig(oneTFile+"/"+EHistOut)

    plt.close()

    ### test normal distribution for mean U
    #block mean
    USelectedAll=UVec

    def meanPerBlock(length):
        blockNum=int(np.floor(len(USelectedAll)/length))
        UMeanBlock=[]
        for blkNum in range(0,blockNum):
            blkU=USelectedAll[blkNum*length:(blkNum+1)*length]
            UMeanBlock.append(np.mean(blkU))
        return UMeanBlock

    fig=plt.figure(figsize=(20,20))
    fig.tight_layout(pad=5.0)
    lengthVals=[2,5,7,10]
    for i in range(0,len(lengthVals)):
        l=lengthVals[i]
        UMeanBlk=meanPerBlock(l)
        ax=fig.add_subplot(2,2,i+1)
        (n,_,_)=ax.hist(UMeanBlk,bins=100,color="aqua")
        xPosTextBlk=(np.max(UMeanBlk)-np.min(UMeanBlk))*1/7+np.min(UMeanBlk)
        yPosTextBlk=np.max(n)*3/4
        meanTmp=np.mean(UMeanBlk)
        meanTmp=np.round(meanTmp,3)
        sdTmp=np.sqrt(np.var(UMeanBlk))
        sdTmp=np.round(sdTmp,3)
        ax.set_title("Bin Length="+str(l))
        ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
    fig.suptitle("T="+str(TVal))
    plt.savefig(oneTFile+"/T"+str(TVal)+"UBlk.png")

    #dist part
    distArray=np.array(df.iloc[:,1:])
    nRow,nCol=distArray.shape
    xA_array=distArray[:,:unitCellNum]
    xB_array=distArray[:,unitCellNum:]
    LVec=xB_array[:,-1]-xA_array[:,-1]
    LLength=len(LVec)
    LMean=np.mean(LVec)
    LVar=np.var(LVec,ddof=1)
    LConfHalfLength=np.sqrt(LVar/len(LVec))
    print("E(L)="+str(LMean))

    d1Array=np.zeros((nRow,unitCellNum),dtype=float)
    d2Array=np.zeros((nRow,unitCellNum-1),dtype=float)

    for j in range(0,unitCellNum):
        d1Array[:,j]=xB_array[:,j]-xA_array[:,j]

    for j in range(0,unitCellNum-1):
        d2Array[:,j]=xA_array[:,j+1]-xB_array[:,j]


    d1Mean=np.mean(d1Array,axis=0)
    d2Mean=np.mean(d2Array,axis=0)




    d1Var=np.var(d1Array,axis=0,ddof=1)
    d2Var=np.var(d2Array,axis=0,ddof=1)


    d1ConfHalfInterval=np.sqrt(d1Var/LLength)
    d2ConfHalfInterval=np.sqrt(d2Var/LLength)

    return [meanU,varU,UConfHalfLength,
            LMean,LVar,LConfHalfLength,
            d1Mean,d1ConfHalfInterval,
            d2Mean,d2ConfHalfInterval
            ]

UMeanValsAll=[]
UVarValsAll=[]
UConfHalfLengthAll=[]

LMeanValsAll=[]
LVarValsAll=[]
LConfHalfLengthAll=[]



d1MeanVecsAll=[]
d1HalfIntervalAll=[]

d2MeanVecsAll=[]
d2HalfIntervalAll=[]

tStatsStart=datetime.now()

for k in range(0,len(sortedTFiles)):
    oneTFile=sortedTFiles[k]
    meanU,varU,UConfHalfLength,LMean,LVar,LConfHalfLength,d1Mean,d1ConfHalfInterval,d2Mean,d2ConfHalfInterval=pltU_dist(oneTFile)

    UMeanValsAll.append(meanU)
    UVarValsAll.append(varU)
    UConfHalfLengthAll.append(UConfHalfLength)

    LMeanValsAll.append(LMean)
    LVarValsAll.append(LVar)
    LConfHalfLengthAll.append(LConfHalfLength)



    d1MeanVecsAll.append(d1Mean)
    d1HalfIntervalAll.append(d1ConfHalfInterval)

    d2MeanVecsAll.append(d2Mean)
    d2HalfIntervalAll.append(d2ConfHalfInterval)


UMeanValsAll=np.array(UMeanValsAll)
UVarValsAll=np.array(UVarValsAll)
UConfHalfLengthAll=np.array(UConfHalfLengthAll)

LMeanValsAll=np.array(LMeanValsAll)
LVarValsAll=np.array(LVarValsAll)
LConfHalfLengthAll=np.array(LConfHalfLengthAll)



d1MeanVecsAll=np.array(d1MeanVecsAll)
d1HalfIntervalAll=np.array(d1HalfIntervalAll)

d2MeanVecsAll=np.array(d2MeanVecsAll)
d2HalfIntervalAll=np.array(d2HalfIntervalAll)

tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)
sortedTVals=np.array(sortedTVals)

TInds=np.where(sortedTVals<100)
TToPlt=sortedTVals[TInds]

######################################################
#plt E(U)
fig, ax = plt.subplots()
ax.errorbar(TToPlt,UMeanValsAll[TInds],yerr=UConfHalfLengthAll,fmt='o',color="black", ecolor='r', capsize=5,label='mc')
# EVVals=[EV(T) for T in interpolatedTVals]
# ax.plot(interpolatedTVals,EVVals,color="green",label="theory")
ax.set_xlabel('$T$')
ax.set_ylabel("E(V)")
ax.set_title("E(V)")
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/EV.png")
plt.close()
#######################################################

#######################################################
#plot var(U)

plt.figure()
plt.scatter(TToPlt,UVarValsAll[TInds],color="violet",label="mc")
# varVVals=[varV(T) for T in interpolatedTVals]
# plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
plt.title("var(V)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/varV.png")
plt.close()


#######################################################


#######################################################
#d1 vals
d1ToPlt=d1MeanVecsAll.T
d1ToPltNRow,d1ToPltNCol=d1ToPlt.shape
d1HalfInterToPlt=d1HalfIntervalAll.T

fig,axes=plt.subplots(d1ToPltNRow,1,figsize=(10,10))
for i in range(0,d1ToPltNRow):
    axes[i].errorbar(TToPlt, d1ToPlt[i,:], yerr=d1HalfInterToPlt[i,:], fmt='o', color="black", ecolor='r', capsize=5, label="d"+str(i)+"A"+str(i)+"B")
    axes[i].set_xlabel('$T$')
    axes[i].set_ylabel('E('+"d"+str(i)+"A"+str(i)+"B"+")")
    axes[i].set_title("Intracell distance between "+str(i)+"A and "+str(i)+"B")
    axes[i].legend()

plt.savefig(csvDataFolderRoot+"/d1.png")
plt.close()
#######################################################

#######################################################
# E(L)
fig, ax = plt.subplots()
ax.errorbar(TToPlt,LMeanValsAll[TInds],yerr=LConfHalfLengthAll,fmt='o',color="black", ecolor='r', capsize=5,label='mc')
ax.set_xlabel('$T$')
ax.set_ylabel("E(L)")
ax.set_title("E(L)")
# ELVals=[EL(T) for T in interpolatedTVals]
# ELVals=np.array(ELVals)
# ax.plot(interpolatedTVals,ELVals,color="blue",label="theory")
ax.legend()
plt.tight_layout()
plt.savefig(csvDataFolderRoot+"/EL.png")
plt.close()
#######################################################
#######################################################
# var(L)
plt.figure()
plt.scatter(TToPlt,LVarValsAll[TInds],color="magenta",label="mc")

# varLVals=[varL(T) for T in interpolatedTVals]
# plt.plot(interpolatedTVals,varLVals,color="green",label="theory")
plt.title("var(L)")
plt.ylabel("var(L)")
plt.xlabel("$T$")
plt.legend(loc="best")
# plt.ylim((0,5.5))
plt.savefig(csvDataFolderRoot+"/varL.png")
plt.close()

#######################################################
