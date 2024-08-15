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
rowNum=0#int(sys.argv[1])
unitCellNum=int(sys.argv[1])

csvDataFolderRoot="../dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/csvOutAll/"

inCsvFile="../V_inv_12_6Params.csv"

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
    meanU2=meanU**2

    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)
    UConfHalfLength=np.sqrt(varU/len(UVec))
    nbins=500
    fig=plt.figure()
    axU=fig.add_subplot()
    (n0,_,_)=axU.hist(UVec,bins=nbins)

    meanUStr=str(np.round(meanU,4))
    print("T="+str(TVal)+", E(U)="+meanUStr)
    sigmaUStr=str(np.round(sigmaU,4))

    axU.set_title("T="+str(TVal))
    axU.set_xlabel("$U$")
    axU.set_ylabel("#")
    xPosUText=(np.max(UVec)-np.min(UVec))*1/2+np.min(UVec)
    yPosUText=np.max(n0)*2/3
    axU.text(xPosUText,yPosUText,"mean="+meanUStr+"\nsd="+sigmaUStr)
    plt.axvline(x=meanU,color="red",label="mean")
    axU.text(meanU*1.1,0.5*np.max(n0),str(meanU)+"$\\pm$"+str(sigmaU),color="red")
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
    LVec=xB_array[:,-1]-xA_array[:,0]
    LLength=len(LVec)
    LMean=np.mean(LVec)
    LVar=np.var(LVec,ddof=1)
    LConfHalfLength=np.sqrt(LVar/len(LVec))
    print("E(L)="+str(LMean))

    LVMean=np.mean(LVec*UVec)

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

    return [meanU,varU,UConfHalfLength,meanU2,LVMean,
            LMean,LVar,LConfHalfLength,
            d1Mean,d1ConfHalfInterval,
            d2Mean,d2ConfHalfInterval
            ]

UMeanValsAll=[]
U2MeanValsAll=[]# E*(V^{2})
UVarValsAll=[]
UConfHalfLengthAll=[]
LVMeanAll=[]

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
    meanU,varU,UConfHalfLength,meanU2,LVMean,LMean,LVar,LConfHalfLength,d1Mean,d1ConfHalfInterval,d2Mean,d2ConfHalfInterval=pltU_dist(oneTFile)

    UMeanValsAll.append(meanU)
    UVarValsAll.append(varU)
    UConfHalfLengthAll.append(UConfHalfLength)
    U2MeanValsAll.append(meanU2)
    LVMeanAll.append(LVMean)

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
U2MeanValsAll=np.array(U2MeanValsAll)
LVMeanAll=np.array(LVMeanAll)

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
ax.set_ylabel("E(U)")
ax.set_title("E(U)")
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/EU.png")
plt.close()
#######################################################

#######################################################
#plot var(U)

plt.figure()
plt.scatter(TToPlt,UVarValsAll[TInds],color="violet",label="mc")
# varVVals=[varV(T) for T in interpolatedTVals]
# plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
plt.title("var(U)")
plt.xlabel("$T$")
plt.legend(loc="best")
plt.savefig(csvDataFolderRoot+"/varU.png")
plt.close()


#######################################################
#######################################################
# C



#
# CValsAll=[]
# for j in range(0,len(UMeanValsAll)):
#     TTmp=sortedTVals[j]
#     CValsAll.append((UVarValsAll[j])/TTmp**2)
#
# CValsAll=np.array(CValsAll)
# #plot C
# plt.figure()
# plt.scatter(TToPlt,CValsAll[TInds]/unitCellNum,color="violet",label="mc")
# # varVVals=[varV(T) for T in interpolatedTVals]
# # plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
# plt.title("C per unit cell, unit cell number="+str(unitCellNum))
# plt.xlabel("$T$")
# plt.legend(loc="best")
# plt.savefig(csvDataFolderRoot+"/bkp.C.png")
# plt.close()

#######################################################

#######################################################
#alpha

# LVDiffTmp=LVMeanAll-LMeanValsAll*UMeanValsAll
# alphaValsAll=[]
# for j in range(0,len(sortedTVals)):
#     TTmp=sortedTVals[j]
#     alphaValsAll.append(LVDiffTmp[j]/(TTmp**2*LMeanValsAll[j]))
# alphaValsAll=np.array(alphaValsAll)/unitCellNum
#
# #plot alpha
# # print(TToPlt)
# # print(alphaValsAll[TInds])
# plt.figure()
# plt.scatter(TToPlt,alphaValsAll[TInds],color="violet",label="mc")
# # varVVals=[varV(T) for T in interpolatedTVals]
# # plt.plot(interpolatedTVals,varVVals,color="navy",label="theory")
# plt.title("alpha per unit cell, unit cell number="+str(unitCellNum))
# plt.xlabel("$T$")
# plt.legend(loc="best")
# plt.savefig(csvDataFolderRoot+"/bkp.alpha.png")
# plt.close()
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



#######################################################
#plot d1, d2 using lattice

d1d2InterleavedArray=[]
# print(d1MeanVecsAll)
# print(d2MeanVecsAll)
for n in range(0,len(d1MeanVecsAll)):
    d1VecOneTemp=d1MeanVecsAll[n]
    d2VecOneTemp=d2MeanVecsAll[n]
    # d1VecOneTemp=np.abs(d1VecOneTemp)
    # d2VecOneTemp=np.abs(d2VecOneTemp)
    rst=[]
    for i ,(x,y) in enumerate(zip(d1VecOneTemp,d2VecOneTemp)):
        rst.append(x)
        rst.append(y)
    rst.extend(d1VecOneTemp[len(d2VecOneTemp):])
    rst=np.array(rst)
    rst=np.insert(rst,0,0)

    d1d2InterleavedArray.append(rst)
d1d2InterleavedArray=np.array(d1d2InterleavedArray)

positions=np.cumsum(d1d2InterleavedArray,axis=1)
xANames=["x"+str(i)+"A" for i in  range(0,unitCellNum)]
xBNames=["x"+str(i)+"B" for i in  range(0,unitCellNum)]
xNames=[]
for i,(x,y) in enumerate(zip(xANames,xBNames)):
    xNames.append(x)
    xNames.append(y)
xNames.extend(xANames[len(xBNames):])
# print(len(sortedTVals))
# fig,axes=plt.subplots(len(sortedTVals),1,figsize=(30,10))
# for n in  range(0,len(sortedTVals)):
#     positionOneTemp=positions[n,:]
#     for i , pos in enumerate(positionOneTemp):
#         color='red' if i % 2 == 0 else 'blue'
#         axes[n].plot(pos,0,"o",color=color)
#     axes[n].set_xticks(positionOneTemp)
#     axes[n].set_xticklabels(xNames)
#     # Hide the y-axis
#     axes[n].get_yaxis().set_visible(False)
#     # Draw a horizontal line representing the axis
#     axes[n].axhline(0, color='black', linewidth=0.5)
#
#     d1d2InterTmp=d1d2InterleavedArray[n,:]
#     textPos=[]
#     for j in range(0,len(positionOneTemp)-1):
#         posTmp=(positionOneTemp[j+1]-positionOneTemp[j])/2+positionOneTemp[j]
#         textPos.append(posTmp)
#
#     for j,val in enumerate(textPos):
#         axes[n].text(val+0.1,-0.01,np.round(d1d2InterTmp[j+1],2),ha='center', va='bottom', fontsize=12)
#
#     axes[n].set_title("T="+str(sortedTVals[n]))
#
#
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# fig.suptitle('OBC, average distances')
# plt.savefig(csvDataFolderRoot+"/lattice.png")

def plot_positions(ax, positions, x_names, d1d2_interleaved_array, title):
    # Plot positions with alternating colors
    for i, pos in enumerate(positions):
        color = 'red' if i % 2 == 0 else 'blue'
        ax.plot(pos, 0, "o", color=color)

    # Set x-ticks and labels
    ax.set_xticks(positions)
    ax.set_xticklabels(x_names)

    # Hide the y-axis and draw a horizontal line
    ax.get_yaxis().set_visible(False)
    ax.axhline(0, color='black', linewidth=0.5)

    # Calculate and annotate text positions
    text_positions = [(positions[j+1] + positions[j]) / 2 for j in range(len(positions) - 1)]
    for j, val in enumerate(text_positions):
        ax.text(val + 0.1, -0.01, np.round(d1d2_interleaved_array[j+1], 2),
                ha='center', va='bottom', fontsize=12)

    # Set the subplot title
    ax.set_title(title)
width=40
if len(sortedTVals) == 1:
    fig, ax = plt.subplots(figsize=(width, 10))
    plot_positions(ax, positions[0, :], xNames, d1d2InterleavedArray[0, :], f"T={sortedTVals[0]}")
else:
    fig, axes = plt.subplots(len(sortedTVals), 1, figsize=(width, 10))
    for n in range(len(sortedTVals)):
        plot_positions(axes[n], positions[n, :], xNames, d1d2InterleavedArray[n, :], f"T={sortedTVals[n]}")

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.suptitle('OBC, average distances')
plt.savefig(csvDataFolderRoot + "/lattice.png")
#######################################################
