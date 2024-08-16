


#when checking equilibrium, if lag is large, it may be because the equilibrium is far
#or, the data are around equilibrium, but h is too small
#for the lowest temperature, the value of h needs some trial and run to guess, once this value is obtained,
#for higher temperatures, one increases the h value
#also, one may change the number of files to check in statictics

#to diagnose the data, run :
python pkl2csv.py T N
###stats

1. python oneStat_U_L.py T
2. python combinedStats.py
3. python stats2_plt.py

###fit
1. python fit_U_L.py T  #for all T
2. python fit_ULCombined.py