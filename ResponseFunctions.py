from LatticeHandler import *
from ResultsList import *

#Script to get and plot the response functions for the desired temperatures
valuesForMean=10000

#Simulation temperatures
iniT=30    #Initial temperature
finT=300   #Final temperature  
tempIncrement=10    #Increment of temperature in each simulation
#Pressure (S.I.)
P=1e5   
#Free volumes (Random units)
fbv=0.2 #Free volume of + state
fsv=1 #Free volume of - state
#Volumes (S.I)
vb=4.151444703e-29  #V+
vs=3.321078553e-29  #V-
#Energies   (S.I.)
eb=-1.660577881e-21 #E+
es=0    #E-
#Lattice length (Unit cells)
l=50



resultList=ResultsList() #ResultList object creation
#Calculating response functions for each temperature
for T in range(iniT, finT+tempIncrement, tempIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, FVb, FVs, Eb, Es, P)
    print(T)
     #Calculating response functions with N/meanSteps measurements of the lattice parameters to calculate
    #every response function
    alpha_p=LatticeHandler().thermalExpansivity(T, valuesForMean)
    beta_t=LatticeHandler().isothermalCompressibility(T, valuesForMean)
    c_p=LatticeHandler().heatCapacity(T, valuesForMean)
    k_s=LatticeHandler().isentropicCompressibility(T, valuesForMean)
    #Adding results to list
    resultList.addResult(Result(l, T, P, vb, vs, fbv, fsv, eb, es, alpha_p, beta_t, c_p, k_s))
#Saving results in Results.txt file using the list of inputdatas
resultList.saveResults(iniT, finT, tempIncrement)
#Plotting graphs
resultList.plotAlpha_P()
resultList.plotBeta_T()
resultList.plotC_P()
resultList.plotK_S()