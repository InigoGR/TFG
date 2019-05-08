#Class to test code

from InputData import *
from LatticeHandler import *
from Lattice import *
from Result import *
from ResultsList import *
import numpy as np


#Creating result list to store the result of the simulation
resultList=ResultsList()
#List of the used inputData objects
inputdata=[]
#Iteration counter to choose the inputData from the list
j=0
#Number of steps to get record system state
meanSteps=100
#Number of measurements to do the mean
valuesForMean=100
#Simulation temperatures
iniT=200
finT=300
tempIncrement=10

for T in range(iniT, finT+tempIncrement, tempIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, FVb, FVs Eb, Es, P)
    inputdata=InputData(20, T, 1000000, 1000,
                  4.151444703e-29, 3.321155762e-29, 1, 5, 0, 1.660577881e-21, 1.2e8)
    print(T)
    #Creation of the lattice with the parameters given by the inputData
    lattice=Lattice(inputdata)
    #Setting equilibrium, steps defined by the inputData
    for i in range(0, inputdata.getEq()):
        LatticeHandler().changeVol(lattice, inputdata)
    #Getting measurements of the evolution of volume, intermolecular energy and enthalpy
    measurements=LatticeHandler().getSystemEvolution(lattice, inputdata, meanSteps)
    #Calculating response functions with N/meanSteps measurements of the lattice parameters to calculate
    #every response function
    alpha_p=LatticeHandler().thermalExpansivity(T, valuesForMean)
    beta_t=LatticeHandler().isothermalCompressibility(T, valuesForMean)
    c_p=LatticeHandler().heatCapacity(T, valuesForMean)
    k_s=LatticeHandler().isentropicCompressibility(T, valuesForMean)
    #Adding results to list
    resultList.addResult(Result(T, alpha_p, beta_t, c_p, k_s))
    #Next iteration
    j+=1
#Saving results in Results.txt file using the list of inputdatas
resultList.saveResults(iniT, finT, tempIncrement)
#Plotting graphs
resultList.plotAlpha_P()
resultList.plotBeta_T()
resultList.plotC_P()
resultList.plotK_S()