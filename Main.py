#Class to run the simulations

from InputData import *
from LatticeHandler import *
from Lattice import *
from Result import *
from ResultsList import *
import numpy as np


#Creating result list to store the result of the simulation
resultList=ResultsList()
#Number of steps to get record system state
meanSteps=100
#Number of measurements to do the mean
valuesForMean=100
#Simulation temperatures
iniT=200
finT=300
tempIncrement=10
#simulation steps
mcSteps=10000000
Eq=30000000


for T in range(iniT, finT+tempIncrement, tempIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, FVb, FVs, Eb, Es, P)
    inputdata=InputData(50, T, mcSteps, Eq,
                  3.653271338e-29, 3.321078553e-29, 1, 5, -1.660577881e-21, 0, 5e7)
    print(T)
    #Creation of the lattice with the parameters given by the inputData
    lattice=Lattice(inputdata)
    print("Thermalization")
    counter=1
    #Setting equilibrium, steps defined by the inputData
    for i in range(0, inputdata.getEq()):
        if math.fmod(counter/Eq*100,1.0)==0:
            print(str(counter/(Eq)*100)+"%")
        counter+=1
        LatticeHandler().changeVol(lattice, inputdata)

    #Getting measurements of the evolution of volume, intermolecular energy and enthalpy
    measurements=LatticeHandler().getSystemEvolution(lattice, inputdata, meanSteps)
   