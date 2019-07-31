#Class to get the Volume as a function of the Pressure of the system

from InputData import *
from LatticeHandler import *
from Lattice import *
import numpy as np
import os

#Temperature
T=250
#Number of steps to get record system state
meanSteps=100
#Number of measurements to do the mean
valuesForMean=100
#Simulation pressures
iniP=int(1e5)
finP=int(5.01e8)
pressureIncrement=int(5e7)
#simulation steps
mcSteps=10000000
initialEq=30000000
Eq=8000000

inputdata=InputData(50, T, mcSteps, Eq,
                  3.653271338e-29, 3.321078553e-29, 1, 5, -1.660577881e-21, 0, iniP)
lattice=Lattice(inputdata)

#Initial equilibrium steps
print("Initial thermalization")
for i in range(0, initialEq):
        if math.fmod(i/initialEq*100,1.0)==0:
            print(str(i/(initialEq)*100)+"%")
        LatticeHandler().changeVol(lattice, lattice.getInputData())


for P in range(iniP, finP+pressureIncrement, pressureIncrement):
    print(P)
    #Changing pressure of the simulation
    lattice.changeInputData(T, P)
    print("Thermalization")
    #Setting new thermalization for new parameters
    for i in range(0, inputdata.getEq()):
        if math.fmod(i/Eq*100,1.0)==0:
            print(str(i/(Eq)*100)+"%")
        LatticeHandler().changeVol(lattice, lattice.getInputData())

    #Getting measurements of the evolution of volume, intermolecular energy and enthalpy
    measurements=LatticeHandler().getSystemEvolution(lattice, lattice.getInputData(), meanSteps)
    #Renaming data files 
    newFileName="Measurements_P"+str(P)
    oldFileName="Measurements_T"+str(T)
    os.rename(oldFileName, newFileName)