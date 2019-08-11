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
iniP=int(1e5) #Initial pressure
finP=int(5.01e8)    #Final pressure
pressureIncrement=int(5e7)  #Increment of pressure in each simulation
#simulation steps
mcSteps=10000000        #MonteCarlo steps during measurement phase 
initialEq=30000000  #Steps to reach the equilibrium before the simulation
Eq=8000000  #Steps to reach equilibrium after changing the temperature of the system

inputdata=InputData(50, T, mcSteps, Eq,
    4.151444703e-29, 3.321078553e-29, 1, 5, -1.660577881e-21, 0, iniP)    #Creating Inputdata object containing the initial simulation parameters
lattice=Lattice(inputdata)  #Creating Lattice object using the simulation parameters contained in the Inputdata object

#Initial equilibrium steps
print("Initial thermalization")
for i in range(0, initialEq):
        if math.fmod(i/initialEq*100,1.0)==0:   #Checking progress of the equilibrium steps
            print(str(i/(initialEq)*100)+"%")
        LatticeHandler().changeVol(lattice, lattice.getInputData()) #Attempting to change state of one cell

#Removing histogram data file to create a new one
if os.path.isfile("Neighbor_Histogram"):
    os.remove("Neighbor_Histogram")

#Repeating simulation for the specified temperatures
for P in range(iniP, finP+pressureIncrement, pressureIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, FVb, FVs, Eb, Es, P)
    print(P)
    #Changing pressure of the simulation
    lattice.changeInputData(T, P)   #Changing simulation parameters
    print("Thermalization")
    #Setting new thermalization for new parameters
    for i in range(0, inputdata.getEq()):
        if math.fmod(i/Eq*100,1.0)==0:  #Checking progress of the equilibrium steps
            print(str(i/(Eq)*100)+"%")
        LatticeHandler().changeVol(lattice, lattice.getInputData()) #Attempting to change state of one cell
    lattice.resetHistogram()
    #Getting measurements of the evolution of volume, intermolecular energy and enthalpy
    measurements=LatticeHandler().getSystemEvolution(lattice, lattice.getInputData(), meanSteps)
    #Renaming data files 
    newFileName="Measurements_P"+str(P)
    oldFileName="Measurements_T"+str(T)
    os.rename(oldFileName, newFileName)
    lattice.saveHistogram() #Saving neighbor histogram