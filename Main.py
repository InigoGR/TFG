#Class to run the simulations

from InputData import *
from LatticeHandler import *
from Lattice import *
from Result import *
from ResultsList import *
import numpy as np
import os


#Creating result list to store the result of the simulation
resultList=ResultsList()
#Number of steps to get record system state
meanSteps=100
#Number of measurements to do the mean
valuesForMean=100
#Simulation temperatures
iniT=200    #Initial temperature
finT=300   #Final temperature  
tempIncrement=10    #Increment of temperature in each simulation
#simulation steps
mcSteps=10000000    #MonteCarlo steps during measurement phase 
initialEq=30000000  #Steps to reach the equilibrium before the simulation
Eq=10000000 #Steps to reach equilibrium after changing the temperature of the system
#Pressure (S.I.)
P=5e7   
#Free volumes (Random units)
fbv=0.2 #Free volume of + state
fsv=1 #Free volume of - state
#Volumes (S.I)
Vb=4.151444703e-29  #V+
Vs=3.321078553e-29  #V-
#Energies   (S.I.)
Eb=-1.660577881e-21 #E+
Es=0    #E-
#Lattice length (Unit cells)
l=50


#Removing probability data file to create a new one
if os.path.isfile("probabilities"):
    os.remove("probabilities")

inputdata=InputData(l, iniT, mcSteps, Eq,
                  Vb, Vs, fbv, fsv, Eb, Es, P)   #Creating Inputdata object containing the initial simulation parameters
lattice=Lattice(inputdata)  #Creating Lattice object using the simulation parameters contained in the Inputdata object

#Initial equilibrium steps
print("Initial thermalization")
for i in range(0, initialEq):
        if math.fmod(i/initialEq*100,1.0)==0: #Checking progress of the equilibrium steps
            print(str(i/(initialEq)*100)+"%")
        LatticeHandler().changeVol(lattice, lattice.getInputData()) #Attempting to change state of one cell

#Removing histogram data file to create a new one
if os.path.isfile("Neighbor_Histogram"):
    os.remove("Neighbor_Histogram")

#Repeating simulation for the specified temperatures
for T in range(iniT, finT+tempIncrement, tempIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, FVb, FVs, Eb, Es, P)
    print(T)
    #Changing temperature of the simulation
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
    lattice.saveHistogram() #Saving neighbor histogram