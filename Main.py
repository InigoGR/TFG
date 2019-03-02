'''
Created on 17 dic. 2018

@author: Inigo Gonzalez
@version 27/02/19

'''
#Class to run the simulation with the specified input parameters

from InputData import *
from LatticeHandler import *
from Lattice import *
from Result import *
from ResultsList import *
import numpy as np


#Creating result list to store the result of the simulation
resultList=ResultsList()
tempIncrement=50
#List of the used inputData objects
inputdata=[]
#Iteration counter to choose the inputData from the list
j=0

for T in range(100, 300+tempIncrement, tempIncrement):
    #Creation of the inputData with the desired parameters(L, T, MCs, Therm, Vb, Vs, Eb, Es, P)
    inputdata.append(InputData(4, T, 1000000, 1000,
                  2, 1, 2, 1, 1 ))
    print(inputdata[j].getT())
    #Creation of the lattice with the parameters given by the inputData
    lattice=Lattice(inputdata[j])
    #Setting equilibrium, steps defined by the inputData
    for i in range(0, inputdata[j].getEq()):
        LatticeHandler().changeVol(lattice, inputdata[j])
    #Getting measurements of the evolution of volume, intermolecular energy and enthalpy
    measurements=LatticeHandler().getSystemEvolution(lattice, inputdata[j], 1000)
    #Calculating response functions with 100 measurements of the lattice parameters to calculate
    #every response function
    alpha_p=LatticeHandler().thermalExpansivity(measurements[0], measurements[2], inputdata[j], 100)
    beta_t=LatticeHandler().isothermalCompressibility(measurements[0], inputdata[j], 100)
    c_p=LatticeHandler().heatCapacity(measurements[0], measurements[1], measurements[2], inputdata[j], 100)
    #Adding results to list
    resultList.addResult(Result(inputdata[j].getT(), alpha_p, beta_t, c_p))
    #Next iteration
    j+=1
#Saving results in Results.txt file using the list of inputdatas
resultList.saveResults(inputdata)
#Plotting graphs
resultList.plotAlpha_P()
resultList.plotBeta_T()
resultList.plotC_P()

