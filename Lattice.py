'''
Created on 15 feb. 2019

@author: Inigo Gonzalez
@version 15/02/19

'''
#Lattice.py
#Class that creates and makes changes to Lattice objects caracterized by their cells
#the total volume and intermolecular energy

from InputData import *
from Constant import *
import random
import math
import numpy as np

class Lattice:
    #Defining attributes of the lattice class
    lattice=0   #Matrix containing the states of the cells (V+=1, V-=0)
    volume=0    #Volume of the lattice
    energy=0    #Energy of the lattice
    neighbor=0  #Matrix which contains for each cell the coordinates of its six neighbors
    probArrayPlus=[0, 0, 0, 0, 0, 0, 0] #Probability of transition form V+ to V- for the seven possible cases
    probArrayMinus=[0, 0, 0, 0, 0, 0, 0]    #Probability of transition form V- to V+ for the seven possible cases
    nVPlus=0    #Number of cells in the V+ state
    nVMinus=0   #Number of ccells in the V- state
    inputData=0 #Object containing the parameters of the simulation
    neighborHistogram=[0, 0, 0, 0, 0, 0, 0] #Array to store the statistics of neighbors state
    #Constructor
    def __init__(self, inputData):
        #Assigning inputData object as attribute
        self.inputData=inputData
        #Creating random lattice
        self.lattice=np.zeros((inputData.getL(),inputData.getL(),inputData.getL()),\
                         dtype=int)
        #Randomizing state of the cells
        for x in range(0, inputData.getL()):
            for y in range(0, inputData.getL()):
                for z in range(0, inputData.getL()):
                    self.lattice[x,y,z]=random.randint(0,1)

        #Calculating total volume
        for column in self.lattice:
            for row in column:
                for cell in row:
                    if cell==1:
                        self.nVPlus+=1
                    else:
                        self.nVMinus+=1
        self.volume=self.nVMinus*inputData.getVs()+self.nVPlus*inputData.getVb()

        #Getting neighbor lattice
        self.neighbor=np.zeros((inputData.getL(),inputData.getL(),inputData.getL(),6,3),\
                         dtype=int)
        #Assigning neighbors to each cell
        for x in range(0, inputData.getL()):
            for y in range(0, inputData.getL()):
                for z in range(0, inputData.getL()):
                    #Creating in each cell a list of the neighbor
                    #with their coordinatess taking into account periodic
                    #boundary conditions

                    #Default neighbor index
                    Xindex1=x+1
                    Xindex2=x-1
                    Yindex1=y+1
                    Yindex2=y-1
                    Zindex1=z+1
                    Zindex2=z-1
                    #Checking for boundary condition
                    if x==0:
                        Xindex2=inputData.getL()-1
                    elif x==inputData.getL()-1:
                        Xindex1=0

                    if y==0:
                        Yindex2=inputData.getL()-1
                    elif y==inputData.getL()-1:
                        Yindex1=0

                    if z==0:
                        Zindex2=inputData.getL()-1
                    elif z==inputData.getL()-1:
                        Zindex1=0
                    #Coordinates of each neighbor
                    self.neighbor[x,y,z]=[[Xindex1, y, z], [Xindex2, y, z],\
                                            [x, Yindex1, z], [x, Yindex2, z],\
                                            [x, y, Zindex1], [x, y, Zindex2]]

        #Calculating intermolecular energy
        energy=0
        #Loop to get the initial intermolecular energy
        for x in range(0, inputData.getL()):
            for y in range(0, inputData.getL()):
                for z in range(0, inputData.getL()):
                    #If the cell has the small volume all the neighbor interactions are the same
                    if self.lattice[x,y,z]==0:
                        energy=energy+6*inputData.getEs()
                    else:
                        bigVol=0
                        #Getting the number of neighbors with big volume
                        for neighborCoord in self.neighbor[x,y,z]:
                            bigVol=bigVol+self.lattice[neighborCoord[0], neighborCoord[1],\
                                 neighborCoord[2]]
                        #Calculating the intermolecular energy of neighbors
                        energy=energy+(bigVol*inputData.getEb()+(6-bigVol)*inputData.getEs())
        self.energy=energy/2 #The interaction between molecules was accounted for twice before

        #Calculating deltaE and deltaV
        deltaE=inputData.getEs()-inputData.getEb()
        deltaV=inputData.getVb()-inputData.getVs()
        #Probabilities for energy change from V+ to V-
        self.probArrayPlus[0]=math.pow(inputData.getLambda(),-1)*math.exp(-((0*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[1]=math.pow(inputData.getLambda(),-1)*math.exp(-((1*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[2]=math.pow(inputData.getLambda(),-1)*math.exp(-((2*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[3]=math.pow(inputData.getLambda(),-1)*math.exp(-((3*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[4]=math.pow(inputData.getLambda(),-1)*math.exp(-((4*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[5]=math.pow(inputData.getLambda(),-1)*math.exp(-((5*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayPlus[6]=math.pow(inputData.getLambda(),-1)*math.exp(-((6*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        #Probabilities for energy change from V- to V+
        self.probArrayMinus[0]=math.pow(inputData.getLambda(),1)*math.exp(((0*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[1]=math.pow(inputData.getLambda(),1)*math.exp(((1*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[2]=math.pow(inputData.getLambda(),1)*math.exp(((2*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[3]=math.pow(inputData.getLambda(),1)*math.exp(((3*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[4]=math.pow(inputData.getLambda(),1)*math.exp(((4*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[5]=math.pow(inputData.getLambda(),1)*math.exp(((5*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
        self.probArrayMinus[6]=math.pow(inputData.getLambda(),1)*math.exp(((6*deltaE-inputData.getP()*deltaV)/Constant().K()/inputData.getT()))
    
    #Method to change the volume of a cell
    def changeCellVolume(self, x, y, z):
        self.lattice[x,y,z]=1-self.lattice[x,y,z]
    #Method to change the energy
    def changeEnergy(self, newEnergy):
        self.energy=newEnergy
    #Method to change the volume
    def changeVolume(self, newVolume):
        self.volume=newVolume
    #Method to register the transition of a minus cell to a plus cell
    def changeMinusToPlus(self):
        self.nVMinus-=1
        self.nVPlus+=1
    #Method to register the transition of a plus cell to a minus cell
    def changePlusToMinus(self):
        self.nVMinus+=1
        self.nVPlus-=1
    #Method to change the inputData attribute, recalculates the change probabilities
    def changeInputData(self, T, P):
        self.inputData.changeT(T)
        self.inputData.changeP(P)

       #Calculating deltaE and deltaV
        deltaE=self.inputData.getEs()-self.inputData.getEb()
        deltaV=self.inputData.getVb()-self.inputData.getVs()
        #Probabilities for energy change from V+ to V-
        self.probArrayPlus[0]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((0*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[1]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((1*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[2]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((2*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[3]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((3*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[4]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((4*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[5]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((5*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayPlus[6]=math.pow(self.inputData.getLambda(),-1)*math.exp(-((6*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        #Probabilities for energy change from V- to V+
        self.probArrayMinus[0]=math.pow(self.inputData.getLambda(),1)*math.exp(((0*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[1]=math.pow(self.inputData.getLambda(),1)*math.exp(((1*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[2]=math.pow(self.inputData.getLambda(),1)*math.exp(((2*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[3]=math.pow(self.inputData.getLambda(),1)*math.exp(((3*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[4]=math.pow(self.inputData.getLambda(),1)*math.exp(((4*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[5]=math.pow(self.inputData.getLambda(),1)*math.exp(((5*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
        self.probArrayMinus[6]=math.pow(self.inputData.getLambda(),1)*math.exp(((6*deltaE-self.inputData.getP()*deltaV)/Constant().K()/self.inputData.getT()))
    #Method to change the histogram array components, add +1 to the n component
    def changeHistogram(self, n):
        self.neighborHistogram[n]+=1
    #Method to reset the histogram values to 0
    def resetHistogram(self):
        self.neighborHistogram=[0, 0, 0, 0, 0, 0, 0]
    #Method to save the histogram values
    def saveHistogram(self):
         #File to write simulation mean energy values
        writeFileName="Neighbor_Histogram"
        fw=open(writeFileName, "a")
        #Header of the file
        fw.write("T="+str(self.inputData.getT())+"\n")
        fw.write("nV+neighbors"+"\t"+"nCells"+"\n")
        for i in range(0,7):
            fw.write(str(i)+"\t"+str(self.neighborHistogram[i])+"\n")
    #Method to get the simulation parameters
    def getInputData(self):
        return self.inputData
    #Method to return the lattice
    def getLattice(self):
        return self.lattice
    #Method to return the neighbor lattice
    def getNeighbors(self):
        return self.neighbor
    #Method to return the volume
    def getVolume(self):
        return self.volume
    #Method to get the energy
    def getEnergy(self):
        return self.energy
    #Method to get the number of V+ cells
    def getNVPlus(self):
        return self.nVPlus
    #Method to get the number of V- cells
    def getNVMinus(self):
        return self.nVMinus
