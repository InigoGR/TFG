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
    lattice=0
    volume=0
    energy=0
    neighbor=0
    probArrayPlus=[0, 0, 0, 0, 0, 0, 0]
    probArrayMinus=[0, 0, 0, 0, 0, 0, 0]
    #Constructor
    def __init__(self, inputData):
        #Creating random lattice
        self.lattice=np.zeros((inputData.getL(),inputData.getL(),inputData.getL()),\
                         dtype=int)
        #Randomizing cells
        for x in range(0, inputData.getL()):
            for y in range(0, inputData.getL()):
                for z in range(0, inputData.getL()):
                    self.lattice[x,y,z]=random.randint(0,1)

        #Calculating total volume
        nVPlus=0
        nVMinus=0
        for column in self.lattice:
            for row in column:
                for cell in row:
                    if cell==1:
                        nVPlus+=1
                    else:
                        nVMinus+=1
        self.volume=nVMinus*inputData.getVs()+nVPlus*inputData.getVb()

        #Getting neighbor lattice
        self.neighbor=np.zeros((inputData.getL(),inputData.getL(),inputData.getL(),6,3),\
                         dtype=int)
        #Assigning neighbors to each cell
        for x in range(0, inputData.getL()):
            for y in range(0, inputData.getL()):
                for z in range(0, inputData.getL()):
                    #Creating in each cell a list of the neighbor
                    #with their coords taking into account periodic
                    #boundary conditions

                    #Default neighbor index
                    Xindex1=x+1
                    Xindex2=x-1
                    Yindex1=y+1
                    Yindex2=y-1
                    Zindex1=z+1
                    Zindex2=z-1
                    #Checking for boundary condiction
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
        self.energy=energy/2
        #Volume difference in 1 cell
        vChange=inputData.getVb()-inputData.getVs()
        #Energy difference in 1 volume change
        eChange=inputData.getEb()-inputData.getEs()
        #Probabilities for energy change from V+ to V-
        self.probArrayPlus[0]=1
        self.probArrayPlus[1]=math.pow(inputData.getLambda(),-1)*math.exp(-((eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayPlus[2]=math.pow(inputData.getLambda(),-1)*math.exp(-((2*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayPlus[3]=math.pow(inputData.getLambda(),-1)*math.exp(-((3*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayPlus[4]=math.pow(inputData.getLambda(),-1)*math.exp(-((4*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayPlus[5]=math.pow(inputData.getLambda(),-1)*math.exp(-((5*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayPlus[6]=math.pow(inputData.getLambda(),-1)*math.exp(-((6*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        #Probabilities for energy change from V- to V+
        self.probArrayMinus[0]=1
        self.probArrayMinus[1]=math.pow(inputData.getLambda(),1)*math.exp(((eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayMinus[2]=math.pow(inputData.getLambda(),1)*math.exp(((2*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayMinus[3]=math.pow(inputData.getLambda(),1)*math.exp(((3*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayMinus[4]=math.pow(inputData.getLambda(),1)*math.exp(((4*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayMinus[5]=math.pow(inputData.getLambda(),1)*math.exp(((5*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
        self.probArrayMinus[6]=math.pow(inputData.getLambda(),1)*math.exp(((6*eChange+inputData.getP()*vChange)/Constant().K()/inputData.getT()))
    #Method to change the volume of a cell
    def changeCellVolume(self, x, y, z):
        self.lattice[x,y,z]=1-self.lattice[x,y,z]
    #Method to change the energy
    def changeEnergy(self, newEnergy):
        self.energy=newEnergy
    #Method to change the volume
    def changeVolume(self, newVolume):
        self.volume=newVolume

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