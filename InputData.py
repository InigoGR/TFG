'''
Created on 17 dic. 2018

@author: Inigo Gonzalez
@version: 20/12/2018
'''

#InputData.py
#Creates an object with standard input arguments for the Main class

class InputData:
    #Input object constructor with some default values
    #Parameters: lattice size, temperature, montecarlo steps, equlibrium steps, big cell volume
    #            small cell volume, big free volume, small free volume, big energy, small energy, pressure
    def __init__( self, l=4, t=200, n=1000000, nEquil=1000,
                  vb=2, vs=1, fbv=1.0, fsv=0.2, eb=1, es=2, p=1 ) :
        self.length=l   #Length of the lattice in cells
        self.temperature=t  #Temperature of the simulation
        self.steps=n    #MonteCarlo steps for the measurement phase
        self.equilibrium=nEquil #MonteCarlo steps for the thermalization phase 
        self.bigVolume=vb   #Volume of the V+ state
        self.smallVolume=vs #Volume of the V- state
        self.bigEnergy=eb   #Energy of the V+ V+ interaction
        self.smallEnergy=es #Energy of the V+ V- or V- V- interaction
        self.pressure=p #Pressure of the simulation
        self.freeBVolume=fbv    #Free volume for the V+ state
        self.freeSVolume=fsv    #Free volume for the V- state
        self.freeLambda=fbv/fsv #Lambda parameter

    #Method to change the temperature of the simulation
    def changeT(self, T):
        self.temperature=T
    #Method to change the pressure of the simulation
    def changeP(self,P):
        self.pressure=P

    #Operations to get the data for the simulation
    def getT(self):
        return self.temperature
    
    def getL(self):
        return self.length
    
    def getN(self):
        return self.steps
    
    def getEq(self):
        return self.equilibrium
    
    def getVb(self):
        return self.bigVolume
    
    def getVs(self):
        return self.smallVolume
    
    def getEb(self):
        return self.bigEnergy

    def getEs(self):
        return self.smallEnergy
    
    def getP(self):
        return self.pressure

    def getFbv(self):
        return self.freeBVolume

    def getFsv(self):
        return self.freeSVolume

    def getLambda(self):
        return self.freeLambda