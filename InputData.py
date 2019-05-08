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
                  vb=2, vs=1, fbv=1, fsv=0.2, eb=1, es=2, p=1 ) :
        self.length=l
        self.temperature=t
        self.steps=n
        self.equilibrium=nEquil
        self.bigVolume=vb
        self.smallVolume=vs
        self.bigEnergy=eb
        self.smallEnergy=es
        self.pressure=p
        self.freeBVolume=fbv
        self.freeSVolume=fsv
        self.freeLambda=fbv/fsv

    
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