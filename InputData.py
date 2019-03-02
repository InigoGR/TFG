'''
Created on 17 dic. 2018

@author: Inigo Gonzalez
@version: 20/12/2018
'''

#InputData.py
#Creates an object with standard input arguments for the Main class

class InputData:
    #Input object constructor with some default values
    #Parameters: lattice size, temperature, montecarlo steps, equlibrium steps, big volume
    #            small volume, big energy, small energy, pressure
    def __init__( self, l=4, t=200, n=1000000, nEquil=1000,
                  vb=2, vs=1, eb=2, es=1, p=1 ) :
        self.length=l
        self.temperature=t
        self.steps=n
        self.equilibrium=nEquil
        self.bigVolume=vb
        self.smallVolume=vs
        self.bigEnergy=eb
        self.smallEnergy=es
        self.pressure=p
    
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