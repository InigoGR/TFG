#'''
#Created on 17 feb. 2019

#@author: Inigo Gonzalez
#'''

#Result.py
#Creates an object in which the results and temperature of a simulation
#will be stored

class Result:
    #Class constructor, stores results and their standard deviation. Alpha_P and the other
    #parameters are expected to be arrays [result, stdDeviation]
    def __init__(self, l, temp, press, vb, vs, fbv, fsv, eb, es, alpha_p, beta_t, c_p, k_s):
        self.l=l
        self.temperature=temp
        self.press=press
        self.vb=vb
        self.vs=vs
        self.lambdaVol=fbv/fsv
        self.eb=eb
        self.es=es
        #Assigning result
        self.alpha_p=alpha_p[0]
        #Assigning error of the result
        self.alphaDev=alpha_p[1]
        self.beta_t=beta_t[0]
        self.betaDev=beta_t[1]
        self.c_p=c_p[0]
        self.cDev=c_p[1]
        self.k_s=k_s[0]
        self.kDev=k_s[1]
    
    #Return methods
    def getL(self):
        return self.l

    def getT(self):
        return self.temperature
    
    def getP(self):
        return self.press

    def getVPlus(self):
        return self.vb

    def getVMinus(self):
        return self.vs

    def getLambda(self):
        return self.lambdaVol
    
    def getEPlus(self):
        return self.eb

    def getEMinus(self):
        return self.es

    def getAlpha_P(self):
        return self.alpha_p

    def getAlphaError(self):
        return self.alphaDev
    
    def getBeta_T(self):
        return self.beta_t

    def getBetaError(self):
        return self.betaDev

    def getC_P(self):
        return self.c_p

    def getCError(self):
        return self.cDev

    def getK_S(self):
        return self.k_s

    def getKError(self):
        return self.kDev


