#'''
#Created on 17 dic. 2018

#@author: Inigo Gonzalez
#'''

#Results.py
#Creates an object in which the results of the simulations
#will be stored

from Result import *
from matplotlib.pyplot import show, errorbar

class ResultsList:
    #Class constructor, creates list to store results
    def __init__(self):
        self.resultList=[]
    
    #Class to get a result from the list
    def getResult(self, index):
        return self.resultList[index]

    #Class to add results to the list
    def addResult(self, result):
        self.resultList.append(result)

    #Class to save results in a file, needs the list of all the inputDatas used
    def saveResults(self, inputDataList):
        #Creating file to store results
        fw=open("Results.txt", "w")
        line_new ="T"+"\t"+"Alpha_P"+"\t"+"Beta_T"+"\t"+"C_P"+"\t"+"K_S"
        fw.write(line_new+"\n")
        #Writing results in a file
        for i in range(0,len(self.resultList)):
            line_new =str(inputDataList[i].getT())+"\t"+str(self.resultList[i].getAlpha_P())+"+-"+\
            str(self.resultList[i].getAlphaError())+"\t"+str(self.resultList[i].getBeta_T())+"+-"+\
            str(self.resultList[i].getBetaError())+"\t"+str(self.resultList[i].getC_P())+"+-"+\
            str(self.resultList[i].getCError())+"\t"+str(self.resultList[i].getK_S())+"+-"+\
            str(self.resultList[i].getKError())
            fw.write(line_new+"\n")
        fw.close()

    #Plotting methods
    #Method to plot the thermal expansivity as a function of temperature
    def plotAlpha_P(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getAlpha_P())
            error.append(result.getAlphaError())
        errorbar(x,y,error,marker='o')
        show()

    #Method to plot the isothermal compressibility as a function of temperature
    def plotBeta_T(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getBeta_T())
            error.append(result.getBetaError())
        errorbar(x,y,error,marker='o')
        show()
    
    #Method to plot the residual heat capacity as a function of temperature
    def plotC_P(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getC_P())
            error.append(result.getCError())
        errorbar(x,y,error,marker='o')
        show()

    #Method to plot the residual isentropic compressibility as a function of temperature
    def plotK_S(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getK_S())
            error.append(result.getKError())
        errorbar(x,y,error,marker='o')
        show()
        
    


