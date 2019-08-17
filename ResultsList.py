#'''
#Created on 17 dic. 2018

#@author: Inigo Gonzalez
#'''

#Results.py
#Creates an object in which the results of the simulations
#will be stored

from Result import *
import matplotlib.pyplot as plt
import math
import numpy as np
from Constant import *

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
    def saveResults(self, iniT, finT, tempIncrement):
        #Creating file to store results
        fw=open("Results.txt", "w")
        line_new ="T"+"\t"+"Alpha_P"+"\t"+"Beta_T"+"\t"+"C_P"+"\t"+"K_S"
        fw.write(line_new+"\n")
        #Creating temperature array
        temperatures=[]
        for T in range(iniT, finT+tempIncrement, tempIncrement):
            temperatures.append(T)
        #Writing results in a file
        for i in range(0,len(self.resultList)):
            line_new =str(temperatures[i])+"\t"+str(self.resultList[i].getAlpha_P())+"+-"+\
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
        
        #Calculating theoretical values
        #Array of volumes
        v=np.linspace(2.14e-5, 2.40e-5, num=10000)
        
        #Asigning necessary values for the theoretical calculation
        vs=self.resultList[0].getVMinus()*6.022e23   
        deltaV=self.resultList[0].getVPlus()*6.022e23-vs
        deltaE=(self.resultList[0].getEMinus()-self.resultList[0].getEPlus())*6.022e23
        p=self.resultList[0].getP()
        lambdaVol=self.resultList[0].getLambda()   

        #Array of temperatures
        T=[]
        for volume in v:
            T.append((p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs)))

        #Turning to cubic meters
        v=v/6.022e23
        #Volume of the whole lattice
        v=np.multiply(v,50*50*50)

        alpha=[]
        plotT=[]

        k=1 #Iteration counter
        while k<(len(T)-1):
            derivative=(v[k+1]-v[k-1])/((T[k+1]-T[k-1]))
            alpha.append(derivative/v[k])
            plotT.append(T[k])
            k+=1

        plt.errorbar(x,y,error, marker='o', label='Sim')
        plt.plot(plotT, alpha, linestyle='dashed', label='Theory')
        plt.legend()
        plt.xlabel('T/K')
        plt.title(r'$\alpha_P$', fontsize=14)
        plt.show()

    #Method to plot the isothermal compressibility as a function of temperature
    def plotBeta_T(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getBeta_T())
            error.append(result.getBetaError())
        
        #Calculating theoretical values
        #Asigning necessary values for the theoretical calculation
        vs=self.resultList[0].getVMinus()*6.022e23   
        deltaV=self.resultList[0].getVPlus()*6.022e23-vs
        deltaE=(self.resultList[0].getEMinus()-self.resultList[0].getEPlus())*6.022e23
        p=self.resultList[0].getP()
        lambdaVol=self.resultList[0].getLambda()
        c=6 #Coordination
        #Calculating 1/V*(dv/dT) at cte pressure for all temperatures
        Bt=[]
        T=[]
        p1=[]
        p2=[]
        p3=[]
        k=0
        #Array of volumes
        v=np.linspace(2.14e-5, 2.40e-5, num=10000)
        #Calculating temperature
        for volume in v:
            tempT=(p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs))
            T.append(tempT)
            p1.append(tempT*8.314/deltaV*math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs)))
            p2.append(c*deltaE/deltaV*(volume-vs)/deltaV)
            p3.append(p1[k]+p2[k])
            k+=1
        
        #Turning to cubic meters
        v=v/6.022e23
        #Volume of the whole lattice
        v=np.multiply(v,50*50*50)      
        
        plotT=[]
        #Calculation of the response function
        k=1 #Iteration counter
        while k<(len(T)-1):
            #Value of the derivative
            derivative=(v[k+1]-v[k-1])/((p3[k+1]-p3[k-1]))
            Bt.append(-derivative/v[k])
            plotT.append(T[k])
            k+=1

        #Plotting results
        plt.errorbar(x,y,error, marker='o', label='Sim')
        plt.plot(plotT, Bt, linewidth=2.0, linestyle='dashed', label='Theory')
        plt.legend()
        plt.title(r'$\beta_T$', fontsize=14)
        plt.show()
        
    
    #Method to plot the residual heat capacity as a function of temperature
    def plotC_P(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getC_P())
            error.append(result.getCError())

        #Calculating theoretical values

        #Asigning necessary values for the theoretical calculation
        vs=self.resultList[0].getVMinus()*6.022e23   
        deltaV=self.resultList[0].getVPlus()*6.022e23-vs
        deltaE=(self.resultList[0].getEMinus()-self.resultList[0].getEPlus())*6.022e23
        p=self.resultList[0].getP()
        lambdaVol=self.resultList[0].getLambda()
        l=self.resultList[0].getL()
        c=6
        #Array of volumes
        v=np.linspace(2.14e-5, 2.40e-5, num=10000)
            

        #Array of temperatures
        T=[]
        p3=[]
        k=0
        for volume in v:
            tempT=(p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs))
            T.append(tempT)
            p1=tempT*8.314/deltaV*math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs))
            p2=c*deltaE/deltaV*(volume-vs)/deltaV
            p3.append(p1+p2)
            k+=1
        

        #Performing unit change
        vs=vs/6.022e23
        v=v/6.022e23   
        deltaV=deltaV/6.022e23  
        deltaE=deltaE/6.022e23    

        #Energy per unit cell
        energyUnitCell=np.multiply(-3*deltaE, np.power((v-vs)/deltaV,2))
        #Energy of the whole lattice
        E=np.multiply(energyUnitCell, math.pow(l,3))

        Cp=[]
        plotT=[]
        #Calculation of the response function
        k=1 #Iteration counter
        while k<(len(T)-1):
            #Calculating alpha
            derivativeAlpha=(v[k+1]-v[k-1])/((T[k+1]-T[k-1]))
            alpha=derivativeAlpha/v[k]
            #Calculating beta
            derivativeBeta=(v[k+1]-v[k-1])/((p3[k+1]-p3[k-1]))
            Bt=-derivativeBeta/v[k]
            #Value of the derivative
            derivative=(E[k+1]-E[k-1])/((T[k+1]-T[k-1]))
            Cp.append(derivative/T[k]+math.pow(l,3)*Constant().K()+T[k]*v[k]*math.pow(alpha,2)/Bt)
            plotT.append(T[k])
            k+=1

        plt.errorbar(x,y,error,marker='o', label='Sim')
        plt.plot(plotT, Cp, linewidth=2.0, linestyle='dashed', label='Theory')
        plt.legend()
        plt.title("$C_p$", fontsize=14)
        plt.show()

    #Method to plot the residual isentropic compressibility as a function of temperature
    def plotK_S(self):
        x=[]
        y=[]
        error=[]
        for result in self.resultList:
            x.append(result.getT())
            y.append(result.getK_S())
            error.append(result.getKError())
        
        #Calculating theoretical values
        #Asigning necessary values for the theoretical calculation
        vs=self.resultList[0].getVMinus()*6.022e23   
        deltaV=self.resultList[0].getVPlus()*6.022e23-vs
        deltaE=(self.resultList[0].getEMinus()-self.resultList[0].getEPlus())*6.022e23
        p=self.resultList[0].getP()
        lambdaVol=self.resultList[0].getLambda()
        l=self.resultList[0].getL()
        c=6 #Coordination
        
        T=[]
        p1=[]
        p2=[]
        p3=[]
        k=0
        #Array of volumes
        v=np.linspace(2.14e-5, 2.40e-5, num=10000)
        #Calculating temperature
        for volume in v:
            tempT=(p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs))
            T.append(tempT)
            p1.append(tempT*8.314/deltaV*math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs)))
            p2.append(c*deltaE/deltaV*(volume-vs)/deltaV)
            p3.append(p1[k]+p2[k])
            k+=1
        
        #Performing unit change
        vs=vs/6.022e23
        v=v/6.022e23   
        deltaV=deltaV/6.022e23  
        deltaE=deltaE/6.022e23    

        #Energy per unit cell
        energyUnitCell=np.multiply(-3*deltaE, np.power((v-vs)/deltaV,2))
        #Energy of the whole lattice
        E=np.multiply(energyUnitCell, math.pow(l,3))

        #Volume of the whole lattice
        v=np.multiply(v,50*50*50)      
        
        Ks=[]
        plotT=[]
        #Calculation of the response function
        k=1 #Iteration counter
        while k<(len(T)-1):
            #Calculating alpha
            derivativeAlpha=(v[k+1]-v[k-1])/((T[k+1]-T[k-1]))
            Alpha=derivativeAlpha/v[k]
            #Calculating beta
            derivativeBeta=(v[k+1]-v[k-1])/((p3[k+1]-p3[k-1]))
            Bt=-derivativeBeta/v[k]
            #Value of the derivative
            derivative=(E[k+1]-E[k-1])/((T[k+1]-T[k-1]))
            Cp=derivative/T[k]+math.pow(l,3)*Constant().K()+T[k]*v[k]*math.pow(Alpha,2)/Bt
            Ks.append(Bt-T[k]*v[k]/math.pow(l,3)*math.pow(Alpha,2)/Cp)
            plotT.append(T[k])
            k+=1

        
        plt.errorbar(x,y,error,marker='o', label='Sim')
        plt.plot(plotT, Ks, linewidth=2.0, linestyle='dashed', label='Theory')
        plt.title("$K_s$", fontsize=14)
        plt.legend()
        plt.show()
        
    
        
            



