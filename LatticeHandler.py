'''
Created on 20 dic. 2018

@author: Iñigo
@version: 20/12/18

'''

#LatticeCounter.py
#It allows to perform various operations with the lattice, getting
#the total volume, the number of cells with V+ or V- 

from InputData import *
from Lattice import *
from Constant import *
from Numericalmethods import *
import numpy as np
import random
import math
import statistics
import matplotlib.pyplot as plt


class LatticeHandler:
    #Method that tries a random cell volume change and checks if the system accepts, making
    #the changes
    def changeVol(self, lattice, inputData):
        x=random.randint(0,inputData.getL()-1)
        y=random.randint(0,inputData.getL()-1)
        z=random.randint(0,inputData.getL()-1)
        volumeChange=0
        tempLattice=lattice.getLattice()
        neighborLattice=lattice.getNeighbors()
        volume=lattice.getVolume()
        energy=lattice.getEnergy()
        #Counting number of neighbors with V+
        n=0
        for neighborCoord in neighborLattice[x,y,z]:
            n=n+tempLattice[neighborCoord[0], neighborCoord[1], neighborCoord[2]]
        #Getting intermolecular energy change
        #Case initialV=V-
        if tempLattice[x,y,z]==0:
            if random.uniform(0,1)<lattice.probArrayMinus[n]:
            
                intEnergyVar=n*(inputData.getEb()-inputData.getEs())
                volumeChange=inputData.getVb()-inputData.getVs()
                #Changing volume
                lattice.changeVolume(volume+volumeChange)
                #Changing cell volume
                lattice.changeCellVolume(x,y,z)
                #Changing energy
                lattice.changeEnergy(energy+intEnergyVar)

        #Case initialV=V+
        else:
            intEnergyVar=n*(inputData.getEs()-inputData.getEb())
            volumeChange=inputData.getVs()-inputData.getVb()      
        #Checking if the change is accepted
            if random.uniform(0,1)<lattice.probArrayPlus[n]:
                #Changing volume
                lattice.changeVolume(lattice.getVolume()+volumeChange)
                #Changing cell volume
                lattice.changeCellVolume(x,y,z)
                #Changing energy
                lattice.changeEnergy(energy+intEnergyVar)
    
    #Method to get the volume, enthalpy and intermolecular energy in intervals separated
    #by the indicated amount of steps with a total of MonteCarlo steps indicated by the inputData.
    #Returns three arrays as shown
    #[[VolumeValues],
    #[EnergyValues],
    #[EnthalpyValues]]
    #Writes in a file the every meanSteps
    def getSystemEvolution(self, lattice, inputData, meanSteps):
        #Opening file and writting header
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "w")
        #Writting lattice size
        line_new="Lattice length="+str(inputData.getL())
        fw.write(line_new+"\n")
        #writting temperature
        line_new ="Temperature="+str(inputData.getT())
        fw.write(line_new+"\n")
        #Writting number of steps
        line_new="Monte Carlo steps="+str(inputData.getN())
        fw.write(line_new+"\n")
        #Writting equilibrium steps
        line_new="Equilibrium steps="+str(inputData.getEq())
        fw.write(line_new+"\n")
        #Writting big cell volume
        line_new="Big cell volume="+str(inputData.getVb())
        fw.write(line_new+"\n")
        #Writting small cell volume
        line_new="Small cell volume="+str(inputData.getVs())
        fw.write(line_new+"\n")
        #Writting big free volume
        line_new="Big free volume="+str(inputData.getFbv())
        fw.write(line_new+"\n")
        #Writting small free volume
        line_new="Small free volume="+str(inputData.getFsv())
        fw.write(line_new+"\n")
        #Writting energy of big cell interaction
        line_new="+ Interaction energy="+str(inputData.getEb())
        fw.write(line_new+"\n")
        #Writting energy of small cell interaction
        line_new="- Interaction energy="+str(inputData.getEs())
        fw.write(line_new+"\n")
        #Writting pressure
        line_new="Pressure="+str(inputData.getP())
        fw.write(line_new+"\n")
        #Writting mean steps
        line_new="Steps used for measurement="+str(meanSteps)
        fw.write(line_new+"\n")
        #Data columns
        line_new="Volume"+"\t"+"Energy"+"\t"+"Enthalpy"
        fw.write(line_new+"\n")
        #Data arrays
        volumeValues=[]
        energyValues=[]
        enthalpyValues=[]
        #The lattice state is stored every meanSteps steps
        counter=1
        for i in range(0, int(inputData.getN()/meanSteps)):
            vol=0
            energy=0
            enthalpy=0
            #Changing cell volumes meanSteps times
            for j in range(0, meanSteps):
                LatticeHandler().changeVol(lattice, inputData)
            #Getting state of the lattice
            vol=lattice.getVolume()
            energy=lattice.getEnergy()
            enthalpy=energy+inputData.getP()*vol
            #Adding measurements to the data arrays
            volumeValues.append(vol)
            energyValues.append(energy)
            enthalpyValues.append(enthalpy)
            #Saving parameters in the measurements file
            line_new=str(vol)+"\t"+str(energy)+"\t"+str(enthalpy)
            fw.write(line_new+"\n")
            if counter/inputData.getN()/meanSteps%0.01==0:
                print(counter/int(inputData.getN()/meanSteps))
            counter+=1
        fw.close()
        return [volumeValues, energyValues, enthalpyValues]
    
    #Method to calculate the mean thermal expansivity for a given temperature from the data files
    #generated in a previous simulation and the number of values to be used when computing the mean values.
    #Returns an array with the mean value and the standard deviation and writes the values 
    #(alpha values) used for the mean in the measurements file.
    def thermalExpansivity(self, T, valuesForMean):
        #Data sets
        volumeSets=[]
        enthalpySets=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(T)
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Reading inputData parameters
        line=lines[2].split("=")
        n=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])

        #Closing file
        fr.close()
        #Number of data sets
        nDataSets=int(n/meanSteps/valuesForMean)
        #Reading all the measurements before the response functions
        nLine=13
        for j in range(0, nDataSets):
            #Filling volume data set
            volumeSet=[]
            enthalpySet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                enthalpySet.append(float(values[2]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            volumeSets.append(volumeSet)
            enthalpySets.append(enthalpySet)
        
        #Array to store the calculated values
        thermExpArray=[]
        product=[]
        #Number of data sets, necessary to get mean values
        nDataSets=len(volumeSets)
        #Product of measurements
        for i in range(0,nDataSets):
            product.append(np.multiply(volumeSets[i],enthalpySets[i]))
        #Calculation
        for i in range(0,nDataSets):
            #Getting the mean of the data set of this iteration made of valuesForMean values
            meanVolume=sum(volumeSets[i])/valuesForMean
            meanEnthalpy=sum(enthalpySets[i])/valuesForMean
            meanProduct=sum(product[i])/valuesForMean
            thermExpArray.append(1/meanVolume/Constant().K()/math.pow(T,2)*\
                (meanProduct-meanEnthalpy*meanVolume))
        thermExp=sum(thermExpArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((thermExpArray[i]-thermExp),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(T)
        fw=open(fileName, "a")
        line_new="Values for mean="+str(valuesForMean)
        fw.write(line_new+"\n")
        line_new ="Alpha_p"+"\n"+str(thermExpArray)
        fw.write(line_new+"\n")
        return [thermExp, stdDeviation]
        

    #Method to calculate the mean isothermal compressibility of a given temperature from data files generated
    #in a previous simulation and the amont of values to be used when calculating the mean values.
    #Returns an array with the mean value and the standard deviation and 
    #writes the values used for the mean in the measurements file. 
    def isothermalCompressibility(self, T, valuesForMean):
        #Data sets
        volumeSets=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(T)
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Reading inputData parameters
        line=lines[2].split("=")
        n=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])

        #Closing file
        fr.close()
        #Number of data sets
        nDataSets=int(n/meanSteps/valuesForMean)
        #Reading all the measurements before the response functions
        nLine=13
        for j in range(0, nDataSets):
            #Filling volume data set
            volumeSet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            volumeSets.append(volumeSet)
        #Squared volume
        volumeSq=[]
        #Number of data sets necessary to get mean values
        nDataSets=len(volumeSets)
        #Filling squared volume array
        for i in range(0,nDataSets):
             volumeSq.append(np.multiply(volumeSets[i],volumeSets[i]))
        #Array to store the calculated values
        isothCompressArray=[]
        #Calculation
        for i in range(0,nDataSets):
            #Getting the mean of the data set of this iteration made of 10 values
            meanVolume=sum(volumeSets[i])/valuesForMean
            meanVolumeSq=sum(volumeSq[i])/valuesForMean
            isothCompressArray.append(1/meanVolume/Constant().K()/T*(meanVolumeSq-math.pow(meanVolume,2)))
        isothCompress=sum(isothCompressArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((isothCompressArray[i]-isothCompress),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(T)
        fw=open(fileName, "a")
        line_new ="Beta_t"+"\n"+str(isothCompressArray)
        fw.write(line_new+"\n")
        return [isothCompress, stdDeviation]
        

    #Method to calculate the residual heat capacity for a given temperature from the data files generated
    #in a previous similation and the amount of values that have to be used for 
    #the mean. Returns an array with the mean value and the standard deviation.
    def heatCapacity(self, T, valuesForMean):
        #Data sets
        volumeSets=[]
        energySets=[]
        enthalpySets=[]
        #Data sets
        volumeSets=[]
        enthalpySets=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(T)
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Reading inputData parameters
        line=lines[0].split("=")
        l=float(line[1])
        line=lines[2].split("=")
        n=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])

        #Closing file
        fr.close()
        #Number of data sets
        nDataSets=int(n/meanSteps/valuesForMean)
        #Reading all the measurements before the response functions
        nLine=13
        for j in range(0, nDataSets):
            #Filling volume data set
            volumeSet=[]
            energySet=[]
            enthalpySet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                energySet.append(float(values[1]))
                enthalpySet.append(float(values[2]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            volumeSets.append(volumeSet)
            energySets.append(energySet)
            enthalpySets.append(enthalpySet)
        #Product arrays
        productE=[]
        productV=[]
        #Number of data sets, necessary to get the mean values
        nDataSets=len(volumeSets)
        #Product of measurements
        for i in range(0,nDataSets):
            productE.append(np.multiply(volumeSets[i],enthalpySets[i]))
            productV.append(np.multiply(energySets[i],enthalpySets[i]))
        #Array to store calculated values
        heatCapArray=[]
        #Calculation of the heat capacity
        for i in range(0,nDataSets):
            #Getting mean values for the iteration
            meanVolume=sum(volumeSets[i])/valuesForMean
            meanEnergy=sum(energySets[i])/valuesForMean
            meanEnthalpy=sum(enthalpySets[i])/valuesForMean
            meanProductE=sum(productE[i])/valuesForMean
            meanProductV=sum(productV[i])/valuesForMean

            heatCapArray.append(1/Constant().K()/math.pow(T,2)*(meanProductE-\
            meanEnergy*meanEnthalpy)+1/Constant().K()/math.pow(T,2)*(meanProductV-\
            meanVolume*meanEnthalpy)-math.pow(l,3)*Constant().K())
        heatCap=sum(heatCapArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((heatCapArray[i]-heatCap),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(T)
        fw=open(fileName, "a")
        line_new ="C_p"+"\n"+str(heatCapArray)
        fw.write(line_new+"\n")
        return [heatCap, stdDeviation]


    #Method to calculate the isentropic compressibility from the data files stored of the system
    #evolution. It needs the response coefficients calculated previously and the amount of values used
    #for the mean. 
    def isentropicCompressibility(self, T, valuesForMean):
        #Creating data arrays
        volumes=[]
        beta_t=[]
        alpha_p=[]
        c_p=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(T)
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Reading steps before measurement
        line=lines[0].split("=")
        l=float(line[1])
        line=lines[2].split("=")
        n=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])
        #Line number, skipping input data lines
        nLine=13
        #Reading all the volume measurements before the response functions
        for j in range(0, int(n/meanSteps/valuesForMean)):
            #Filling volume data set
            volumeSet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            volumes.append(volumeSet)
        #Skipping "values for mean"
        nLine+=1
        #Skipping "Alpha_P" header line, not used
        nLine+=1
        #Reading alpha_p values, removing brackets from both ends of the line
        lines[nLine]=lines[nLine].replace("[","")
        lines[nLine]=lines[nLine].replace("]","")
        #Spliting line into separate values
        values=lines[nLine].split(",")
        #Introducing values into data array of alpha_p
        for i in range(0, len(values)):
            alpha_p.append(float(values[i]))
        nLine+=1

        #Skipping "Beta_T" header line, not used
        nLine+=1
        #Reading beta_t values, removing brackets from both ends of the line
        lines[nLine]=lines[nLine].replace("[","")
        lines[nLine]=lines[nLine].replace("]","")
        #Spliting line into separate values
        values=lines[nLine].split(",")
        #Introducing values into data array of alpha_p
        for i in range(0, len(values)):
            beta_t.append(float(values[i]))
        nLine+=1

        #Skippimg "C_P" header line, not used
        nLine+=1
        #Reading c_p values, removing brackets from both ends of the line
        lines[nLine]=lines[nLine].replace("[","")
        lines[nLine]=lines[nLine].replace("]","")
        #Spliting line into separate values
        values=lines[nLine].split(",")
        #Introducing values into data array of alpha_p
        for i in range(0, len(values)):
            c_p.append(float(values[i]))
        #Closing file
        fr.close()

        #Calculating mean values
        meanVolumes=[]
        #Introducing mean volume of each set into array
        for volumeSet in volumes:
            meanVolumes.append(sum(volumeSet)/len(volumeSet))
        
        #Array to store values of Ks
        isentCompressArray=[]
        #Calculating values of the isentropic compressibility
        for i in range(0, len(meanVolumes)):
            isentCompress=beta_t[i]-T*meanVolumes[i]/math.pow(l,3)*\
                math.pow(alpha_p[i],2)/c_p[i]
            isentCompressArray.append(isentCompress)
        isentropicCompressibility=sum(isentCompressArray)/len(isentCompressArray)
        #Calculating deviation
        residue=0
        for i in range(0,len(isentCompressArray)):
            residue=residue+math.pow((isentCompressArray[i]-isentropicCompressibility),2)
        stdDeviation=math.sqrt(residue/(len(isentCompressArray)-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(T)
        fw=open(fileName, "a")
        line_new ="K_s"+"\n"+str(isentCompressArray)
        fw.write(line_new+"\n")
        return [isentropicCompressibility, stdDeviation]

    

    #Class that shows the evolution of the volume for a given temperature using the corresponding
    #data file 'Measurements_Tx'. It needs the number of values for the mean
    def volEvo(self, T, valuesForMean):
        #Creating data arrays
        volumes=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(T)
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Reading inputData parameters
        line=lines[0].split("=")
        l=float(line[1])
        line=lines[2].split("=")
        n=float(line[1])
        line=lines[4].split("=")
        vb=float(line[1])
        line=lines[5].split("=")
        vs=float(line[1])
        line=lines[6].split("=")
        fvb=float(line[1])
        line=lines[7].split("=")
        fvs=float(line[1])
        line=lines[8].split("=")
        eb=float(line[1])
        line=lines[9].split("=")
        es=float(line[1])
        line=lines[10].split("=")
        p=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])
        lambdaVol=fvb/fvs
        deltaV=vb-vs
        deltaE=es-eb
        #Number of data sets
        nDataSets=int(n/meanSteps/valuesForMean)
        #Reading all the measurements before the response functions
        nLine=13
        for j in range(0, nDataSets):
            #Filling volume data set
            volumeSet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            volumes.append(volumeSet)
        #Plotting mean values of the volume sets
        x=[]
        y=[]
        error=[]

        #Calculating deviation to get error of every volume
        #k counter
        k=0
        for dataSet in volumes:
            #Initial residue 
            residue=0
            #Calculating mean
            meanVolume=statistics.mean(dataSet)
            #Adding mean to data array
            y.append(meanVolume)
            x.append(k)
            #Calculating standard deviation of every data set
            for value in dataSet:
                residue=residue+math.pow((value-meanVolume),2)
            stdDeviation=math.sqrt(residue/(nDataSets-1))
            error.append(stdDeviation)
            k+=1
        
        #Calculating value given by mean-field approximation
        
        #Function to calculate the volume given by the mean field theory
        def fun(v):
            return Constant().K()/deltaV*math.log(lambdaVol*(vb-v)/(v-vs))-(p-6*deltaE/deltaV*(v-vs)/deltaV)/T

        volumeUnitCell=Numericalmethods().newton(fun, 1.1*vs, 10e-9, 10e-32)
        meanFieldVolume=volumeUnitCell*math.pow(l,3)

        x2=x
        y2=[]
        for o in range(0, nDataSets):
            y2.append(meanFieldVolume) 
        plt.plot(x,y,x2,y2,linewidth=2.0)
        plt.title("Volume", fontsize=14)
        plt.xscale('log')
        plt.show()              
       
        


         

        











                    
                
                

        
