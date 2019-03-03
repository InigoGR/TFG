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
import numpy as np
import random
import math


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
            intEnergyVar=n*(inputData.getEb()-inputData.getEs())
            volumeChange=inputData.getVb()-inputData.getVs()
            #Changing volume
            lattice.changeVolume(lattice.getVolume()+volumeChange)
            #Changing cell volume
            lattice.changeCellVolume(x,y,z)
            #Changing energy
            lattice.changeEnergy(energy+intEnergyVar)

        #Case initialV=V+
        else:
            intEnergyVar=n*(inputData.getEs()-inputData.getEb())
            volumeChange=inputData.getVs()-inputData.getVb()      
        #Checking if the change is accepted
            if random.uniform(0,1)<math.pow(inputData.getL(),3)*(volume+volumeChange)/(volume)*\
                lattice.probArray[n]:
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
    #Writes in a file the  every meanSteps
    def getSystemEvolution(self, lattice, inputData, meanSteps):
        #Opening file and writting header
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "w")
        line_new ="T="+str(inputData.getT())
        fw.write(line_new+"\n")
        line_new="Volume"+"\t"+"Energy"+"\t"+"Enthalpy"
        fw.write(line_new+"\n")
        #Data arrays
        volumeValues=[]
        energyValues=[]
        enthalpyValues=[]
        #The lattice state is stored every meanSteps steps
        for i in range(0, int(inputData.getN()/meanSteps)):
            vol=0
            energy=0
            enthalpy=0
            #Changing cell volumes meanSteps times
            for j in range(0, meanSteps):
                LatticeHandler().changeVol(lattice, inputData)
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
        fw.close()
        return [volumeValues, energyValues, enthalpyValues]
    
    #Method to calculate the mean thermal expansivity from an array of mean values of Volume
    #and enthalpy an inputData object and the number of values to be used to calculate the mean
    #that will be used to get the response function. Returns an array with the mean value and the
    #standard deviation and writes the values (alpha values) used for the mean in the measurements file.
    def thermalExpansivity(self, volumes, enthalpy, inputData, valuesForMean):
        #Data sets
        volumeSets=[]
        enthalpySets=[]
        #Number of measurements
        nValues=len(volumes)
        #Dividing the arrays into data sets of nValues measurements
        #j--->Measurement counter, indicates location inside measurement arrays
        j=0
        for i in range(0, int(nValues/valuesForMean)):
            #Creating set arrays inside the main array
            volumeSets.append([])
            enthalpySets.append([])
            for k in range(0, valuesForMean):
                #Filling set arrays 
                volumeSets[i].append(volumes[j])
                enthalpySets[i].append(enthalpy[j])
                #Next component of the parameter array
                j+=1
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
            thermExpArray.append(1/meanVolume/Constant().K()/math.pow(inputData.getT(),2)*\
                (meanProduct-meanEnthalpy*meanVolume))
        thermExp=sum(thermExpArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((thermExpArray[i]-thermExp),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "a")
        line_new ="Alpha_p"+"\n"+str(thermExpArray)
        fw.write(line_new+"\n")
        return [thermExp, stdDeviation]
        

    #Method to calculate the mean isothermal compressibility from an array of mean values of Volume,
    #an inputData object and the number of values that have to be used to get the mean values.
    #Returns an array with the mean value and the standard deviation and writes the values used
    #for the mean in the measurements file. 
    def isothermalCompressibility(self, volumes, inputData, valuesForMean):
        #Data sets
        volumeSets=[]
        #Number of measurements
        nValues=len(volumes)
        #Dividing the arrays into data sets of nValues measurements
        #j--->Measurement counter, indicates location inside measurement arrays
        j=0
        for i in range(0, int(nValues/valuesForMean)):
            #Creating set arrays inside the main array
            volumeSets.append([])
            for k in range(0, valuesForMean):
                #Filling set arrays 
                volumeSets[i].append(volumes[j])
                #Next component of the parameter array
                j+=1
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
            isothCompressArray.append(1/meanVolume/Constant().K()/inputData.getT()*(meanVolumeSq-math.pow(meanVolume,2)))
        isothCompress=sum(isothCompressArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((isothCompressArray[i]-isothCompress),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "a")
        line_new ="Beta_t"+"\n"+str(isothCompressArray)
        fw.write(line_new+"\n")
        return [isothCompress, stdDeviation]
        

    #Method to calculate the residual heat capacity from an array of mean values of Volume,
    #Intermolecular energy and Hamiltonian and the amount of values that have to be used for 
    # the mean. Returns an array with the mean value and the standard deviation.
    def heatCapacity(self, volumes, energies, enthalpy, inputData, valuesForMean):
        #Data sets
        volumeSets=[]
        energySets=[]
        enthalpySets=[]
        #Number of measurements
        nValues=len(volumes)
        #Dividing the arrays into data sets of nValues measurements
        #j--->Measurement counter, indicates location inside measurement arrays
        j=0
        for i in range(0, int(nValues/valuesForMean)):
            #Creating set arrays inside the main array
            volumeSets.append([])
            enthalpySets.append([])
            energySets.append([])
            for k in range(0, valuesForMean):
                #Filling set arrays 
                volumeSets[i].append(volumes[j])
                enthalpySets[i].append(enthalpy[j])
                energySets[i].append(energies[j])
                #Next component of the parameter array
                j+=1
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

            heatCapArray.append(1/Constant().K()/math.pow(inputData.getT(),2)*(meanProductE-\
            meanEnergy*meanEnthalpy)+1/Constant().K()/math.pow(inputData.getT(),2)*(meanProductV-\
            meanVolume*meanEnthalpy)-math.pow(inputData.getL(),3)*Constant().K())
        heatCap=sum(heatCapArray)/nDataSets
        #Calculating deviation
        residue=0
        for i in range(0,nDataSets):
            residue=residue+math.pow((heatCapArray[i]-heatCap),2)
        stdDeviation=math.sqrt(residue/(nDataSets-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "a")
        line_new ="C_p"+"\n"+str(heatCapArray)
        fw.write(line_new+"\n")
        return [heatCap, stdDeviation]

    def isentropicCompressibility(self, inputData, meanSteps, valuesForMean):
        #Creating data arrays
        volumes=[]
        beta_t=[]
        alpha_p=[]
        c_p=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(inputData.getT())
        #Retrieving data from the file
        fr=open(fileName, "r")
        #Getting file lines
        lines=fr.readlines()
        #Line number, skipping two first lines
        nLine=2
        #Reading all the measurements before the response functions
        for j in range(0, int(inputData.getN()/meanSteps/valuesForMean)):
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
            isentCompress=beta_t[i]-inputData.getT()*meanVolumes[i]*math.pow(alpha_p[i],2)/c_p[i]
            isentCompressArray.append(isentCompress)
        isentropicCompressibility=sum(isentCompressArray)/len(isentCompressArray)
        #Calculating deviation
        residue=0
        for i in range(0,len(isentCompressArray)):
            residue=residue+math.pow((isentCompressArray[i]-isentropicCompressibility),2)
        stdDeviation=math.sqrt(residue/(len(isentCompressArray)-1))
        #Writting results in file (multiple values not mean value)
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "a")
        line_new ="K_s"+"\n"+str(isentCompressArray)
        fw.write(line_new+"\n")
        return [isentropicCompressibility, stdDeviation]

        

        


         

        











                    
                
                

        
