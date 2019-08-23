'''
Created on 20 dic. 2018

@author: Iñigo
@version: 04/08/19

'''

#LatticeCounter.py
#Class responsible for the operations on the lattice simulation, theoretical calculations and result management.

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
    
    #Method to get the state of the neighbors of each cell and register it in the neighborHistogram variable
    def neighborStatistics(self, lattice):
        L=lattice.getInputData().getL()
        for x in range(0, L):
            for y in range(0, L):
                for z in range(0, L):
                    neighbors=lattice.getNeighbors()[x,y,z]  #Storing neighbor coords
                    n=0 #Number of neighbors in V+ state
                    for neighborCoord in neighbors:  #Getting coordinates in the lattice for every neighbor  
                        n=n+lattice.getLattice()[neighborCoord[0], neighborCoord[1], neighborCoord[2]]  
                    lattice.changeHistogram(n) 
    
    #Method that tries a random cell volume change and checks if the system accepts, making
    #the changes in volume, energy and configuration when necessary.
    def changeVol(self, lattice, inputData):
        #Choosing random cell in the lattice
        x=random.randint(0,inputData.getL()-1)
        y=random.randint(0,inputData.getL()-1)
        z=random.randint(0,inputData.getL()-1)
        volumeChange=0  #Parameter to store the change in volume
        tempLattice=lattice.getLattice()    #Storing state of the lattice (matrix of 0s and 1s)
        neighbors=lattice.getNeighbors()[x,y,z]  #Storing neighbor coords
        volume=lattice.getVolume()  #Storing volume of the lattice
        energy=lattice.getEnergy()  #Storing energy of the lattice
        #Counting number of neighbors with V+
        n=0 #parameter to count the number of neighbors in V+ state
        #Counting V+ neighbors
        for neighborCoord in neighbors:  #Getting coordinates in the lattice for every neighbor  
            n=n+tempLattice[neighborCoord[0], neighborCoord[1], neighborCoord[2]]   #Cumulative counting (V+=1, V-=0)
        #Getting intermolecular energy change
        #Case initialV=V-
        if tempLattice[x,y,z]==0:
            if random.uniform(0,1)<lattice.probArrayMinus[n]:    #The change of state of the cell is accepted
            
                intEnergyVar=n*(inputData.getEb()-inputData.getEs())    #Calculating variation of the energy
                volumeChange=inputData.getVb()-inputData.getVs()    #Calculating variation of the volume   
                #Changing volume
                lattice.changeVolume(volume+volumeChange)
                #Changing cell volume
                lattice.changeCellVolume(x,y,z)
                #Changing energy
                lattice.changeEnergy(energy+intEnergyVar)
                #Changing cell type counter
                lattice.changeMinusToPlus()

        #Case initialV=V+
        else:      
        #Checking if the change is accepted
            if random.uniform(0,1)<lattice.probArrayPlus[n]: #The change of state of the cell is accepted
                
                intEnergyVar=n*(inputData.getEs()-inputData.getEb())    #Calculating variation of the energy
                volumeChange=inputData.getVs()-inputData.getVb()    #Calculating variation of the volume
                #Changing volume
                lattice.changeVolume(lattice.getVolume()+volumeChange)
                #Changing cell volume
                lattice.changeCellVolume(x,y,z)
                #Changing energy
                lattice.changeEnergy(energy+intEnergyVar)
                #Changing cell type counter
                lattice.changePlusToMinus()
    
    #Method to get the volume, enthalpy and intermolecular energy in intervals separated
    #by the indicated amount of steps with a total of MonteCarlo steps indicated by the inputData.
    #Returns three arrays as shown
    #[[VolumeValues],
    #[EnergyValues],
    #[EnthalpyValues]]
    #Writes in a file the measurements of the system every meanSteps
    def getSystemEvolution(self, lattice, inputData, meanSteps):
        #Opening file and writting header
        fileName="Measurements_T"+str(inputData.getT())
        fw=open(fileName, "w")
        #Opening file to store probabilities
        fileName="probabilities"
        fw2=open(fileName, "a")
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
        #Writting V+ volume
        line_new="Big cell volume="+str(inputData.getVb())
        fw.write(line_new+"\n")
        #Writting V- volume
        line_new="Small cell volume="+str(inputData.getVs())
        fw.write(line_new+"\n")
        #Writting big free volume
        line_new="Big free volume="+str(inputData.getFbv())
        fw.write(line_new+"\n")
        #Writting small free volume
        line_new="Small free volume="+str(inputData.getFsv())
        fw.write(line_new+"\n")
        #Writting energy of V+ V+ interaction
        line_new="+ Interaction energy="+str(inputData.getEb())
        fw.write(line_new+"\n")
        #Writting energy of V+ V- or V- V- interaction
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
        #Number of volumes to be registered (Total amount of steps divided by the steps taken between measurements)
        nVols=int(inputData.getN()/meanSteps)

        print("Taking measurements")

        for i in range(0, nVols):
            vol=0
            energy=0
            enthalpy=0
            #Changing cell volumes meanSteps times before taking a measurement
            for j in range(0, meanSteps):
                LatticeHandler().changeVol(lattice, inputData)  #Try to change state of one cell
            #Getting state of the lattice
            vol=lattice.getVolume() #Getting volume of the lattice
            energy=lattice.getEnergy()  #Getting energy of the lattice
            enthalpy=energy+inputData.getP()*vol    #Getting enthalpy of the lattice
            #Adding measurements to the data arrays
            volumeValues.append(vol)    #Adding measured volume to measurements array
            energyValues.append(energy) #Adding measured energy to measurements array
            enthalpyValues.append(enthalpy) #Adding measured enthalpy to measurements array
            #Saving parameters in the measurements file
            line_new=str(vol)+"\t"+str(energy)+"\t"+str(enthalpy)
            fw.write(line_new+"\n")
            #Simulation progress indicator
            if math.fmod(counter/nVols*100,1.0)==0: #Checking progress of the measurement phase
                print(str(counter/(nVols)*100)+"%")
                LatticeHandler().neighborStatistics(lattice)   #Updating neighbor histogram array
            counter+=1

        fw.close()  #Closing file
        fw2.write("Probabilities for T="+str(inputData.getT())+"\n")
        fw2.write("+ to -"+"\n")
        fw2.write(str(lattice.probArrayPlus)+"\n")
        fw2.write("- to +"+"\n")
        fw2.write(str(lattice.probArrayMinus)+"\n"+"\n")
        fw2.close()
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
            #Introducing data set into data array in molar units
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
        

    #Method to calculate the heat capacity for a given temperature from the data files generated
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
            #Residual heat capacity
            resHeatCap=1/Constant().K()/math.pow(T,2)*(meanProductE-\
            meanEnergy*meanEnthalpy)+1/Constant().K()/math.pow(T,2)*(meanProductV-\
            meanVolume*meanEnthalpy)
            #Turning to molar units
            heatCapArray.append(resHeatCap/math.pow(l,3)*6.022e23+5/2*6.022e23*Constant().K())
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
            isentCompress=beta_t[i]-T*meanVolumes[i]*\
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

    

    #Method that shows the evolution of the volume for a given temperature using the corresponding
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
        #Calculating deltaE and deltaV
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
            return Constant().K()/deltaV*math.log(lambdaVol*(vb-v)/(v-vs))-(p-6*deltaE/deltaV*\
                (v-vs)/deltaV)/T

        volumeUnitCell=Numericalmethods().newton(fun, (vb+vs)/2, 10e-9, 10e-32)
        meanFieldVolume=volumeUnitCell*math.pow(l,3)

        x2=x
        y2=[]
        for o in range(0, nDataSets):
            y2.append(meanFieldVolume) 
        plt.errorbar(x,y,error)
        plt.plot(x2,y2,linewidth=2.0, linestyle='dashed')
        plt.title("Volume-Time", fontsize=14)
        plt.xscale('log')
        plt.xlabel("t/Steps")
        plt.ylabel("V/Jm^3")
        plt.show()              
       
    #Method that shows the evolution of the volume as the temperature increases and compares it to the theoretical value
    def volEvoTemps(self, iniT, finT, tempIncrement):
        Volumes=[]
        error=[]
        temps=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(iniT)
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
        #Number of measurements
        nMeasurements=int(n/meanSteps)
        #File to write simulation mean energy values
        writeFileName="Volume_Temperature"
        fw=open(writeFileName, "w")
        #Header of the file
        fw.write("Pressure="+str(p)+" Pa"+"\n")
        fw.write("L="+str(l)+"\n")
        fw.write("dV="+str(deltaV*6.022e23)+" m^3/mol"+"\n")
        fw.write("dE="+str(deltaE*6.022e23)+" J/mol"+"\n")
        fw.write("Lambda="+str(lambdaVol)+"\n")
        #Header of the file
        fw.write("T/K"+"\t"+"V/m^3"+"\n")
        for T in range(iniT, finT+tempIncrement, tempIncrement):
            temps.append(T)
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
            #Number of measurements
            nMeasurements=int(n/meanSteps)
            #Reading all the measurements before the response functions
            nLine=13
            volumeSet=[]
            for j in range(0, nMeasurements):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                #Changing line
                nLine+=1
            #Changing to molar units
            volumeSet=np.multiply(volumeSet, 6.022e23/math.pow(l,3))
            #Introducing data set into data array
            Volumes.append(statistics.mean(volumeSet))
            error.append(statistics.stdev(volumeSet))
            fw.write(str(T)+"\t"+str(statistics.mean(volumeSet))+"  +-  "+str(statistics.stdev(volumeSet))+"\n")

        #Calculating theoretical values
        #Array of volumes
        v=np.linspace(2.11e-5, 2.4999999e-5, num=10000)
        
        #Performing unit change
        vs=vs*6.022e23   
        deltaV=deltaV*6.022e23    
        deltaE=deltaE*6.022e23    

        #Array of temperatures
        T=[]
        for volume in v:
            T.append((p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs)))



        #Plotting results
        plt.errorbar(temps,Volumes,error,marker='o', linewidth=2.0, label='Sim')
        plt.plot(T, v, linewidth=2.0, linestyle='dashed', label='Theory')
        plt.legend()
        plt.title("Volume-Temperature", fontsize=14)
        plt.xlabel("T/K")
        plt.ylabel(r"$V\ /\ m^3\ mol^-1$")
        plt.show() 


    #Method that shows the evolution of the volume as the pressure increases
    def volEvoPress(self, iniP, finP, pressureIncrement):
        Volumes=[]
        error=[]
        pressures=[]
        #File to write simulation mean volume values
        writeFileName="Volume_Pressure"
        fw=open(writeFileName, "w")
        #Data file corresponding to the needed parameters
        fileName="Measurements_P"+str(iniP)
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
        line=lines[1].split("=")
        T=float(line[1])
        line=lines[11].split("=")
        meanSteps=float(line[1])
        lambdaVol=fvb/fvs
        deltaV=vb-vs
        deltaE=es-eb
        #Header of the file
        #Number of measurements
        nMeasurements=int(n/meanSteps)
        #File to write simulation mean energy values
        writeFileName="Pressure_Volume"
        fw=open(writeFileName, "w")
        #Header of the file
        fw.write("Temperature="+str(T)+" K"+"\n")
        fw.write("L="+str(l)+"\n")
        fw.write("dV="+str(deltaV*6.022e23)+" m^3/mol"+"\n")
        fw.write("dE="+str(deltaE*6.022e23)+" J/mol"+"\n")
        fw.write("Lambda="+str(lambdaVol)+"\n")
        fw.write("P/Pa"+"\t"+"V/m^3 mol^-1"+"\n")
        for P in range(iniP, finP, pressureIncrement):
            pressures.append(P)
            #Data file corresponding to the needed parameters
            fileName="Measurements_P"+str(P)
            #Retrieving data from the file
            fr=open(fileName, "r")
            #Getting file lines
            lines=fr.readlines()
            #Reading inputData parameters
            line=lines[0].split("=")
            l=float(line[1])
            line=lines[1].split("=")
            T=float(line[1])
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
            #Number of measurements
            nMeasurements=int(n/meanSteps)
            #Reading all the measurements before the response functions
            nLine=13
            volumeSet=[]
            for j in range(0, nMeasurements):
                values=lines[nLine].split("\t")
                #Introducing read volume (values[0]) in the data set
                volumeSet.append(float(values[0]))
                #Changing line
                nLine+=1
            #Changing to molar units
            volumeSet=np.multiply(volumeSet, 6.022e23/math.pow(l,3))
            #Introducing data set into data array
            Volumes.append(statistics.mean(volumeSet))
            error.append(statistics.stdev(volumeSet))
            fw.write(str(P)+"\t"+str(statistics.mean(volumeSet))+"  +-  "+str(statistics.stdev(volumeSet))+"\n")
         

        #System parameters
        c=6
        v0=2e-5
        deltaV=0.5e-5
        deltaE=1000
        lambdaVol=0.2
        T=250

        #Array of volumes
        v=np.linspace(2.02e-5, 2.45e-5, num=1000)
        #Array of presures
        p1=[]
        p2=[]
        p3=[]
        k=0
        for volume in v:
            p1.append(T*8.314/deltaV*math.log(lambdaVol*(v0+deltaV-volume)/(volume-v0)))
            p2.append(c*deltaE/deltaV*(volume-v0)/deltaV)
            p3.append(p1[k]+p2[k])
            k+=1

        #Plotting results
        plt.errorbar(pressures,Volumes,error, marker='o', linewidth=2.0, label='sim')
        plt.plot(p3,v, linestyle='dashed', label='theory')
        plt.title("Volume-Pressure", fontsize=14)
        plt.xlabel("P/Pa")
        plt.ylabel(r"$V\ /\ m^3\ mol^-1$")
        plt.legend()
        plt.show() 

    
    #Method that shows the evolution of the energy for a given temperature using the corresponding
    #data file 'Measurements_Tx'. It needs the number of values for the mean
    def energyEvo(self, T, valuesForMean):
        #Creating data arrays
        energies=[]
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
            #Filling energy data set
            energySet=[]
            for i in range(0, valuesForMean):
                values=lines[nLine].split("\t")
                #Introducing read energy (values[1]) in the data set
                energySet.append(float(values[1]))
                #Changing line
                nLine+=1
            #Introducing data set into data array
            energies.append(energySet)
        #Plotting mean values of the energy sets
        x=[]
        y=[]
        error=[]

        #Calculating deviation to get error of every energy
        #k counter
        k=0
        for dataSet in energies:
            #Initial residue 
            residue=0
            #Calculating mean
            meanEnergy=statistics.mean(dataSet)
            #Adding mean to data array
            y.append(meanEnergy)
            x.append(k)
            #Calculating standard deviation of every data set
            for value in dataSet:
                residue=residue+math.pow((value-meanEnergy),2)
            stdDeviation=math.sqrt(residue/(nDataSets-1))
            error.append(stdDeviation)
            k+=1
        
        #Calculating value given by mean-field approximation
        
        #Function to calculate the volume given by the mean field theory
        def fun(v):
            return Constant().K()/deltaV*math.log(lambdaVol*(vb-v)/(v-vs))-(p-6*deltaE/deltaV*\
                (v-vs)/deltaV)/T

        volumeUnitCell=Numericalmethods().newton(fun, (vb+vs)/2, 10e-9, 10e-32)
        energyUnitCell=-3*deltaE*math.pow((volumeUnitCell-vs)/deltaV,2)
        meanFieldEnergy=energyUnitCell*math.pow(l,3)
        


        x2=x
        y2=[]
        for o in range(0, nDataSets):
            y2.append(meanFieldEnergy) 
        plt.errorbar(x,y,error, label='Sim')
        plt.plot(x2,y2,linewidth=2.0, linestyle='dashed', label='Theory')
        plt.title("Energy-Time", fontsize=14)
        plt.xscale('log')
        plt.xlabel("t/Steps")
        plt.ylabel("E/J")
        plt.legend()
        plt.show()
        


    #Method that shows the evolution of the energy as the temperature increases. Creates a file with the simulation results.
    def energyEvoTemps(self, iniT, finT, tempIncrement):
        Energies=[]
        error=[]
        temps=[]
        #Data file corresponding to the needed parameters
        fileName="Measurements_T"+str(iniT)
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
        #Number of measurements
        nMeasurements=int(n/meanSteps)
        #File to write simulation mean energy values
        writeFileName="Energy_Temperature"
        fw=open(writeFileName, "w")
        #Header of the file
        fw.write("Pressure="+str(p)+" Pa"+"\n")
        fw.write("L="+str(l)+"\n")
        fw.write("dV="+str(deltaV*6.022e23)+" m^3/mol"+"\n")
        fw.write("dE="+str(deltaE*6.022e23)+" J/mol"+"\n")
        fw.write("Lambda="+str(lambdaVol)+"\n")
        fw.write("T/K"+"\t"+"Energy/J mol^-1"+"\n")
        #calculating mean energy for every simulation temperature
        for T in range(iniT, finT+tempIncrement, tempIncrement):
            temps.append(T)
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
            #Number of measurements
            nMeasurements=int(n/meanSteps)
            #Reading all the measurements before the response functions and introducing energies into array
            nLine=13
            energySet=[]
            for j in range(0, nMeasurements):
                values=lines[nLine].split("\t")
                #Introducing read energy (values[1]) in the data set
                energySet.append(float(values[1]))
                #Changing line
                nLine+=1
            #Turning to molar units
            energySet=np.multiply(energySet, 6.022e23/math.pow(l,3))
            #Introducing data set mean into data array
            Energies.append(statistics.mean(energySet))
            error.append(statistics.stdev(energySet))
            #Printing mean energies in file
            fw.write(str(T)+"\t"+str(statistics.mean(energySet))+"  +-  "+str(statistics.stdev(energySet))+"\n")
        

        
        #Calculating theoretical values
        #Array of volumes
        v=np.linspace(2.11e-5, 2.4999999e-5, num=10000)
        
        #Performing unit change
        vs=vs*6.022e23   
        deltaV=deltaV*6.022e23  
        deltaE=deltaE*6.022e23    

        #Array of temperatures
        T=[]
        for volume in v:
            T.append((p-6*deltaE/deltaV*(volume-vs)/deltaV)/8.31*deltaV/math.log(lambdaVol*(vs+deltaV-volume)/(volume-vs)))
        

        #Performing unit change
        vs=vs/6.022e23
        v=v/6.022e23   
        deltaV=deltaV/6.022e23  
        deltaE=deltaE/6.022e23    

        #Energy per unit cell
        energyUnitCell=np.multiply(-3*deltaE, np.power((v-vs)/deltaV,2))
        #Energy of the whole lattice
        meanFieldEnergy=np.multiply(energyUnitCell, 6.022e23)
       
        
        #Plotting results
        plt.errorbar(temps,Energies,error, marker='o', linewidth=2.0, label='Sim')
        plt.plot(T,meanFieldEnergy,linewidth=2.0, linestyle='dashed', label='Theory')
        plt.legend()
        plt.title("Energy-Temperature", fontsize=14)
        plt.xlabel("T/K")
        plt.ylabel(r"$E\ /\ J\ mol^-1$")
        plt.show()

    #Method for calculating theoretically the isothermal compressibility
    def isothermalCompressibilityTheory(self, p, iniT=200, finT=300, v0=3.321155762e-29,
        deltaV=8.30288941e-30, deltaE=1.660577881e-21, lambdaVol=0.2):
        
        #Calculating 1/V*(dv/dT) at cte pressure for all temperatures
        Kt=[]
        TArray=[]
        for T in range(iniT, finT):
            #Function to calculate the volume given by the mean field theory for a given temperature
            def volumeForP(p):
                def fun(v):
                    return Constant().K()/deltaV*math.log(lambdaVol*(v0+deltaV-v)/(v-v0))-(p-6*\
                        deltaE/deltaV*(v-v0)/deltaV)/T

                volumeUnitCell=Numericalmethods().newton(fun, 1.1*v0, 10e-9, 10e-34)
                return volumeUnitCell
            #Volume for the established pressure
            v=volumeForP(p)
            #Value of the derivative
            derivative=Numericalmethods().df(volumeForP, p, 1)
            #Calculation of the response function
            KtTemp=-derivative/v
            Kt.append(KtTemp)
            TArray.append(T)

        #Plotting results
        plt.plot(TArray, Kt,linewidth=2.0)
        plt.title(r'$\beta_T$', fontsize=14)
        plt.show()  
    
    #Method for calculating theoretically the thermal expansivity
    def thermalExpansivityTheory(self, p, iniT=200, finT=300, v0=3.321155762e-29,
        deltaV=8.30288941e-30, deltaE=1.660577881e-21, lambdaVol=0.2):
        #Function to calculate the volume given by the mean field theory for a given temperature
        def volumeForT(T):
            def fun(v):
                return Constant().K()/deltaV*math.log(lambdaVol*(v0+deltaV-v)/(v-v0))-(p-6*\
                    deltaE/deltaV*(v-v0)/deltaV)/T

            volumeUnitCell=Numericalmethods().newton(fun, (v0+deltaV)/2, 10e-7, 10e-37)
            return volumeUnitCell
        
        #Calculating 1/V*(dv/dT) at cte pressure for all temperatures
        alpha=[]
        TArray=[]
        for T in range(iniT, finT):
            v=volumeForT(T)
            derivative=Numericalmethods().df(volumeForT, T, 1)
            alphaTemp=derivative/v
            alpha.append(alphaTemp)
            TArray.append(T)

        #Plotting results
        plt.plot(TArray, alpha,linewidth=2.0)
        plt.title(r'$\alpha_P$', fontsize=14)
        plt.show()  

    #Method to plot the neighbor histogram
    def plotNeighborHist(self, T):
        fileName="Neighbor_Histogram"
        #Retrieving data from the file
        fr=open(fileName, "r")
        lines=fr.readlines()
        nline=0 #Number of the line being read
        #Looking for part of the file containing the data corresponding to the desired temperature
        for line in lines:
            if line == ("T="+str(T)+"\n"):  #Comparing read line to expected header line
                x=[]    #Number of V+ neighbors
                y=[]    #Number of cells in that situation
                for i in range(0,7):    
                    line=lines[nline+i+2].split("\t")   #+2 since there are two lines before the data
                    line=line[1].replace("\n", "")
                    x.append(i)
                    y.append(int(line))
                plt.bar(x, y)
                plt.title("T="+str(T))
                plt.show()
            nline+=1

    
    


        
            


         

        











                    
                
                

        
