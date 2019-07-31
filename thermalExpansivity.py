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
        line_new="Values for mean="+str(valuesForMean)
        fw.write(line_new+"\n")
        line_new ="Alpha_p"+"\n"+str(thermExpArray)
        fw.write(line_new+"\n")
        return [thermExp, stdDeviation]