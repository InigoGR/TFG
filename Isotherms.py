import math
import numpy as np
from Constant import *
import matplotlib.pyplot as plt

#System parameters
c=6
v0=2e-5
deltaV=0.5e-5
deltaE=1000
lambdaVol=0.2
T=250

#Array of volumes
v=np.linspace(2.001e-5, 2.499e-5, num=1000)
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
print(len(v))
print(len(p3))
#Turning to cubic meters
v=v/6.022e23
#Volume of the whole lattice
v=np.multiply(v,50*50*50)

#plt.plot(v,p1)
#plt.plot(v,p2)
plt.plot(p3,v)
plt.xlabel("P/Pa")
plt.ylabel("V/m^3")
plt.show()

