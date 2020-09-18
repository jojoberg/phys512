import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate

#Import Data and Assign Voltage and Temperature Lists
mydata=np.loadtxt("./PS1/lakeshore.txt")
temp=mydata[:,0]
volt=mydata[:,1]

#Plot Data
plt.clf();
plt.plot(volt,temp,"*")

vInterp=np.linspace(volt[1],volt[-1],2001)

#Interpolate using Polynomial
def poly_interp(volt,temp,x):
    poly=0.0
    for n in range(0,len(volt),4):
        min=n
        max=n+4
        if n==(len(volt)-4): max=n+3
        if x<=volt[min] and x>=volt[max]:
            for i in range(min,max):
                x_use=np.append(volt[min:i],volt[i+1:max])
                x0=volt[i]
                y0=temp[i]
                denom=np.prod((x0-x_use))
                num=1.0
                for xi in x_use:
                   num=num*(x-xi)
                poly=poly+y0*num/denom
        else:
            continue
    return poly

#Plot Polynomial Interpolation
polyPlot=np.zeros(len(vInterp))
for j in range(len(vInterp)):
   polyPlot[j]=poly_interp(volt,temp,vInterp[j])
plt.plot(vInterp,polyPlot)

#Interpolate Using Cubic Spline to compare
spln=interpolate.splrep(volt[::-1],temp[::-1])
cubic=interpolate.splev(vInterp,spln)

#Ask user for arbitrary voltage value and return temperature interpolation
vUser=float(input("Please input the desired voltage:\n"))
if vUser <= volt[-1] or vUser >= volt[0]:
    print('Warning: desired voltage lies outside of known data set.')
tOutP=poly_interp(volt,temp,vUser)
print("The temperature calculated using the polynomial fit is",tOutP,"Kelvin.")
tOutC=interpolate.splev(vUser,spln)
print("And the temperature calculated using the cubic spline fit is",tOutC, "Kelvin.")

#Calculate approximate error by comparing even and odd points
tEven=temp[::2]
vEven=volt[::2]
tOdd=temp[1::2]
vOdd=volt[1::2]

#Error in Spline
splnO=interpolate.splrep(vOdd[::-1],tOdd[::-1])
errorC=np.mean(np.absolute(interpolate.splev(vEven,splnO)-tEven))
plt.plot(vOdd,interpolate.splev(vOdd,splnO))

#Error in Polynomial
polyOdd=np.zeros(len(vEven))
for i in range(len(vEven)):
    polyOdd[i]=poly_interp(vOdd,tOdd,vEven[i])
plt.plot(vEven,polyOdd)
errorP=np.mean(np.abs(polyOdd-tEven))

print("The approximate error in the polynomial fit is",errorP)
print("The approximate error in the cubic spline fit is",errorC)

plt.savefig('lakeshore_out.png')




