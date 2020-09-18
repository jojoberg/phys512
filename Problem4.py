import numpy as np
from matplotlib import pyplot as plt
import scipy.constants as const
from scipy import integrate

#Define Ring Field Equation  (from Griffiths Solutions - using u substitution u=cos(theta))
def ring(u,z):
    num=(z-bigR*u)
    denom=np.sqrt((bigR**2+z**2-2.0*bigR*z*u)**3)
    return num/denom


#Define my integrator (Simpsons)
def myintegrate(z,fun,xmin,xmax,tol):
    x=np.linspace(xmin,xmax,5)
    y=fun(x,z)
    dx=x[1]-x[0]
    area1=(y[0]+y[-1]+2*np.sum(y[1:-1]))*dx/2
    area2=(y[0]+4*np.sum(y[1::2])+2*np.sum(y[2:-1:2])+y[-1])*dx/3
    myerr=np.abs(area1-area2)
    if myerr<tol:
        return area2
    else:
        xm=0.5*(xmin+xmax)
        a1=myintegrate(z,fun,xmin,xm,tol/2)
        a2=myintegrate(z,fun,xm,xmax,tol/2)
        return a1+a2

#Set desired constants
bigR=1.0 #define radius of sphere
sigma=0.1 #define charge density

const=(2*np.pi*bigR**2*sigma)/(4.0*np.pi*const.epsilon_0)

#Integrate both ways 
z=np.linspace(0,10,501)
Ez1=np.zeros(len(z))
Ez2=np.zeros(len(z))
Etrue=np.zeros(len(z))

xmin=-1
xmax=1
tol=0.001
for zi in z: 
    if (zi-bigR)!=0:
        ans=myintegrate(zi,ring,xmin,xmax,tol)
        index=np.where(z==zi)[0][0]
        Ez1[index]=const*ans
    else: #my integrator gets stuck in an infinite loop if I include zi=bigR because there is a singularity there
        continue
    Ez2[index]=const*integrate.quad(ring,-1,1,args=(zi))[0]
    if zi>bigR: 
        Etrue[index]=2*const*(1/zi**2)
    else:
        Etrue[index]=0.0

#Print Errors in Fits
myerror=np.mean(np.abs(Ez1-Etrue)) #where is this large error coming from? graphs seem to agree

scipyerror=np.mean(np.abs(Ez2-Etrue))

print('The error in my integrator was approximately',myerror)
print('And the error in the scipy integrator was approximately',scipyerror)


#Plot results
plt.clf;
plt.plot(z,Ez1)
plt.plot(z,Ez2)
plt.plot(z,Etrue)
plt.savefig("integrate.png")




