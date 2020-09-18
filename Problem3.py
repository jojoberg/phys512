import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import splev, splrep
import math

#Define constants and clear plot
pi=math.pi
plt.clf();

#Cosine Function
x=np.linspace(-pi/2,pi/2,12)
y=np.cos(x)
xx=np.linspace(x[0],x[-1],1001)
ytruecos=np.cos(xx)

#Polynomial Interpolation
poly=0.0
for i in range(len(x)):
    x_use=np.append(x[:i],x[i+1:])
    x0=x[i]
    y0=y[i]
    denom=np.prod((x0-x_use))
    num=1.0
    for xi in x_use:
        num=num*(xx-xi)
    poly=poly+y0*num/denom
    
print('The error in the polynomial interpolation of cosine is',np.mean(np.abs(poly-ytruecos)))
plt.plot(xx,poly)

#Cubic Spline Interpolation
spln=splrep(x,y)
cubic=splev(xx,spln)

print('The error in the cubic spline interpolation of cosine is',np.mean(np.abs(cubic-ytruecos)))
plt.plot(xx,cubic)

#Rational Function Interpolation
def rat_eval(p,q,x):
    top=0
    for i in range(len(p)):
        top=top+p[i]*x**i
    bot=1
    for i in range(len(q)):
        bot=bot+q[i]*x**(i+1)
    return top/bot

def rat_fit(x,y,n,m):
    mat=np.zeros([n+m-1,n+m-1])
    for i in range(n): #0 to n
        mat[:,i]=x**i #define ith element for every row
    for i in range(1,m): #1 to m because q starts at 1
        mat[:,i-1+n]=-y*x**i #define i-1+nth element for every row
    pars=np.dot(np.linalg.inv(mat),y)
    p=pars[:n]
    q=pars[n:]
    return p,q

n=6  #n+m must equal 13 (number of points + 1)
m=7
p,q=rat_fit(x,y,n,m)
rational=rat_eval(p,q,xx)

print('The error in the rational function interpolation of cosine is',np.mean(np.abs(rational-ytruecos)))
plt.plot(xx,rational)

#Plot Points and Save Figure
plt.plot(x,y,"*")
plt.savefig('cosine_out.png')


#Lorentzian Function
def lorentz(x):
    return 1/(1+x**2)

plt.clf();

x2=np.linspace(-1,1,12)
xx2=np.linspace(x2[0],x2[-1],1001)
y2=lorentz(x2)
ytruelorentz=lorentz(xx2)

#Polynomial Interpolation (take 2)
poly2=0.0
for i in range(len(x2)):
    x_use2=np.append(x2[:i],x2[i+1:])
    x02=x2[i]
    y02=y2[i]
    denom=np.prod((x02-x_use2))
    num=1.0
    for xi in x_use2:
        num=num*(xx2-xi)
    poly2=poly2+y02*num/denom

print('\nThe error in the polynomial interpolation of the lorentzian is',np.mean(np.abs(poly2-ytruelorentz)))
plt.plot(xx2,poly2)

#Cubic Spline (again)
spln2=splrep(x2,y2)
cubic2=splev(xx2,spln2)

print('The error in the cubic spline interpolation of the lorentzian is',np.mean(np.abs(cubic2-ytruelorentz)))
plt.plot(xx2,cubic2)

#Rational Function - np.inv(mat)

n=6  #n+m must equal 13 (number of points + 1)
m=7
p2,q2=rat_fit(x2,y2,n,m)
rational2=rat_eval(p2,q2,xx2)

#Rational Function - np.pinv(mat)
def rat_fit2(x,y,n,m):
    mat=np.zeros([n+m-1,n+m-1])
    for i in range(n): #0 to n
        mat[:,i]=x**i #define ith element for every row
    for i in range(1,m): #1 to m because q starts at 1
        mat[:,i-1+n]=-y*x**i #define i-1+nth element for every row
    pars=np.dot(np.linalg.pinv(mat),y)
    p=pars[:n]
    q=pars[n:]
    return p,q

p3,q3=rat_fit2(x2,y2,n,m)
rational3=rat_eval(p3,q3,xx2)

print('The error in the rational function interpolation of the lorentzian is',np.mean(np.abs(rational2-ytruelorentz)))

print('The error in the rational function interpolation with np.pinv() is', np.mean(np.abs(rational3-ytruelorentz)))
plt.plot(xx2,rational3)


#Plot Lorentzian Points and Save Figure
plt.plot(x2,y2,"*")
plt.savefig('lorentzian_out.png')

#Plot Lorentzian with bad rational function fit
plt.clf();
plt.plot(x2,y2,'*')
plt.plot(xx2,rational2)
plt.savefig('lorentz_rational_inv.png')

#Trying Lorentzian rational function fit with four points instead of 12
x3=np.linspace(-1,1,4)
xx3=np.linspace(x3[0],x3[-1],1001)
y3=lorentz(x3)
ytruelorentz2=lorentz(xx3)

n=2
m=3

p4,q4=rat_fit(x3,y3,n,m)
rational4=rat_eval(p4,q4,xx3)
print('\nThe error in the four-point rational function interpolation of the lorentzian is',np.mean(np.abs(rational4-ytruelorentz2)))

p5,q5=rat_fit2(x3,y3,n,m)
rational5=rat_eval(p5,q5,xx3)
print('The error in the four-point rational function interpolation with np.pinv() is', np.mean(np.abs(rational4-ytruelorentz2)))

plt.clf();
plt.plot(x3,y3,'*')
plt.plot(xx3,rational4)
plt.plot(xx3,rational5)
plt.savefig('lorentz_ratfit_4.png')





