import numpy as np

#Exp(x)
expvals=np.linspace(-5,-1,10)
x0=1
#d(exp(x))/dx = exp(x)
#true derivative is just exp(x0)
truth=np.exp(x0)

print("\ndx, approximate derivative, and error in exp(x)")
for myexp in expvals:
    #exp(x) error vs dx
    dx=10**myexp
    f1=np.exp(x0+dx)
    fm=np.exp(x0-dx)
    f2=np.exp(x0+2*dx)
    fm2=np.exp(x0-2*dx)
    deriv1=(f1-fm)/(2*dx)  #make the derivative from (f(x+dx)-f(x-dx))/2dx
    deriv2=(f2-fm2)/(4*dx) #make the derivative out of (f(x+2dx)-f(x-2dx))/4dx
    derivf=(4*deriv1-deriv2)/3.0
    print(myexp,'\t',derivf,'\t',np.abs(derivf-truth))

#Exp(0.01x)
expvals2=np.linspace(-5,-1,10)
x1=0.01
#d(exp(0.01x))/dx = 0.01exp(0.01x)
#true derivative is just 0.01exp(0.01x0)
truth2=(x1)*np.exp(x1)
print("dx, approximate derivative, and error in exp(0.01x)")
for myexp2 in expvals2:
    #exp(0.01x) error vs dx
    dx2=10**myexp2
    fn1=np.exp(x1+dx2)
    fnm=np.exp(x1-dx2)
    fn2=np.exp(x1+2*dx2)
    fnm2=np.exp(x1-2*dx2)
    derivn1=(fn1-fnm)/(2*dx2)  #make the derivative from (f(x+dx)-f(x-dx))/2dx
    derivn2=(fn2-fnm2)/(4*dx2) #make the derivative out of (f(x+2dx)-f(x-2dx))/4dx
    derivnf=(4*derivn1-derivn2)/3.0
    print(myexp2,'\t',derivnf,'\t',np.abs(derivnf-truth2))



