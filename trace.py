import numpy as np
import mylib 
from scipy.integrate import odeint

a1 = 1.0 / 4.0;
b1 = 3.0 / 32.0;
b2 = 9.0 / 32.0;
c1 = 1932.0 / 2197.0;
c2 = -7200.0 / 2197.0;
c3 = 7296.0 / 2197.0;
d1 = 439.0 / 216.0;
d2 = -8.0;
d3 = 3680.0 / 513.0;
d4 = -845.0 / 4104.0;
e1 = -8.0 / 27.0;
e2 = 2.0;
e3 = -3544.0 / 2565.0;
e4 = 1859.0 / 4104.0;
e5 = -11.0 / 40.0;
x1 = 25.0 / 216.0;
x2 = 0.0;
x3 = 1408.0 / 2565.0;
x4 = 2197.0 / 4104.0;
x5 = -1.0 / 5.0;
z1 = 16.0 / 135.0;
z2 = 0.0;
z3 = 6656.0 / 12825.0;
z4 = 28561.0 / 56430.0;
z5 = -9.0 / 50.0;
z6 = 2.0 / 55.0;

var0=[0,0,0,0,0,0]
h=1e-2
errmax=1e-1
errmin=1e-2
var=var0
t=0
mytime=[]
myvar=[]

fp=open('data.dat','w')

for ind in range(1000): 
    check=0
    vartmp=var
    diff=mylib.equation(vartmp,t)
    k1=h*diff
    vartmp = var + a1*k1

    diff=mylib.equation(vartmp,t)
    k2=h*diff
    vartmp = var + b1*k1 + b2*k2

    diff=mylib.equation(vartmp,t)
    k3=h*diff
    vartmp=var+c1*k1+c2*k2+c3*k3

    diff=mylib.equation(vartmp,t)
    k4=h*diff
    vartmp=var+d1*k1+d2*k2+d3*k3+d4*k4

    diff=mylib.equation(vartmp,t)
    k5=h*diff
    vartmp=var+e1*k1+e2*k2+e3*k3+e4*k4+e5*k5

    diff=mylib.equation(vartmp,t)
    k6=h*diff


    y  = var  + x1*k1  + x2*k2  + x3*k3  + x4*k4  + x5*k5 ;
    z  = var  + z1*k1  + z2*k2  + z3*k3  + z4*k4  + z5*k5  + z6*k6 ;
    print(y)
    print(z)
    for i in np.arange(6):
        if var[i]==0 and y[i]==0:
            err=errmin
        elif var[i]<10:
            err= np.abs(y[i]  - z [i])
        else:
            err= np.abs((y[i]  - z [i]) / max(np.abs(var[i] ), np.abs(y[i] ) ) );
        print(err)
        if err>errmax:
            check=1
        elif err<errmin and check!=1:
            check=-1


    if check==1:
        h=h/1.1
    elif check==-1:
        t=t+h
        var=y
        h=h*1.1
        print([ind,t,var])
        mytime.append(t)
        myvar.append(var)
        fp.write('%f %f %f %f %f %f %f\n'%(t,var[0],var[1],var[2],var[3],var[4],var[5]))

    else:
        t=t+h
        var=y
        print([ind,t,var])
        mytime.append(t)
        myvar.append(var)
        fp.write('%f %f %f %f %f %f %f\n'%(t,var[0],var[1],var[2],var[3],var[4],var[5]))
    if var[5]>6e7:
        break
