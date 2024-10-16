import csv
import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import numpy as np
import ast

d = 0.0005

xf = 0.2

n = int(7) #int(input("How many frequencies?"))
Name = 'Diff21_MT_' #input("What is the designate for the graphs?")

RI = float(70) #float(input("What is the internal resistivity? (Ohm-centimetre.)"))
CM = float(0.0021) #float(input("What is the DC capacitance per unit area of membrane? (micro F-centimetre squared.) "))
CMin = float(0.00021) #float(input("What is the minimum capacitance of the membrane? "))

def mtcond(f):
    a = 69*10**(-10)
    b = -2*69*10**(-5)
    return a*f**2 + b*f +70

def cap_mem(f,CM2,CINF2):
    c_inf = CINF2
    c_del = CM2 - c_inf
    a = c_inf + c_del/(1+2*np.pi*f*2e-3)
    return a*10**(-6)

def constants1(CM1,CINF,RI1,f):
    ri = 4*mtcond(f)/(np.pi*d**2)
    tau  = np.sqrt(np.pi*f*ri*cap_mem(f,CM1,CINF))
    return(tau, ri)

def diffsol(f,x,t,tau,R):
    t2 = np.exp(-tau*x)
    t3 = np.sin(2*np.pi*f*t - tau*x)
    return t2*t3

for i in range(3,n):
    a = ['11(a)','11(b)','11(c)','11(d)','(e)']
    f = 10**(i)
    lines = ["k-","k--","k-.","k:"]
    C = constants1(CM,CMin,RI,f)
    tau1 = C[0]
    R1 = C[1]
    tf = 0.002
    def solution5(w,x,t):
        S1 = diffsol(w,x,t,tau1,R1)
        return S1
    xspace = np.linspace(0,xf,10001)
    I1 = abs(solution5(f,xspace,tf))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(xspace,I1,lines[0], label = 'Cable equation')
    plt.xlabel("Distance from source", fontsize = 12)  # add X-axis label
    plt.ylabel("Current", fontsize = 12)  # add Y-axis label
    plt.title(a[i-3]+" Current at time "+str(tf)+"; frequency = "+str(f), fontsize = 12)
    plt.legend(loc='best')
    plt.savefig(Name+"Freq="+str(f)+"xf="+str(xf)+"_tf="+str(tf)+".png")
    plt.show()