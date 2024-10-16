import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import cmath
from itertools import cycle
import csv

name = input("Designate for graphs: ")
d = 0.0005 #float(input("What is the axon diameter (in cm)? "))
h = 0.0005 #float(input("What is the thickness of the extracellular layer (in cm)? "))
prop = 0 #float(input("What is the induction proportionality constant? "))
RE = float(35) #float(input("What is the value of the external resistance? "))
RM = float(2500) #float(input("What is the DC value of the membrane resistance? "))
RMin = float(1500) #float(input("What is the minimum membrance resistance? "))
RI = float(70) #float(input("What is the internal resistance? "))
IND = float(0.0021) #float(input("What is the value of inductance? "))
CDC = 1.1 #float(input("What is the DC capacitance (in microF)? "))
CINF = 0.55 # float(input("What is the capacitance at infinity? "))
alpha = float(10) #float(input("What is the decay rate?" ))

xf = 0.04
tf = 0.002

def constants1():
    re = RE/(np.pi*h*d)
    rm = RM/(np.pi*d)
    ri = 4*RI/(np.pi*d**2)
    g = 1/rm
    R = re+ri
    return(R, g, re, rm, ri)

def cap_mem(f):
    c_del = CDC - CINF
    a = CINF + c_del/(1+f*2e-6)
    return a*np.pi*d*10**(-6)

def cm(f):
    return CDC*np.pi*d*10**(-6)

def LC(f):
    return (1+prop*f)*IND*cm(f)

def LC1(f):
    return (1+prop*f)*IND*cap_mem(f)

Const = constants1()
R = Const[0]
G = Const[1]

def gam1(f):
    g = LC(f)*(-2*np.pi*f*1j - alpha)**2 + (R*cm(f)+(1+prop*f)*IND*G)*(-2*np.pi*f*1j-alpha) +R*G
    return cmath.sqrt(g)

def gam2(f):
    h = LC1(f)*(-2*np.pi*f*1j - alpha)**2 + (R*cap_mem(f)+(1+prop*f)*IND*G)*(-2*np.pi*f*1j-alpha) +R*G
    return cmath.sqrt(h)

def sol1(f,x,t,I,R):
    t1 = np.exp(-alpha*t - x*R)
    t2 = np.sin(2*np.pi*f*t + x*I)
    return t1*t2

def var_shunt(f):
    r_del = RM - RMin
    a = RMin + r_del/(1+f*5e-5)
    return np.pi*d/a

def constants4(f):
    tau = cap_mem(f)*var_res(f)
    re = RE/(np.pi*2*h*d)
    rm = RM/(np.pi*d)
    ri = 4*RI/(np.pi*d**2)
    lam = np.sqrt(rm/(re+ri))
    return (tau, lam, re, rm, ri) 

def gam4(f):
    h = LC1(f)*(-2*np.pi*f*1j - alpha)**2 + (R*cap_mem(f)+(1+prop*f)*IND*var_shunt(f))*(-2*np.pi*f*1j-alpha) 
    +R*var_shunt(f)
    return cmath.sqrt(h)

results_dic1 = {}
results_dic2 = {}
results_dic3 = {}
for i in range(3,8):
    f = 10**i
    t = tf 
    a = ['4(a)','4(b)','4(c)','4(d)','(f)']
    G1 = gam1(f)
    G2 = gam2(f)
    G4 = gam4(f)
    R1 = G1.real
    I1 = G1.imag
    R2 = G2.real
    I2 = G2.imag
    R4 = G4.real
    I4 = G4.imag
    def solution(x,t,w):
        S1 = sol1(w,x,t,I1,R1)
        return S1 
    def solution1(x,t,w):
        S2 = sol1(w,x,t,I2,R2)
        return S2
    def solution4(x,t,w):
        S4 = sol1(w,x,t,I4,R4)
        return S4
    xin1 = -t*2*np.pi*f/I1
    xin2 = -t*2*np.pi*f/I2
    xin4 = -t*2*np.pi*f/I4
    XSpace1 = np.linspace(0,xf,100001)
    XSpace2 = np.linspace(0,xf,100001)
    XSpace4 = np.linspace(0,xf,100001)
    II1 = abs(solution(XSpace1,t,f))
    III1 = II1.tolist()
    results_dic1[i-1] = III1
    II2 = solution1(XSpace2,t,f)
    III2 = II2.tolist()
    results_dic2[i-1] = III2
    II4 = solution4(XSpace4,t,f)
    III4 = II4.tolist()
    results_dic3[i-1] = III4
    lines = ["k-","k--","k:"]
    plt.plot(XSpace1,II1,lines[0], label = 'Telegrapher') 
    plt.axhline(0, color='black', linewidth=.5)
    plt.axvline(0, color='black', linewidth=.5)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.xlabel("Distance from source", fontsize = 12)  # add X-axis label
    plt.ylabel("Current", fontsize =12)  # add Y-axis label
    plt.title(a[i-3]+" Current at time "+str(tf)+"; frequency = "+str(f), fontsize =12)
    plt.legend(loc='best')
    plt.savefig(name+"Var_cap_Freq="+str(f)+"xf="+str(xf)+"_tf="+str(tf)+".png")
    plt.show()
