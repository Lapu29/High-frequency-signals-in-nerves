import csv
import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt
import numpy as np
import ast
import sys

csv.field_size_limit(sys.maxsize)

d = 0.0005
h = 0.0005

xf = 0.04
tf = 0.002

n = int(6) #int(input("How many frequencies?"))
Name = 'No_ind_BW_' #input("What is the designate for the graphs?")
#FR = float(input("What is the maximum frequncy (as a power of ten)? "))
RE = float(35) #float(input("What is the external resistance (Ohm-centimetre)? "))
RM = float(2500) #float(input("What is the DC resistance per unit area across membrane? (Ohm-centimetre squared.) "))
RMin = float(1500) #float(input("What is the minimum resistance of the membrane? "))
RI = float(70) #float(input("What is the internal resistivity? (Ohm-centimetre.)"))
CM = float(1.1) #float(input("What is the DC capacitance per unit area of membrane? (micro F-centimetre squared.) "))
CMin = float(0.55) #float(input("What is the minimum capacitance of the membrane? "))

#reader = csv.DictReader(open('myfile.csv'))

with open('Cytrial/Test_Int1_xf=0.04_tf=0.002.csv', 'r') as ff:
    reader1 = csv.reader(ff)
    Z1 = list(reader1)

Y1 = ast.literal_eval(Z1[0][0])    

    
with open('Cytrial/Test_Int2_xf=0.04_tf=0.002.csv', newline='') as gg:
    reader2 = csv.reader(gg)
    Z2 = list(reader2)

Y2 = ast.literal_eval(Z2[0][0])       
    
with open('Cytrial/Test_Int4_xf=0.04_tf=0.002.csv', newline='') as hh:
    reader3 = csv.reader(hh)
    Z4 = list(reader3)

Y4 = ast.literal_eval(Z4[0][0])    

def cap_mem(f,CM2,CINF2):
    c_inf = CINF2
    c_del = CM2 - c_inf
    a = c_inf + c_del/(1+2*np.pi*f*2e-3)
    return a*10**(-6)

def constants1(CM1,CINF,RM1,RE1,RI1,f):
    tau = cap_mem(f,CM1,CINF)*RM1/(10**6)
    re = RE1/(np.pi*h*d)
    rm = RM1/(np.pi*d)
    ri = 4*RI/(np.pi*d**2)
    lam = np.sqrt(rm/ri)
    return(tau, lam, re, rm, ri)

def constants2(f):
    tau = cap_mem(f)*RM
    re = RE/(np.pi*h*d)
    rm = RM/(np.pi*d)
    ri = 4*RI/(np.pi*d**2)
    lam = np.sqrt(rm/(re+ri))
    return (tau, lam, re, rm, ri) 

def var_res(f):
    r_del = RM - RMin
    a = RMin + r_del/(1+f*5e-5)
    return a

def constants4(f):
    tau = cap_mem(f)*var_res(f)
    re = RE/(np.pi*2*h*d)
    rm = RM/(np.pi*d)
    #rm = var_res(f)/(np.pi*d)
    ri = 4*RI/(np.pi*d**2)
    lam = np.sqrt(rm/(re+ri))
    return (tau, lam, re, rm, ri) 

for i in range(2,n):
    a = ['3(a)','3(b)','3(c)','3(d)','(e)']
    f = 10**(i+1)
    C1 = constants1(CM,CMin,RM,RE,RI,f)
    Lam = C1[1]
    Xf1 = xf/Lam
    XSpace1 = np.linspace(0,Xf1,1033)
    xspace1 = [x*Lam for x in XSpace1]
    C4 = constants1(CM,CMin,RM,RE,RI,f)
    Lam2 = C4[1]
    Xf2 = xf/Lam2
    XSpace2 = np.linspace(0,Xf2,1033)
    xspace2 = [x*Lam2 for x in XSpace2]
    I1 = np.array([abs(x) for x in Y1[i]])
    I2 = np.array([abs(x) for x in Y2[i]])
    I4 = np.array([abs(x) for x in Y4[i]])
    lines = ["k-","k--","k-.","k:"]
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.plot(xspace1,I1,lines[0], label = 'Constant capacitance')
    plt.plot(xspace1,I2,lines[1], label = 'Variable capacitance')
    plt.plot(xspace1,I4,lines[3], label = 'Variable membrane resistance')
    plt.xlabel("Distance from source", fontsize = 12)  # add X-axis label
    plt.ylabel("Current", fontsize = 12)  # add Y-axis label
    plt.title(a[i-2]+" Current at time "+str(tf)+"; frequency = "+str(f), fontsize = 12)
    plt.legend(loc='best')
    plt.savefig(Name+"Frequency="+str(f)+"xf="+str(xf)+"_tf="+str(tf)+".png")
    plt.show()