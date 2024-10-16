import numpy as np
import scipy.integrate as integrate
from matplotlib import pyplot as plt
cimport numpy as np
import csv
from multiprocessing import Pool
from multiprocessing import Process
from concurrent.futures import ProcessPoolExecutor

#n = int
#RE = float
#RM = float
#RMin = float
#RI = float
#CM = float 
#alpha = float

cdef np.float_t d = 0.0005
cdef np.float_t h = 0.0005

#def all_inputs():
#    global n
#    n = int(6)#int(input("How many frequencies? "))
    #Name = input("What is the designate for the graphs? ")
    #FR = float(input("What is the maximum frequncy (as a power of ten)? "))
#    global RE
#    RE = float(0)#float(input("What is the external resistance (Ohm-centimetre)? "))
#    global RM
#    RM = float(2500)#float(input("What is the DC resistance per unit area across membrane? (Ohm-centimetre squared.) "))
#    global RMin
#    RMin = float(1500) #float(input("What is the minimum resistance of the membrane? "))
#    global RI
#    RI = float(35) #float(input("What is the internal resistivilty? (Ohm-centimetre.)"))
#    global CM
#    CM = float(0.1)#float(input("What is the DC capacitance per unit area of membrane? (micro F-centimetre squared.) "))
#    global alpha
#    alpha =float(10) #float(input("What is the decay constant (alpha)? "))

cdef np.float_t xf = 0.04
cdef np.float_t tf = 0.00001

def constants1(CM1,RM1,RE1,RI1):
    cdef double tau = CM1*RM1/(10**6)
    cdef double re = RE1/(np.pi*h*d)
    cdef double rm = RM1/(np.pi*d)
    cdef double ri = 4*RI1/(np.pi*d**2)
    cdef double lam = np.sqrt(rm/(re+ri))
    return(tau, lam, re, rm, ri)

# xf and tf have to stay the same throughout 

def XTerm(X,  gam, omega):
    first = 0.5*X*np.exp(-X**2/(4*gam)-gam)
    second = np.sqrt(np.pi)*gam**(1.5)
    return -first/second

def sinterm(T, Tau, gam, omega, alpha):
    first = np.exp(-alpha*Tau*(T-gam))
    second = np.sin(2*np.pi*omega*Tau*(T-gam))
    return first*second

def current(X,T,Tau,gam,omeg,alpha):
    return XTerm(X,gam,omeg)*sinterm(T,Tau,gam,omeg,alpha)           

def int1(n,CM0,RM0,RE0,RI0,alpha0):
    #cdef double alpha = alpha0
    alpha = alpha0
    freq = [10**i for i in range(1,n+1)]
    omega = [2*np.pi*f for f in freq]
    #print(omega)
    #print(Tf, Xf)
    #TSpace = np.linspace(0,Tf,11)
    results_dic1 = {}
    for j in range(0,len(freq)):
        Const1 = constants1(CM0,RM0,RE0,RI0)
        Tau = Const1[0]
        Lam = Const1[1]
        print("Tau1 = ", Tau)
        print("Lam1 = ", Lam)
        Tf = tf/Tau
        Xf = xf/Lam
        XSpace = np.linspace(0,Xf,1033)
        w = freq[j]
        results_dic1[j] = []
        #print(omeg)
        T = Tf
        X = Xf
        for X in XSpace:
            Sum = 0
            k = 1
            def cur(gamma):
                return current(X,T,Tau,gamma,w,alpha)
            while k <= 2*w*Tau*T -1:
                #T = TSpace[-1]
                gamtimes = np.linspace(k/(2*w*Tau),(k+1)/(2*w*Tau),8)
                cur_array = cur(gamtimes)
                I = np.trapz(cur_array,gamtimes)
                #I = integrate.quad(cur, 0.00001, T)
                Sum = Sum + I
                #y = np.array(results_dic[t])
                k += 1
            gamtimes1 = np.linspace(1/(2*w*Tau*10**5),1/(2*w*Tau),10**5)
            cur_array1 = cur(gamtimes1)
            I1 = np.trapz(cur_array1,gamtimes1)
            Sum = Sum + I1
            results_dic1[j].append(Sum)
    return results_dic1

# Now for variable capacitance



def cap_mem(f,CM2, CINF2):
    cdef double c_inf = CINF2
    cdef double c_del = CM2 - c_inf
    cdef double a = c_inf + c_del/(1+f*2e-3)
    return a*10**(-6)

def constants2(f,CM2,CINF2,RM2,RE2,RI2):
    cdef double tau = cap_mem(f,CM2,CINF2)*RM2
    cdef double re = RE2/(np.pi*h*d)
    cdef double rm = RM2/(np.pi*d)
    cdef double ri = 4*RI2/(np.pi*d**2)
    cdef double lam = np.sqrt(rm/(re+ri))
    return (tau, lam, re, rm, ri) 
    
def XTerm2(X,Tau,gam,omega):
    first = 0.5*X*np.exp(-X**2/(4*gam)-gam)
    second = np.sqrt(np.pi)*gam**(1.5)
    return first/second

def sinterm2(T,Tau,gam,omega,alpha):
    first = np.exp(-alpha*Tau*(T-gam))
    second = np.sin(2*np.pi*omega*Tau*(T-gam))
    return first*second

def current2(X,T,Tau,gam,omeg,alpha):
    return XTerm2(X,Tau, gam, omeg)*sinterm2(T, Tau, gam,omeg,alpha)           

def int2(n,CM0,CINF2,RM0,RE0,RI0,alpha0):
    freq = [10**i for i in range(1,n+1)]
    omega = [2*np.pi*f for f in freq]
    results_dic2 = {}
    cdef double alpha = alpha0
    for j in range(len(freq)):
        w = freq[j]
        #CM = cap_mem(w)
        #print(CM)
        Const = constants2(w,CM0,CINF2,RM0,RE0,RI0)
        #print(Const)
        Tau = Const[0]
        print("Tau2 = ", Tau)
        Lam = Const[1]
        print("Lam2 = ", Lam)
        Tf = tf/Tau
        Xf = xf/Lam
        #print(Tf,Xf)
        #TSpace = np.linspace(0,Tf,11)
        XSpace = np.linspace(0,Xf,1033)
        results_dic2[j] = []
        T = Tf
        X = Xf
        for X in XSpace:
            Sum = 0
            k = 1
            def cur2(gamma):
                return current2(X,T,Tau,gamma,w,alpha)
            while k <= 2*w*Tau*T - 1:
                gamtimes2 = np.linspace(k/(2*w*Tau),(k+1)/(2*w*Tau),8)
                cur_array2 = cur2(gamtimes2)
                I = np.trapz(cur_array2,gamtimes2)
                #I = integrate.quad(cur, 0.0001, T)
                Sum = Sum + I
                k += 1
            gamtimes2 = np.linspace(1/(2*w*Tau*10**5),1/(2*w*Tau),100000)
            cur_array2 = cur2(gamtimes2)
            I2 = np.trapz(cur_array2,gamtimes2)
            Sum = Sum + I2
            results_dic2[j].append(Sum)
    return results_dic2

# Now for membrane conductivity, added to variable capacitance


#def constants3(f):
    #print("Note that this function takes resistivity as arguments.")
    #r_del = RE - RMin
#    tau = cap_mem(f)*var_res(f)
#    re = RE/(np.pi*h*d)
#    rm = var_res(f)/(np.pi*d)
#    ri = 4*RI/(np.pi*d**2)
#    lam = np.sqrt(rm/(re+ri))
#    return (tau, lam, re, rm, ri) 
    
#def XTerm3(X, Tau, gam, omega):
#    first = 0.5*X*np.exp(-X**2/(4*gam)-gam)
#    second = np.sqrt(np.pi*gam**(1.5))
#    return -first/second

#def sinterm3(T, Tau, gam, omega):
#    first = np.exp(-alpha*Tau*(T-gam))
#    second = np.sin(2*np.pi*omega*Tau*(T-gam))
#    return first*second

#def current3(X,T,Tau, gam,omeg):
#    return XTerm2(X,Tau, gam, omeg)*sinterm2(T, Tau, gam,omeg)           

    
#def int3(): 
#    freq = [10**i for i in range(1,n+1)]
#    omega = [2*np.pi*f for f in freq]
#    results_dic3 = {}
#    for j in range(len(omega)):
#        w = omega[j]
#        CM = cap_mem(w)
        #print(CM)
#        Const = constants3(omega[j])
#        print("Tau3 = ",Const[0])
        #print(Const)
#        Tau = Const[0]
        #print("Tau = ", Tau)
#        Lam = Const[1]
#        print("Lam3 = ", Const[1] )
#        Tf = tf/Tau
#        Xf = xf/Lam
        #print(Tf,Xf)
        #TSpace = np.linspace(0,Tf,11)
#        XSpace = np.linspace(0,Xf,101)
#        results_dic3[j] = []
#        T = Tf
#        for X in XSpace:
#            Sum = 0
#            k = 1
#            while k <= 2*w*Tau*T - 1: 
#                gamtimes3 = np.linspace(k/(2*w*Tau),(k+1)/(2*w*Tau),4)
#                def cur3(gamma):
#                    return current3(X,T,Tau,gamma,w)
#                cur_array3 = cur3(gamtimes3)
#                I = np.trapz(cur_array3,gamtimes3)
                #I = integrate.quad(cur, 0.0001, T)
#               Sum = Sum + I
#               k += 1
#            results_dic3[j].append(Sum)
#    return results_dic3

# Now add increased conductance of extracellular medium, together with all the others

#Ext_Res = 15
#Min_Res = 20
#Res_del = Ext_Res - Min_Res

#def var_ext_res(f):
#    a = Ext_Res + Min_Res*1000/(1000+f)
#    return a

def var_res(f,RM,RMin):
    cdef double r_del = RM - RMin
    cdef double a = RMin + r_del/(1+f*2e-4)
    return a

def constants4(f,CM3,CINF2,RM3,RMin3,RE3,RI3):
    #print("Note that this function takes resistivity as arguments.")
    cdef double tau = cap_mem(f,CM3,CINF2)*var_res(f,RM3,RMin3)
    cdef double re = RE3/(np.pi*2*h*d)
    cdef double rm = var_res(f,RM3,RMin3)/(np.pi*d)
    #rm = var_res(f)/(np.pi*d)
    cdef double ri = 4*RI3/(np.pi*d**2)
    cdef double lam = np.sqrt(rm/(re+ri))
    return (tau, lam, re, rm, ri) 
    
def XTerm4(X,Tau,gam,omega):
    first = 0.5*X*np.exp(-X**2/(4*gam)-gam)
    second = np.sqrt(np.pi)*gam**(1.5)
    return first/second

def sinterm4(T,Tau,gam,omega,alpha):
    first = np.exp(-alpha*Tau*(T-gam))
    second = np.sin(2*np.pi*omega*Tau*(T-gam))
    return first*second

def current4(X,T,Tau,gam,omeg,alpha):
    return XTerm2(X,Tau, gam, omeg)*sinterm2(T, Tau, gam,omeg,alpha)           

    
def int4(n,CM0,CINF2,RM0,RMin,RE0,RI0,alpha0): 
    freq = [10**i for i in range(1,n+1)]
    omega = [2*np.pi*f for f in freq]
    results_dic4 = {}
    cdef double alpha = alpha0
    for j in range(len(freq)):
        w = freq[j]
        #CM = cap_mem(w,CM0,CINF2)
        #print(CM)
        Const = constants4(w,CM0,CINF2,RM0,RMin,RE0,RI0)
        print("Tau4 = ", Const[0])
        #print(Const)
        Tau = Const[0]
        Lam = Const[1]
        print("Lam4 = ", Lam)
        Tf = tf/Tau
        Xf = xf/Lam
        #print(Tf,Xf)
        XSpace = np.linspace(0,Tf,1033)
        #XSpace = np.linspace(0,Xf,101)
        results_dic4[j] = []
        T = Tf
        X = Xf
        for X in XSpace:
            Sum = 0
            k = 1
            def cur4(gamma):
                return current4(X,T,Tau,gamma,w,alpha)
            while k < 2*w*Tau*T - 1:
                gamtimes4 = np.linspace(k/(2*w*Tau),(k+1)/(2*w*Tau),8)
                cur_array4 = cur4(gamtimes4)
                I = np.trapz(cur_array4,gamtimes4)
                #I = integrate.quad(cur, 0.0001, T)
                Sum = Sum + I
                k += 1
            gamtimes4 = np.linspace(1/(2*w*Tau*10**5),1/(2*w*Tau),100000)
            cur_array4 = cur4(gamtimes4)
            I4 = np.trapz(cur_array4,gamtimes4)
            Sum = Sum + I4
            results_dic4[j].append(Sum)
    return results_dic4

#AI = all_inputs()
