import os
import csv
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

def csv_list(s):
    with open(s) as file:
        F = csv.reader(file)
        L = []
        for f in F:
            L.append(f)
    return L


Datos = {} 
#Datos[med_n] = {R:medn_R,...,Frec:frecs}
for f in os.listdir():
    if f.endswith('.csv'):
        i = f.find('med')
        s = f[i+3:]
        j = s.find('_')
        n = s[:j]
        A = s[j:][1]
        medn = csv_list(f)[0]
        x0,dx,x = list(map(int,f[:i-1].split('_')))
        frecs = np.arange(x0,x+dx,dx)
        m = 'med_{}'.format(n)
        #Datos[m] = {'Frec': frecs}
        try:
            Datos[m][A] = medn
        except:
            Datos[m] = {A: medn}
        Datos[m]['Frec'] = frecs


##Correr esto para ver si las long estan bien
#for i in Datos.keys():
#    for j in Datos[i].keys():
#        if not len(Datos[i][j][0]) == len(Datos[i][j][1]):
#            print('Datos[\'{}\'][\'{}\'] anda mal'.format(i,j))

Datos['med_5']['Muestra'] = ('muestra maciza','cobre',{'radio':12.5,'error':0.5})
Datos['med_6']['Muestra'] = ('muestra maciza','cobre',{'radio':9.5,'error':0.5})
Datos['med_7']['Muestra'] = ('muestra maciza','alumninio',{'radio':11,'error':0.5})
Datos['med_9']['Muestra'] = ('muestra hueca','cobre',{'radio externo':16.44,'espesor':0.7,'error':0.02})
Datos['med_10']['Muestra'] = ('muestra hueca','cobre',{'radio externo':16.44,'espesor':0.7,'error':0.02})
Datos['med_11']['Muestra'] = ('muestra hueca','cobre',{'radio externo':13,'espesor':3,'error':0.02})
Datos['med_12']['Muestra'] = ('muestra hueca','cobre',{'radio externo':13,'espesor':3,'error':0.02,'Vpp':8})



'''
#Para reproducir lo que hicieron Andy y Nehuen
f10 = list(Datos['med_10']['Frec'])
f11 = list(Datos['med_11']['Frec'])
r10 = list(map(float,Datos['med_10']['R']))
r11 = list(map(float,Datos['med_11']['R']))
plt.plot(f10,r10,f11,r11)
plt.show()
'''

#mu es constante
mu = 1

#C es constante de proporcionalidad entre V y chi
C = 1/100000

#chi(rho,omega,a)
def chi(rho,omega,a):
    try:
        delta = np.sqrt((2*rho)/(mu*omega))
        k = (1+1j)/delta
        J0 = sp.jv(0,k*a)
        J1 = sp.jv(1,k*a)
        X = (2*J1)/(k*a*J0-1)
        return X
    except:
        return np.NaN

'''

m = 'med_7'
W = Datos[m]['Frec']
X_list = list(map(float,Datos[m]['X']))
Y_list = list(map(float,Datos[m]['Y']))
R_list = list(map(float,Datos[m]['R']))
T_list = list(map(float,Datos[m]['T']))
x_list = [r*np.cos(t) for r,t in zip(R_list,T_list)]
y_list = [r*np.sin(t) for r,t in zip(R_list,T_list)]

plt.subplot(2,1,1)
plt.plot(W,X_list,label='x(X)')
plt.plot(W,x_list,label='x(R,T)')
plt.legend()
plt.subplot(2,1,2)
plt.plot(W,Y_list,label='y(Y)')
plt.plot(W,y_list,label='y(R,T)')
plt.legend()
plt.show()

'''

def difer(x):
    x = [t for t in x if t!=10 and type(t)!=str]
    try:
        return max(max(np.abs(a-b) for a in x) for b in x)
    except:
        return 0

#Solo algunas med tienen "a"
OMRHO = {}
for m in ['med_5','med_6','med_7']:
    a = Datos[m]['Muestra'][2]['radio']
    N = len(Datos[m]['Frec'])
    RHO_RT,RHO_XY,OMEGA = [],[],[]
    for i in range(N):
        omega = Datos[m]['Frec'][i]
        r = float(Datos[m]['R'][i])
        t = float(Datos[m]['T'][i])
        x = r*np.cos(t)
        y = r*np.sin(t)
        chi_Re = -y/(omega*C)
        chi_Im = x/(omega*C)
        Chi_Re_rho = lambda rho: chi(rho,omega,a).real -chi_Re
        rho_guess_re = 10.0
        rho_r = fsolve(Chi_Re_rho, rho_guess_re)
        Chi_Im_rho = lambda rho: chi(rho,omega,a).imag -chi_Im
        rho_guess_im = 10.0
        rho_i = fsolve(Chi_Im_rho, rho_guess_im)
        RHO_RT.append((rho_r,rho_i))
        x = float(Datos[m]['X'][i])
        y = float(Datos[m]['Y'][i])
        chi_Re = -y/(omega*C)
        chi_Im = x/(omega*C)
        Chi_Re_rho = lambda rho: chi(rho,omega,a).real -chi_Re
        rho_guess_re = 10.0
        rho_r = fsolve(Chi_Re_rho, rho_guess_re)
        Chi_Im_rho = lambda rho: chi(rho,omega,a).imag -chi_Im
        rho_guess_im = 10.0
        rho_i = fsolve(Chi_Im_rho, rho_guess_im)
        RHO_XY.append((rho_r,rho_i))
        OMEGA.append(omega)
    RHO = [('omega: {}'.format(om),float(rt[0]),float(rt[1]),float(xy[0]),float(xy[1]))
           for om,rt,xy in zip(OMEGA,RHO_RT,RHO_XY)]
    RHO.sort(key=difer)
    #for r in RHO: print(r)
    om_rho = {}
    for x in RHO:
        if x.count(10) not in [3,4]:
            if difer(x) < np.inf:
                try:
                    om_rho['omega'].append(float(x[0][7:]))
                    om_rho['rho'].append(np.mean([t for t in x if t!=10 and type(t)!=str]))
                except:
                    om_rho['omega'] = [float(x[0][7:])]
                    om_rho['rho'] = [np.mean([t for t in x if t!=10 and type(t)!=str])]
    om_rho = list(zip(om_rho['omega'],om_rho['rho']))
    om_rho.sort(key=lambda x:x[0])
    om_plot = [x[0] for x in om_rho]
    rho_plot = [x[1] for x in om_rho]
    OMRHO[m] = {'omega':om_plot,'rho':rho_plot}

D = OMRHO['med_7']
#plt.plot(D['omega'],D['rho'],'*',label='Rho')
Delta = [np.sqrt((2*rho)/(mu*omega)) for rho in D['rho']]
plt.plot(D['omega'],Delta,'*',label='Delta')
plt.legend()
plt.grid()
plt.show()