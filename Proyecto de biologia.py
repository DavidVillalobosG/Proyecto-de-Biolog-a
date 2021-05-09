import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp

#Resuelve el sistema de ecuaciones diferenciales
def EcuacionDiferencial (y,t):
    global Cm, gk, gNa, gl, Vk, VNa, Vl
    dy = np.zeros((4,)) #Crea una lista de ceros, donde cada entrada representará una ecuación diferencial
    Vm = y[0]
    h = y[1]
    m = y[2]
    n = y[3]

    dy[0] = Estimulo(t)/Cm - gk*n**4*(Vm-Vk)/Cm - gNa*m**3*h*(Vm-VNa)/Cm - gl*(Vm-Vl)/Cm
    dy[1] = alfaH(Vm)*(1-h)-betaH(Vm)*h
    dy[2] = alfaM(Vm)*(1-m)-betaM(Vm)*m
    dy[3] = alfaN(Vm)*(1-n)-betaN(Vm)*n

    return dy

#Esta función simula un estímulo de 200mV por un 1ms cada 5 ms.
def Estimulo (t):
    if 0<t<1:
        return 200
    if 5<t<6:
        return 200
    if 10<t<11:
        return 200
    if 15<t<16:
        return 200
    else:
        return 0


def alfaH (Vm):
    return 0.07*np.exp(-Vm/20)

def betaH (Vm):
    return 1/(np.exp((30-Vm)/10)+1)

def alfaM (Vm):
    return (25-Vm)/(10*(np.exp((25-Vm)/10)-1))

def betaM (Vm):
    return 4*np.exp(-Vm/18)

def alfaN (Vm):
    return (10-Vm)/(100*(np.exp(10-Vm)-1))

def betaN (Vm):
    return 0.125*np.exp(-Vm/80)

#Espacio de tiempo a trabajar

t = np.linspace(0,20,100000)

#Parámetros para calcular las EDO
Cm = 1
gk = 35 #(mS/cm^2)
Vk = -77 #(mV)
gNa = 40 #(mS/cm^2)
VNa = 55 #(mV)
gl = 0.3 #(mS/cm^2)
Vl = -65 #(mS/cm^2)

#Se establecen las condiciones para el estado estacionario

V0 = 0
h0 = alfaH(0)/(alfaH(0)+betaH(0))
m0 = alfaM(0)/(alfaM(0)+betaM(0))
n0 = alfaN(0)/(alfaN(0)+betaH(0))

#Se establece la condición inicial según el estado estacionario

y = np.array([V0,h0,m0,n0])

#Se soluciona la EDO

Resultado = scp.odeint(EcuacionDiferencial,y,t)

#Se grafican los resultados

fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(t, Resultado[:, 0])
ax.set_xlabel('Tiempo (ms)')
ax.set_ylabel('Vm (mV)')
ax.set_title('Potencial en la Neurona')
plt.grid()

plt.show()

#http://www.red-mat.unam.mx/foro/volumenes/vol024/TesisRicardalast-f.pdf
#https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/