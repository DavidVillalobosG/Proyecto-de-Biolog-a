import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp


def EcuacionDiferencial (y,t,i):
    """
    Estable el sistema de ecuaciones diferenciales a resolver
    :param y: Arreglo con los valores de Vm, h, m y n
    :param t: Arreglo con los valores de t
    :param i: Amplitud de la corriente de estímulo
    :return: Lista con las ecuaciones diferenciales evaluadas
    """

    # Se definen los parámetros globales según la literatura
    global Cm, gk, gNa, gl, Vk, VNa, Vl
    dy = np.zeros((4,)) # Crea una lista de ceros, donde cada entrada representará una ecuación diferencial
    Vm = y[0]
    h = y[1]
    m = y[2]
    n = y[3]

    dy[0] = Estimulo(t,i)/Cm - gk*n**4*(Vm-Vk)/Cm - gNa*m**3*h*(Vm-VNa)/Cm - gl*(Vm-Vl)/Cm  # Ecuación 13, dV/dt
    dy[1] = alfaH(Vm)*(1-h)-betaH(Vm)*h  # Ecuación 10, dh/dt
    dy[2] = alfaM(Vm)*(1-m)-betaM(Vm)*m  # Ecuación 7, dm/dt
    dy[3] = alfaN(Vm)*(1-n)-betaN(Vm)*n  # Ecuación 4, dn/dt

    return dy


#Esta función simula un estímulo de corriente i microA/cm^2 por un 1ms en t=5ms y t=10ms
#La amplitud i tiene un valor default de 150microA/cm^2
def Estimulo (t,i=150):
    if 5<t<6:
        return i
    if 15<t<16:
        return i
    else:
        return 0


# Definición de variables dependientes del potencial de la membrana

def alfaH (Vm):
    return 0.07*np.exp(-Vm/20)  # Ecuación 11

def betaH (Vm):
    return 1/(np.exp((30-Vm)/10)+1)  # Ecuación 12

def alfaM (Vm):
    return (25-Vm)/(10*(np.exp((25-Vm)/10)-1))  # Ecuación 8

def betaM (Vm):
    return 4*np.exp(-Vm/18)  # Ecuación 9

def alfaN (Vm):
    return (10-Vm)/(100*(np.exp(10-Vm)-1))  # Ecuación 5

def betaN (Vm):
    return 0.125*np.exp(-Vm/80)  # Ecuación 6

#Espacio de tiempo a trabajar
t = np.linspace(0,20,100000)

#Parámetros para calcular las EDO
Cm = 1  #(microF/cm^2)
gk = 35  #(mS/cm^2)
Vk = -77  #(mV)
gNa = 40  #(mS/cm^2)
VNa = 55  #(mV)
gl = 0.3  #(mS/cm^2)
Vl = -65  #(mS/cm^2)

#Se establecen las condiciones para el estado estacionario

V0 = -70  # Potencial de reposo
h0 = alfaH(V0)/(alfaH(V0)+betaH(V0))
m0 = alfaM(V0)/(alfaM(V0)+betaM(V0))
n0 = alfaN(V0)/(alfaN(V0)+betaH(V0))

#Se almacenan las condiciones iniciales, según el estado estacionario, en un arreglo
y0 = np.array([V0,h0,m0,n0])

#Se establecen amplitudes para la corriente de estímulo
i1 = 25
i2 = 50
i3 = 150

#Se soluciona la EDO
resultado1 = scp.odeint(EcuacionDiferencial, y0, t, args=(i1,))
resultado2 = scp.odeint(EcuacionDiferencial, y0, t, args=(i2,))
resultado3 = scp.odeint(EcuacionDiferencial, y0, t, args=(i3,))

#Se grafican los resultados
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(t, resultado1[:, 0])
ax.plot(t, resultado2[:, 0])
ax.plot(t, resultado3[:, 0])
ax.set_xlabel('Tiempo (ms)')
ax.set_ylabel('Vm (mV)')
ax.set_title('Potencial en la Neurona')
plt.legend(['i1', 'i2', 'i3'])
plt.grid()
plt.show()

#http://www.red-mat.unam.mx/foro/volumenes/vol024/TesisRicardalast-f.pdf
#https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1392413/