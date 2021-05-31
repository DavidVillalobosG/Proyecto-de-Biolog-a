# Importación de librerías

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp



# Definición de parámetros globales

Cm = 1  #(microF/cm^2)
gk = 36  #(mS/cm^2)
Vk = -12  #(mV)
gNa = 120  #(mS/cm^2)
VNa = 115  #(mV)
gl = 0.3  #(mS/cm^2)
Vl = 10.613  #(mV)



# Definición de funciones

def EcuacionDiferencial (y,t,i):
    """
    Establece el sistema de ecuaciones diferenciales a resolver
    :param y: Arreglo con los valores de Vm, h, m y n
    :param t: Arreglo con los valores de t
    :param i: Amplitud de la corriente de estímulo
    :return: Lista con las ecuaciones diferenciales evaluadas
    """

    dy = np.zeros((4,)) # Crea una lista de ceros, donde cada entrada representará una ecuación diferencial
    Vm = y[0]
    h = y[1]
    m = y[2]
    n = y[3]

    dy[0] = Estimulo(t, i)/Cm - gk*n**4*(Vm-Vk)/Cm - gNa*m**3*h*(Vm-VNa)/Cm - gl*(Vm-Vl)/Cm  # Ecuación 13, dV/dt
    dy[1] = alfaH(Vm)*(1-h)-betaH(Vm)*h  # Ecuación 10, dh/dt
    dy[2] = alfaM(Vm)*(1-m)-betaM(Vm)*m  # Ecuación 7, dm/dt
    dy[3] = alfaN(Vm)*(1-n)-betaN(Vm)*n  # Ecuación 4, dn/dt

    return dy


# Función estímulo

def Estimulo (t,i=150): # Se agrega t como parámetro para poder generar estímulos dependientes del tiempo
    """
    Esta función genera un estímulo constante de amplitud i (microA/cm^2)
    :param t: Float, tiempo t a evaluar
    :param i: Float, amplitud de la corriente de estímulo
    :return: Corriente de estímulo para el t dado
    """
    # Se muestra en forma de comentario un ejemplo de aplicación de estímulos dependientes del tiempo
    # if 5<t<6:
        # return i
    # if 15<t<16:
        # return i
    # else:
        # return 0

    # Sin embargo, para la implementación se utiliza una corriente constante en el tiempo
    return i


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
    return (10-Vm)/(100*(np.exp((10-Vm)/10)-1))  # Ecuación 5

def betaN (Vm):
    return 0.125*np.exp(-Vm/80)  # Ecuación 6


def CalcularFrecuencia(t,V):
    """
    Calcula la frecuencia de una función utilizando los primeros dos máximos locales
    :param t: Arreglo con los valores de t (eje x)
    :param resultado: Arreglo con los valores de V (eje y)
    :return: Float, frecuencia de la función
    """
    # Se determina la posición del primer máximo local
    max1 = -70
    for i in range(len(t) - 1):
        if V[i] > max1:
            max1 = V[i]
            indiceMax1 = i

    # Se determina la posición del mínimo local más próximo
    min = max1
    for i in range(indiceMax1, len(t) - 1):
        if V[i] < min:
            min = V[i]
            indiceMin = i

    # Se determina la posición del segundo máximo local
    max2 = -70
    for i in range(indiceMin, len(t) - 1):
        if V[i] > max2:
            max2 = V[i]
            indiceMax2 = i

    # Se calcula el período como la diferencia en t entre los dos máximos locales encontrados
    periodo = t[indiceMax2] - t[indiceMax1]
    frecuencia = 1000/periodo # Hz

    return frecuencia


# Evalúa los valores de corriente en el modelo teórico
def RespuestaFrecuenciaTeorica(I):
    cita = 27*np.log(I+1)
    return cita

# Evalua los valores de corriente en el modelo de ajuste
def RespuestaFrecuenciaAjuste(I):
    citaAj = 26.904*np.log(I)+1.3729
    return citaAj


# Main

# Espacio de tiempo a trabajar, en ms
t = np.linspace(0,100,1000)

# Estado inicial
V0 = 0
h0 = alfaH(V0)/(alfaH(V0)+betaH(V0))
m0 = alfaM(V0)/(alfaM(V0)+betaM(V0))
n0 = alfaN(V0)/(alfaN(V0)+betaN(V0))

# Se almacenan las condiciones iniciales en un arreglo
y0 = np.array([V0,h0,m0,n0])

# Se establecen amplitudes para la corriente de estímulo
i1 = 2
i2 = 8
i3 = 14
i4 = 20

# Se soluciona la EDO
# Se resta -70 para graficar la tensión absoluta en vez de la diferencia con respecto al potencial de reposo
resultado1 = scp.odeint(EcuacionDiferencial, y0, t, args=(i1,)) - 70
resultado2 = scp.odeint(EcuacionDiferencial, y0, t, args=(i2,)) - 70
resultado3 = scp.odeint(EcuacionDiferencial, y0, t, args=(i3,)) - 70
resultado4 = scp.odeint(EcuacionDiferencial, y0, t, args=(i4,)) - 70



# Determinación de la relación entre frecuencia de respuesta y amplitud de la corriente de estímulo

# Se genera una lista con amplitudes de corriente a evaluar
listaAmplitudes = np.linspace(15,400,800)
# Se inicializa una lista para almacenar las respuestas obtenidas para cada amplitud de corriente
listaRespuestas = []
# Se inicializan listas para almacenar las frecuencias para cada respuesta
listaFrecuencias = []
listaFrecuenciasTeorica = []
listaFrecuenciasAjuste = []

# Se calcula y almacena la frecuencia de respuesta para cada amplitud de corriente
for i in listaAmplitudes:
    resultado = scp.odeint(EcuacionDiferencial, y0, t, args=(i,)) - 70
    listaRespuestas.append(resultado[:, 0])
for irespuesta in listaRespuestas:
    listaFrecuencias.append(CalcularFrecuencia(t, irespuesta))

# Se calcula y almacena la frecuencia teórica y la frecuencia obtenida con el modelo de ajuste
for i in range(0, len(listaAmplitudes)):
    teoria = RespuestaFrecuenciaTeorica(listaAmplitudes[i])
    ajuste = RespuestaFrecuenciaAjuste(listaAmplitudes[i])
    listaFrecuenciasTeorica.append(teoria)
    listaFrecuenciasAjuste.append(ajuste)


# Se grafican los resultados

# Potencial en la neurona para diversas amplitudes de corriente de estímulo
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(t, resultado1[:, 0])
ax.plot(t, resultado2[:, 0])
ax.plot(t, resultado3[:, 0])
ax.plot(t, resultado4[:, 0])
ax.set_xlabel('Tiempo (ms)')
ax.set_ylabel('Vm (mV)')
ax.set_title('Potencial en la Neurona')
plt.legend(['2 microA/cm^2', '8 microA/cm^2', '14 microA/cm^2', '20 microA/cm^2'])
plt.grid()
plt.show()

# Frecuencia de respuesta en función de la amplitud de la corriente de estímulo, obtenida por simulación
fig, ax = plt.subplots(figsize=(12, 7))
ax.plot(listaAmplitudes, listaFrecuencias)
ax.set_xlabel('Amplitud de la corriente de estímulo (microA/cm^2)')
ax.set_ylabel('Frecuencia de respuesta (Hz)')
ax.set_title('Relación entre la amplitud del estímulo y la frecuencia de respuesta')
plt.grid()
plt.show()

# Comparación entre la relación teórica y la modelada
fig, ax = plt.subplots(figsize = (12,7))
ax.plot(listaAmplitudes, listaFrecuenciasTeorica)
ax.plot(listaAmplitudes, listaFrecuenciasAjuste)
ax.set_xlabel('Amplitud de la corriente de estímulo (microA/cm^2)')
ax.set_ylabel('Frecuencia de respuesta (Hz)')
ax.set_title('Comparación entre modelo de ajuste propuesto y el modelo teórico')
plt.legend(['Modelo teórico','Modelo de ajuste'])
plt.grid()
plt.show()