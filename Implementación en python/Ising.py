import numpy as np
from functools import reduce

# Definimos las matrices de Pauli y la identidad
I = np.eye(2) #Función numpy que aplica una matriz identidad 2x2
sx = np.array([[0, 1], [1, 0]])  #sigma_x
sz = np.array([[1, 0], [0, -1]]) #sigma_z

# Función para hacer el producto tensorial de una lista de matrices
def producto_tensorial(matrices):
    return reduce(np.kron, matrices) #np.kron aplica el producto tensorial delta de kronecker a una lista de matrices 
                                     #el "reduce" lo aplica iterativamente
    
# Construcción de la matriz hamiltoniana para N espines
def hamiltoniano(N, J, g):
    H = np.zeros((2**N, 2**N)) #Dimensiones

    for i in range(N - 1): #términos de sigmaz*sigmaz
        op_z = [I] * N # creamos una lista con N matrices identidad (uno por cada espín).
        op_z[i] = sz #operador sz
        op_z[i + 1] = sz #operador sz+1
        H -= J * producto_tensorial(op_z) #calcula el producto tensorial entre las sz

    # Término de campo en x: sigma_x en cada sitio
    for i in range(N):
        op_x = [I] * N
        op_x[i] = sx
        H -= g * producto_tensorial(op_x) 

    return H

#####################

N = 4      #espines
J = 1.0    #acoplamiento
g = 1.0    #componente del campo transversal
print(hamiltoniano(4, 1.0, 1.0))

######################

import matplotlib.pyplot as plt
from scipy.linalg import expm

# Estado inicial, todos los espines arriba 
def estado_inicial(N):
    arriba = np.array([[1], [0]])  # estado |0> (spin arriba)
    psi0 = producto_tensorial([arriba] * N)
    return psi0

# Operador σ^z en el sitio i
def sigma_z_en_i(N, i):
    ops = [I] * N
    ops[i] = sz
    return producto_tensorial(ops)

# Evolución temporal
def evolucion(psi0, H, tiempos):
    estados = []
    for t in tiempos:
        U = expm(-1j * H * t)
        psi_t = U @ psi0 #@ funciona para producto matricial
        estados.append(psi_t)
    return estados

# Observables: <σ^z_i>(t) para cada i
def magnetizacion_z(estados, N):
    magnetizaciones = np.zeros((len(estados), N))
    for t_idx, psi in enumerate(estados):
        for i in range(N):
            Sz_i = sigma_z_en_i(N, i)
            # Corrige la advertencia usando .item()
            magnetizaciones[t_idx, i] = np.real((psi.conj().T @ Sz_i @ psi).item())
    return magnetizaciones

### Gráfico

N = 4       # Número de espines
J = 1.0     # Acoplamiento
g = 1.0     # Campo transversal

# Estado inicial y evolución
psi0 = estado_inicial(N)
tiempos = np.linspace(0, 10, 200)  # tiempos de evolución
estados = evolucion(psi0, hamiltoniano(N, J, g), tiempos)

# Calcular magnetización en z para cada sitio
magz = magnetizacion_z(estados, N)

# Visualizar
plt.figure(figsize=(10, 6))
for i in range(N):
    plt.plot(tiempos, magz[:, i], label=f'Sitio {i}')
plt.xlabel('Tiempo')
plt.ylabel(r'$\langle \sigma^z_i \rangle$')
plt.title('Evolución temporal del modelo de Ising')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
