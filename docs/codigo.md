# Código del modelo de Ising cuántico en Python

Este apartado presenta la implementación completa del **modelo de Ising cuántico unidimensional con campo transversal**, un sistema fundamental en la física de muchos cuerpos. A continuación, se explica cada parte del código que simula su dinámica cuántica.

---

## Hamiltoniano del modelo

El Hamiltoniano está dado por:

$$
H = -J \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - g \sum_{i=1}^N \sigma^x_i
$$

donde:
- \( \sigma^x \) y \( \sigma^z \) son matrices de Pauli,  
- \( J \) es la constante de acoplamiento entre espines vecinos,  
- \( g \) es la intensidad del campo transversal,  
- \( N \) es el número total de espines en la cadena.

---

## Implementación del modelo en Python

```python
import numpy as np
from functools import reduce
from scipy.linalg import expm
```

### Matrices de Pauli y matriz identidad

```python
I = np.eye(2)
sx = np.array([[0, 1], [1, 0]])   # σ_x
sz = np.array([[1, 0], [0, -1]])  # σ_z
```

### Producto tensorial

```python
def producto_tensorial(matrices):
    return reduce(np.kron, matrices)
```

### Construcción del Hamiltoniano

```python
def hamiltoniano(N, J, g):
    H = np.zeros((2**N, 2**N))
    
    # Término de interacción σ^z_i σ^z_{i+1}
    for i in range(N - 1):
        op_z = [I] * N
        op_z[i] = sz
        op_z[i + 1] = sz
        H -= J * producto_tensorial(op_z)
    
    # Término de campo transversal σ^x_i
    for i in range(N):
        op_x = [I] * N
        op_x[i] = sx
        H -= g * producto_tensorial(op_x)
    
    return H
```

### Estado inicial: todos los espines arriba

```python
def estado_inicial(N):
    arriba = np.array([[1], [0]])  # Estado |0⟩
    return producto_tensorial([arriba] * N)
```

### Evolución temporal del sistema

```python
def evolucion(psi0, H, tiempos):
    estados = []
    for t in tiempos:
        U = expm(-1j * H * t)
        psi_t = U @ psi0
        estados.append(psi_t)
    return estados
```

### Operador σ^z en el sitio i

```python
def sigma_z_en_i(N, i):
    ops = [I] * N
    ops[i] = sz
    return producto_tensorial(ops)
```

### Magnetización \( \langle \sigma^z_i \rangle (t) \)

```python
def magnetizacion_z(estados, N):
    magnetizaciones = np.zeros((len(estados), N))
    for t_idx, psi in enumerate(estados):
        for i in range(N):
            Sz_i = sigma_z_en_i(N, i)
            magnetizaciones[t_idx, i] = np.real((psi.conj().T @ Sz_i @ psi).item())
    return magnetizaciones
```

### Parámetros y ejecución de la simulación

```python
N = 4
J = 1.0
g = 1.0

H = hamiltoniano(N, J, g)
psi0 = estado_inicial(N)
tiempos = np.linspace(0, 5, 100)
estados = evolucion(psi0, H, tiempos)
magnetizaciones = magnetizacion_z(estados, N)
```

