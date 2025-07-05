#Bienvenido

# Modelo de Ising Cuántico

El modelo de Ising cuántico unidimensional es un sistema fundamental en física de muchos cuerpos y teoría cuántica de espines. Este modelo describe una cadena de \( N \) espines interactuando con sus vecinos más cercanos y con un campo magnético transversal.

## Objetivo

Estudiar la evolución temporal del estado cuántico \( \mid \psi(t) \rangle \) gobernado por la ecuación de Schrödinger:

$$
\mid \psi(t) \rangle = e^{-iHt} \mid \psi(0) \rangle
$$

donde \( H \) es el Hamiltoniano del sistema.

## Hamiltoniano del modelo

El Hamiltoniano del modelo de Ising cuántico unidimensional con campo transversal se escribe como:

$$
H = -J \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - h \sum_{i=1}^{N} \sigma^x_i
$$

donde:

- \( \sigma^z_i \) y \( \sigma^x_i \) son matrices de Pauli que actúan sobre el espín \( i \),
- \( J \) representa la constante de acoplamiento (interacción entre espines vecinos),
- \( h \) representa la intensidad del campo magnético transversal.

## Matrices de Pauli

Las matrices de Pauli son operadores \( 2 \times 2 \) fundamentales que actúan sobre espines individuales:

$$
\sigma^x = 
\begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix},
\quad
\sigma^y = 
\begin{pmatrix}
0 & -i \\
i & 0
\end{pmatrix},
\quad
\sigma^z = 
\begin{pmatrix}
1 & 0 \\
0 & -1
\end{pmatrix}
$$

También se utiliza la **matriz identidad**:

$$
I = 
\begin{pmatrix}
1 & 0 \\
0 & 1
\end{pmatrix}
$$

## Observables relevantes

Podemos estudiar:

### Magnetización total

$$
\langle \psi(t) \mid \sum_i \sigma^z_i \mid \psi(t) \rangle
$$

### Correlación entre espines

$$
\langle \psi(t) \mid \sigma^z_i \sigma^z_j \mid \psi(t) \rangle
$$

---

## Implementación numérica

Para implementar este modelo en código (por ejemplo, en C++), es necesario:

- Representar las matrices de Pauli y la Identidad como matrices `2x2`.
- Usar el producto tensorial para construir las matrices que actúan sobre el espacio completo de dimensión \( 2^N \).
- Construir el Hamiltoniano sumando los términos de interacción y campo.
- Diagonalizar el Hamiltoniano para obtener autovalores y autovectores.
- Evolucionar el estado en el tiempo con \( e^{-iHt} \mid \psi(0) \rangle \).

Este procedimiento es **computacionalmente costoso** para \( N > 10 \), por lo que se aplica principalmente a sistemas pequeños.
