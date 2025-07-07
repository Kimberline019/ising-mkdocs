#Bienvenido
# Modelo de Ising cuántico unidimensional en una grilla de \( N \) espines

El **modelo de Ising cuántico** es una extensión natural del **modelo de Ising clásico**, uno de los más simples y famosos para describir **materiales magnéticos**, propuesto por Ernst Ising en 1925.

Para el modelo clásico, se tiene una red donde cada punto contiene un espín que puede tomar solo dos valores \( s_i = \pm1 \). El sistema busca minimizar su energía total dependiendo de cómo interactúan los espines. Para medir la energía total de una configuración se usa el siguiente **Hamiltoniano clásico**:

\[
H = -J \sum_{\langle i, j\rangle} s_i s_j - h \sum_{i} s_i
\]

Donde:

- \( s_i = \pm 1 \) es el espín clásico en el sitio \( i \),
- \( J \) es la constante de acoplamiento:
  - \( J > 0 \): favorece que los espines se alineen (**ferromagnetismo**),
  - \( J < 0 \): favorece que se desalineen (**antiferromagnetismo**),
- \( h \) es el campo magnético externo clásico.

---

## Hamiltoniano del modelo cuántico

El **modelo de Ising cuántico** permite estudiar el comportamiento colectivo de muchos espines cuánticos en una red, donde los espines pueden interactuar entre sí \( \text{(alinearse o desalinearse)} \). El **Hamiltoniano del modelo cuántico 1D con campo transversal** (perpendicular al eje de los espines) es:

\[
\hat{H} = -J\sum_{i=1}^{N-1} \hat{\sigma}_{i}^{z} \hat{\sigma}_{i + 1}^z - g \sum_{i = 1}^N \hat{\sigma}_{i}^x
\]

Donde:

- \( J \): escala energética de la interacción ferromagnética,
- \( g \): intensidad del campo transversal cuántico,
- \( \hat{\sigma}_{i}^{z} \): matriz de Pauli-Z que mide si un espín está en \( \uparrow \) o \( \downarrow \),
- \( \hat{\sigma}_{i}^{x} \): matriz de Pauli-X que representa la posibilidad de cambio de estado \( \left(\uparrow \leftrightarrow \downarrow \right) \).

---

## Matrices de Pauli

\[
\sigma_x =
\begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix}, \quad
\sigma_y =
\begin{pmatrix}
0 & -i \\
i & 0
\end{pmatrix}, \quad
\sigma_z =
\begin{pmatrix}
1 & 0 \\
0 & -1
\end{pmatrix}, \quad
I =
\begin{pmatrix}
1 & 0 \\
0 & 1
\end{pmatrix}
\]

---

## Interpretación física del Hamiltoniano

**Primer término** \( \hat{\sigma}_{i}^{z} \hat{\sigma}_{i+1}^{z} \):
- Representa la interacción entre espines vecinos.
- Si están alineados \( (\uparrow\uparrow \text{ o } \downarrow\downarrow) \), da baja energía.
- Si están opuestos \( (\uparrow\downarrow \text{ o } \downarrow\uparrow) \), da mayor energía.

**Segundo término** \( \hat{\sigma}_{i}^x \):
- Representa la acción de un campo en la dirección \( x \),
- Introduce **fluctuaciones cuánticas** que permiten que los espines cambien de estado.

---

##  Competencia entre orden y fluctuación

- Cuando \( h = 0 \): sólo hay interacción entre espines → tienden a estar **alineados** → estado **ferromagnético**.
- Cuando \( J = 0 \): los espines siguen solo el campo externo → actúan **independientemente** → estado **paramagnético**.
- Cuando \( h \sim J \): hay un equilibrio → aparece una **transición cuántica de fase**.

---

##  Representación en mecánica cuántica

Cada espín vive en un espacio de dimensión 2 ( \( \uparrow \) o \( \downarrow \) ). Una cadena de \( N \) espines tiene un espacio de estados de dimensión \( 2^N \). Por ejemplo:

- \( N = 1 \): estados = \( \{ \uparrow, \downarrow \} \)
- \( N = 2 \): estados = \( \{ \uparrow\uparrow, \uparrow\downarrow, \downarrow\uparrow, \downarrow\downarrow \} \)
- \( N = 3 \): estados = \( \{ \uparrow\uparrow\uparrow, \uparrow\uparrow\downarrow, \ldots, \downarrow\downarrow\downarrow \} \Rightarrow 8 \) posibles

El Hamiltoniano es entonces una **matriz de dimensión \( 2^N \times 2^N \)** que describe todas las interacciones posibles del sistema.

---

##  Dinámica cuántica

A diferencia del modelo clásico, en el modelo cuántico el Hamiltoniano también **describe la evolución temporal** del estado del sistema:

$$
\frac{\partial |\psi(t)\rangle}{\partial t} = -i\hat{H} |\psi(t)\rangle
$$



La solución de esta ecuación es:


$$
|\psi(t)\rangle = e^{-i \hat{H}(t - t_0)} |\psi(t = t_0)\rangle
$$



Este formalismo permite estudiar:

- Dinámica de espines y correlaciones en el tiempo,
- Entrelazamiento cuántico,
- Simulación cuántica de materiales y algoritmos cuánticos.



---

## Implementación numérica

Para implementar este modelo en código (por ejemplo, en C++), es necesario:

- Representar las matrices de Pauli y la Identidad como matrices `2x2`.
- Usar el producto tensorial para construir las matrices que actúan sobre el espacio completo de dimensión \( 2^N \).
- Construir el Hamiltoniano sumando los términos de interacción y campo.
- Diagonalizar el Hamiltoniano para obtener autovalores y autovectores.
- Evolucionar el estado en el tiempo con \( e^{-iHt} \mid \psi(0) \rangle \).

Este procedimiento es **computacionalmente costoso** para \( N > 10 \), por lo que se aplica principalmente a sistemas pequeños.
