import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
from scipy.linalg import expm

# Matrices de Pauli y la identidad
I = np.eye(2)
sx = np.array([[0, 1], [1, 0]])  # sigma_x
sz = np.array([[1, 0], [0, -1]]) # sigma_z

def producto_tensorial(matrices):
    """
    Calcula el producto tensorial de una lista de matrices.

    Args:
        matrices (list of np.ndarray): Lista de matrices 2D de numpy.

    Returns:
        np.ndarray: Matriz resultante del producto tensorial.

    Example:
        >>> producto_tensorial([I, sx]).shape
        (4, 4)
    """
    return reduce(np.kron, matrices)

def hamiltoniano(N, J, g):
    """
    Construye el Hamiltoniano del modelo con campo transversal.

    Args:
        N (int): Número de espines.
        J (float): Constante de energía.
        g (float): Intensidad del campo transversal.

    Returns:
        np.ndarray: Matriz Hamiltoniana de dimensión 2^N x 2^N.

    Example:
        >>> H = hamiltoniano(4, 1.0, 1.0)
    """
    H = np.zeros((2**N, 2**N), dtype=complex)

    for i in range(N - 1):
        op_z = [I] * N
        op_z[i] = sz
        op_z[i + 1] = sz
        H -= J * producto_tensorial(op_z)

    for i in range(N):
        op_x = [I] * N
        op_x[i] = sx
        H -= g * producto_tensorial(op_x)

    return H

def estado_inicial(N):
    """
    Construye el estado inicial con todos los espines hacia arriba.

    Args:
        N (int): Número de espines.

    Returns:
        np.ndarray: Vector columna de dimensión 2^N representando |0...0⟩.

    Example:
        >>> psi0 = estado_inicial(3)
    """
    arriba = np.array([[1], [0]])
    psi0 = producto_tensorial([arriba] * N)
    return psi0

def sigma_z_en_i(N, i):
    """
    Construye el operador sigma z actuando en el sitio i de una cadena de N espines.

    Args:
        N (int): Número de espines.
        i (int): Índice del sitio.

    Returns:
        np.ndarray: Matriz 2^N x 2^N del operador sigma_z en el sitio i.

    Example:
        >>> Sz_0 = sigma_z_en_i(3, 0)
        >>> Sz_0.shape
        (8, 8)
    """
    ops = [I] * N
    ops[i] = sz
    return producto_tensorial(ops)

def evolucion(psi0, H, tiempos):
    """
    Calcula la evolución temporal del estado cuántico bajo el Hamiltoniano dado.

    Args:
        psi0 (np.ndarray): Estado inicial (vector columna).
        H (np.ndarray): Hamiltoniano.
        tiempos (np.ndarray): Array de tiempos en los que evaluar.

    Returns:
        list of np.ndarray: Lista de estados evolucionados |ψ(t)⟩.

    Example:
        >>> estados = evolucion(psi0, H, np.linspace(0, 1, 5))
    """
    estados = []
    for t in tiempos:
        U = expm(-1j * H * t)
        psi_t = U @ psi0
        estados.append(psi_t)
    return estados

def magnetizacion_z(estados, N):
    """
    Calcula la magnetización en z sigma:_z(t) para cada espín y cada instante.

    Args:
        estados (list of np.ndarray): Lista de estados phi(t).
        N (int): Número de espines.

    Returns:
        np.ndarray: Matriz (len(estados), N) con ⟨sigma_z_i⟩(t).

    Example:
        >>> magz = magnetizacion_z(estados, 4)
    """
    magnetizaciones = np.zeros((len(estados), N))
    for t_idx, psi in enumerate(estados):
        for i in range(N):
            Sz_i = sigma_z_en_i(N, i)
            magnetizaciones[t_idx, i] = np.real((psi.conj().T @ Sz_i @ psi).item())
    return magnetizaciones

