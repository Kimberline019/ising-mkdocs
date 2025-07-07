#ifndef PAULI_H
#define PAULI_H
#include "Matriz.hpp"
/**
 * @brief Define las matrices de Pauli usando la clase Matriz como constantes para su fácil acceso.
 *
 */
class Pauli {
public:
    static Matriz I();
    static Matriz X();
    static Matriz Z();
};

#endif

