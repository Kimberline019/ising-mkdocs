#ifndef PAULI_H
#define PAULI_H
#include "Matriz.hpp"

class Pauli {
public:
    static Matriz I();
    static Matriz X();
    static Matriz Z();
};

#endif

