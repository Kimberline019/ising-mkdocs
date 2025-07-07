#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP
#include "Matriz.hpp"
#include <vector>
#include <complex>
class Hamiltoniano {
private:
    Matriz H;            
    int N;            
public:
    Hamiltoniano(int N_);
    void construir(double J, double g); 
    const Matriz& matriz() const;
};

#endif

