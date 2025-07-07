#ifndef RK4_HPP
#define RK4_HPP

#include <vector>
#include <complex>
#include "Matriz.hpp"

class rk4 {
private:
    Matriz H;
    double h;
public:
    rk4(const Matriz& Hamiltoniano, double paso);
    Matriz aplicar(const Matriz &phi0) const;
    Matriz direct(const Matriz &phi0, int &n,double &t) const;
};

#endif
