#include "hamiltoniano.hpp"
#include "pauli.hpp" 
#include <iostream>

Hamiltoniano::Hamiltoniano(int N_): N(N_),H(1 << N_, 1 << N_) {
    int dim =1<<N; 
    H=Matriz(dim, dim);
}
void Hamiltoniano::construir(double J, double g) {
    int dim =1 <<N;
    H=Matriz(dim,dim);
    for (int i=0;i<N;++i){
        int ip =(i+1)%N;
        Matriz term_zz=Matriz(1,1,{std::complex<double>(1,0)});
        for (int site=0;site<N;++site){
            if (site == i)
                term_zz=term_zz.tensor(Pauli::Z());
            else if (site==ip)
                term_zz=term_zz.tensor(Pauli::Z());
            else
                term_zz = term_zz.tensor(Pauli::I());
        }
        H=H+term_zz*(-J);
    }
    for (int i=0;i<N;++i){
        Matriz term_x=Matriz(1,1,{std::complex<double>(1,0)});
        for (int site=0;site<N;++site){
            if (site==i)
                term_x=term_x.tensor(Pauli::X());
            else
                term_x=term_x.tensor(Pauli::I());
        }
        H=H+term_x*(-g);
    }
}

const Matriz& Hamiltoniano::matriz() const {
    return H;
}

