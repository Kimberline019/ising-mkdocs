#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include "Matriz.hpp"
#include "rk4.hpp"

rk4::rk4(const Matriz& Hamiltoniano, double paso):H(Hamiltoniano),h(paso) {}
Matriz rk4::aplicar(const Matriz &phi0) const{
    std::complex<double> z(0.0, -1.0);
    std::complex<double> coef= z*h;
    Matriz k1(phi0.sizef(),1),k2(phi0.sizef(),1),k3(phi0.sizef(),1),k4(phi0.sizef(),1);
    k1=H*phi0*coef;
    k2=H*(phi0+k1*0.5)*coef;
    k3=H*(phi0+k2*0.5)*coef;
    k4=H*(phi0+k3)*coef;
    return phi0+(k1+k2*2+k3*2+k4)*(1.0/6.0);
}
Matriz rk4::direct(const Matriz& phi0,int &n,double &t) const{
    	std::complex<double> z(0.0, -1.0);
   	 Matriz phi1(phi0.sizef(), 1);
    	 Matriz identidad(H.sizef(), H.sizec());
    	for (int i = 0; i < H.sizef(); ++i) {
        	identidad(i, i) = {1.0, 0.0};
    	}
    	Matriz suma = identidad;
    	Matriz potenciaH = identidad;
    	std::complex<double> coef = 1.0;  // Coeficiente factorial
    	for (int k=1;k<=n;++k){
        	potenciaH = potenciaH*H;
        	coef /= static_cast<double>(k);
        	suma =suma+potenciaH*(coef*std::pow(z*t,k));
    	}

    	phi1 = suma*phi0;
    	return phi1;
}
