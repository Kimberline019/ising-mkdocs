#include <iostream>
#include <vector>
#include <complex>
#include "Matriz.hpp"
#include "hamiltoniano.hpp"
#include "rk4.hpp"
#include "pauli.hpp"
#include <omp.h>
#include <sys/time.h>

double seconds() {
    struct timeval tmp;
    gettimeofday(&tmp, nullptr);
    return tmp.tv_sec + tmp.tv_usec * 1e-6;
}


int main(){
	double J=1.0;
	double g=1.0;
	double tF=2.0;
	double h=0.0001;
	int spin = 3;
	int n=100;

   int num_threads;
    std::cout << "Ingrese el número de hilos a usar: ";
    std::cin >> num_threads;
    omp_set_num_threads(num_threads); 


    	int dim = 1 << spin; 
   	std::vector<std::complex<double>> inicial(dim, {1.0, 0.0});
	Matriz phi0(dim,1,inicial);
	Matriz phicopy(dim,1);

	Hamiltoniano H_(spin);
	H_.construir(J,g);
	Matriz H=H_.matriz();
	
	rk4 rk4(H,h);

    double time_start = seconds();

	double contador=0;
	while (contador <=tF){
		phicopy=rk4.aplicar(phi0); 
		std::cout<<"t = "<< contador<< std::endl;
		for (int i=0;i< phicopy.sizef();++i){
    			std::cout << phicopy(i,0)<< std::endl;
			}
		std::cout << std::endl;
		contador+=h;
		phi0=phicopy;
	}
	phicopy=rk4.direct(phi0,n,tF);
	for (int i=0;i<phicopy.sizef();++i){
        	std::cout<<phicopy(i,0)<<std::endl;
        }
        std::cout << std::endl;
    
    double time_end = seconds();
    std::cout << "Tiempo de ejecución: " << (time_end - time_start) << " segundos" << std::endl;
    std::cout << "Número de hilos utilizados: " << num_threads << std::endl;
	return 0;
	}
