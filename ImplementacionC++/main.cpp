#include <iostream>
#include <vector>
#include <complex>
#include "Matriz.hpp"
#include "hamiltoniano.hpp"
#include "rk4.hpp"
#include "pauli.hpp"
/* @brief Da los valores de todas las cantidades relevantes para la resolución del problema y aplica los métodos.
 */

int main(){
	double J=1.0; // Valor de la escala energética
	double g=1.0; //Valor del parámetro energético del campo transversal
	double tF=2.0; //Tiempo final
	double h=0.01; //Tamaño del paso
	int spin = 1;  //Cantidad de espines
	int n=2; //Cantidad de términos para la expansión de Taylor 
    	int dim = 1 << spin; 
   	std::vector<std::complex<double>> inicial(dim, {1.0, 0.0}); //Condiciones iniciales de phi
	Matriz phi0(dim,1,inicial);
	Matriz phi_0(dim,1,inicial);
	Matriz phicopy(dim,1);
	Hamiltoniano H_(spin);
	H_.construir(J,g);
	Matriz H=H_.matriz();
	rk4 rk4(H,h);
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
	std::cout << phi0.norma()<< std::endl;
	phicopy=rk4.direct(phi_0,n,tF);
	for (int i=0;i<phicopy.sizef();++i){
        	std::cout<<phicopy(i,0)<<std::endl;
        }
        std::cout << std::endl;
	std::cout << phicopy.norma()<< std::endl;

	return 0;
	}
