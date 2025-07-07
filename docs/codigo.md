# Código del modelo de Ising cuántico en Python

Este apartado presenta la implementación completa del **modelo de Ising cuántico unidimensional con campo transversal**, un sistema fundamental en la física de muchos cuerpos. A continuación, se explica cada parte del código que simula su dinámica cuántica.

---

## Hamiltoniano del modelo

El Hamiltoniano está dado por:

$$
H = -J \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - g \sum_{i=1}^N \sigma^x_i
$$

donde:
- \( \sigma^x \) y \( \sigma^z \) son matrices de Pauli,  
- \( J \) es la constante de acoplamiento entre espines vecinos,  
- \( g \) es la intensidad del campo transversal,  
- \( N \) es el número total de espines en la cadena.

---

## Implementación del modelo en Python

```python
import numpy as np
from functools import reduce
from scipy.linalg import expm
```

### Matrices de Pauli y matriz identidad

```python
I = np.eye(2)
sx = np.array([[0, 1], [1, 0]])   # σ_x
sz = np.array([[1, 0], [0, -1]])  # σ_z
```

### Producto tensorial

```python
def producto_tensorial(matrices):
    return reduce(np.kron, matrices)
```

### Construcción del Hamiltoniano

```python
def hamiltoniano(N, J, g):
    H = np.zeros((2**N, 2**N))
    
    # Término de interacción σ^z_i σ^z_{i+1}
    for i in range(N - 1):
        op_z = [I] * N
        op_z[i] = sz
        op_z[i + 1] = sz
        H -= J * producto_tensorial(op_z)
    
    # Término de campo transversal σ^x_i
    for i in range(N):
        op_x = [I] * N
        op_x[i] = sx
        H -= g * producto_tensorial(op_x)
    
    return H
```

### Estado inicial: todos los espines arriba

```python
def estado_inicial(N):
    arriba = np.array([[1], [0]])  # Estado |0⟩
    return producto_tensorial([arriba] * N)
```

### Evolución temporal del sistema

```python
def evolucion(psi0, H, tiempos):
    estados = []
    for t in tiempos:
        U = expm(-1j * H * t)
        psi_t = U @ psi0
        estados.append(psi_t)
    return estados
```

### Operador σ^z en el sitio i

```python
def sigma_z_en_i(N, i):
    ops = [I] * N
    ops[i] = sz
    return producto_tensorial(ops)
```

### Magnetización \( \langle \sigma^z_i \rangle (t) \)

```python
def magnetizacion_z(estados, N):
    magnetizaciones = np.zeros((len(estados), N))
    for t_idx, psi in enumerate(estados):
        for i in range(N):
            Sz_i = sigma_z_en_i(N, i)
            magnetizaciones[t_idx, i] = np.real((psi.conj().T @ Sz_i @ psi).item())
    return magnetizaciones
```

### Parámetros y ejecución de la simulación

```python
N = 4
J = 1.0
g = 1.0

H = hamiltoniano(N, J, g)
psi0 = estado_inicial(N)
tiempos = np.linspace(0, 5, 100)
estados = evolucion(psi0, H, tiempos)
magnetizaciones = magnetizacion_z(estados, N)
```


##Implementación en C++

##Matriz.cpp

Este código implementa una clase Matriz que representa matrices de números complejos y define operaciones fundamentales como suma, multiplicación matricial, multiplicación por escalar y producto tensorial, todas optimizadas con paralelización mediante OpenMP; su propósito dentro del modelo de Ising cuántico unidimensional es facilitar la construcción eficiente del Hamiltoniano a través de productos tensoriales de matrices de Pauli, así como la evolución del estado cuántico mediante multiplicaciones repetidas de matrices durante la integración numérica de la ecuación de Schrödinger.

```cpp
#include <vector>
#include <complex>
#include <iostream>
#include "Matriz.hpp"
#include <omp.h>
Matriz::Matriz(int filas_,int columnas_): filas(filas_), columnas(columnas_),data(filas*columnas, {0.0,0.0}){}
Matriz::Matriz(int filas_,int columnas_,const std::vector<std::complex<double>>&valores): filas(filas_),columnas(columnas_), data(valores){}
Matriz Matriz::operator+(const Matriz &obj) const {
	Matriz newMatriz(filas,columnas);
	#pragma omp parallel for
	for (int i=0; i <filas;++i){
		for (int j=0; j <columnas;++j){
		newMatriz(i,j)=(*this)(i,j)+obj(i,j);
		}
	}
	return newMatriz;
}

int Matriz::sizef() const {
	return filas;
}
int Matriz::sizec() const {
	return columnas;
}
std::complex<double>& Matriz::operator()(int i, int j){
	return data[i*columnas+j];
}
const std::complex<double>& Matriz::operator()(int i, int j) const{
        return data[i*columnas+j];
}
const std::vector<std::complex<double>>& Matriz::datos() const {
	return data;
}
Matriz Matriz::operator*(const Matriz &obj) const{
   	 Matriz newMatriz(filas,obj.sizec());
	 #pragma omp parallel for
   	 for (int i = 0; i < filas; ++i){
        	for (int j = 0; j < obj.sizec(); ++j){
            		for (int k = 0; k < columnas; ++k){
                		newMatriz(i, j) += (*this)(i, k) * obj(k, j);
			}
		}
	 }
     	 return newMatriz;
}
Matriz Matriz::operator*(const std::complex<double>& escalar) const{
	Matriz newMatriz(filas,columnas);
	#pragma omp parallel for
	for (int i=0; i < filas; ++i){
		for (int j = 0; j < columnas; ++j){
			newMatriz(i,j)=(*this)(i, j) * escalar;
		}
	}
	return newMatriz;
}
Matriz Matriz::tensor(const Matriz &obj) const {
    int m = filas;
    int n = columnas;
    int p = obj.sizef();
    int q = obj.sizec();
    Matriz newMatriz(m*p,n*q);
    #pragma omp parallel for collapse(2)
    for (int i=0;i<m;++i){
        for (int j=0;j<n;++j){
            for (int k=0;k<p;++k){
                for (int l=0;l<q;++l){
                    newMatriz(i*p+k,j*q+l)=(*this)(i,j)*obj(k,l);
                }
            }
        }
    }

    return newMatriz;
}
```



##Matriz.hpp

```cpp
#ifndef MATRIZ_HPP
#define MATRIZ_HPP

#include <vector>
#include <complex>
#include <iostream>

class Matriz {
	private:
		int filas, columnas;
		std::vector<std::complex<double>> data;
	public:
		Matriz(int filas, int columnas, const std::vector<std::complex<double>>&valores);
		Matriz(int filas, int columnas);
		Matriz operator+(const Matriz &obj) const;
		Matriz operator*(const Matriz &obj) const;
		Matriz operator*(const std::complex<double> &escalar) const;
		Matriz tensor(const Matriz &obj) const;
		int sizef() const;
		int sizec() const;
		std::complex<double>& operator()(int i, int j);
   		const std::complex<double>& operator()(int i, int j) const;
		const std::vector<std::complex<double>>& datos() const;
};

#endif
```


##hamiltoniano.cpp

```cpp
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
```


##hamiltoniano.hpp


```cpp
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
```

##main.cpp


```cpp
#include <iostream>
#include <vector>
#include <complex>
#include "Matriz.hpp"
#include "hamiltoniano.hpp"
#include "rk4.hpp"
#include "pauli.hpp"
int main(){
	double J=1.0;
	double g=1.0;
	double tF=2.0;
	double h=0.0001;
	int spin = 3;
	int n=100;
    	int dim = 1 << spin; 
   	std::vector<std::complex<double>> inicial(dim, {1.0, 0.0});
	Matriz phi0(dim,1,inicial);
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
	phicopy=rk4.direct(phi0,n,tF);
	for (int i=0;i<phicopy.sizef();++i){
        	std::cout<<phicopy(i,0)<<std::endl;
        }
        std::cout << std::endl;

	return 0;
	}
```



##pauli.cpp

```cpp
#include "pauli.hpp"

Matriz Pauli::I(){
	return Matriz(2,2,{std::complex<double>(1,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(1,0)});
}
Matriz Pauli::X(){
return Matriz(2,2,{std::complex<double>(0,0),std::complex<double>(1,0),std::complex<double>(1,0),std::complex<double>(0,0)});
}
Matriz Pauli::Z(){
	return Matriz(2,2,{std::complex<double>(1,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(-1,0)});
```



##pauli.hpp

```cpp

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
```


##rk4.cpp

```cpp

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
```

##rk4.hpp

```cpp

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
```

