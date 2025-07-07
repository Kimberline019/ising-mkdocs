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


## Implementación en C++

 HEAD
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



## Matriz.hpp


```cpp

#ifndef MATRIZ_HPP 
#define MATRIZ_HPP 

#include <vector>
#include <complex>
#include <iostream>
/**
 * @brief Representa matrices. 
 * Simplifica cálculos con matrices mediante operadores, lo que evita for loops explícitos en el resto del código.
 * Almacena las dimensiones de la matriz en variables int  y sus datos en un std::vector<std::complex<double>>
 */
class Matriz {
	private:
		int filas, columnas;
		std::vector<std::complex<double>> data;
	public:
		/**
	 	* @brief Constructor personalizado que representa una matriz con el tamaño y datos dados.
	 	* Ejemplo de uso:
	 	* @code
	 	*  Matriz B(2,2,({1,0},{0,0},{0,1},{0,1});
	 	* @endcode
	 	*
	 	* @param filas  Número de filas.
	 	* @param columnas  Número de columnas.
	 	* @param vector complejo con los datos aplanados de la matriz.
	 	*/
		Matriz(int filas, int columnas, const std::vector<std::complex<double>>&valores);
		/**
	 	* @brief Constructor personalizado que representa una matriz no inicializada con el tamaño dado.
	 	*
	 	* Ejemplo de uso:
	 	* @code
	 	* Matriz A(2,2);
	 	* @endcode
	 	*
	 	* @param filas Número de filas.
	 	* @param columnas Número de columnas.
	 	*/
		Matriz(int filas, int columnas);
		 /**
		 * @brief Realiza la suma de dos matrices.
		 *
		 * Ejemplo de uso:
		 * @code
		 * Matriz C=A+B;
		 * @endcode
		 *
		 * @param &obj Matriz a la derecha de la suma.
		 * @return La suma de ambas matrices.
		 */
		Matriz operator+(const Matriz &obj) const;
                 /**
                 * @brief Realiza el producto entre dos matrices.
                 *
                 * Ejemplo de uso:
                 * @code
                 * Matriz C=A*B;
                 * @endcode
                 *
                 * @param &obj Matriz a la derecha del producto.
                 * @return El producto matricial.
                 */
		Matriz operator*(const Matriz &obj) const;
                 /**
                 * @brief Realiza el producto escalar.
                 *
                 * Ejemplo de uso:
                 * @code
                 * Matriz C=A*b;
                 * @endcode
                 *
                 * @param &escalar escalar por el que se va a multiplicar (debe siempre estar a la derecha de la matriz).
                 * @return Resultado del producto.
                 */
		Matriz operator*(const std::complex<double> &escalar) const;
		 /**
                 * @brief Realiza el producto tensorial.
                 *
                 * Ejemplo de uso:
                 * @code
                 * Matriz C=A.tensor(B);
                 * @endcode
                 *
                 * @param &obj Matriz derecha por la que se realizará el producto tensorial.
                 * @return Matriz resultado del producto tensorial.
                 */
		Matriz tensor(const Matriz &obj) const;
		 /**
                 * @brief Da la cantidad de filas de la matriz.
                 *
                 * Ejemplo de uso:
                 * @code
                 * int b=A.sizef();
		 * //b=2
                 * @endcode
                 *
                 * @return Número de filas de la matriz.
                 */
		int sizef() const;
		 /**
                 * @brief Da la cantidad de columnas de la matriz.
                 *
                 * Ejemplo de uso:
                 * @code
                 * int b=A.sizec();
                 * //b=2
                 * @endcode
                 *
                 * @return Número de columnas de la matriz.
                 */
		int sizec() const;
		 /**
                 * @brief Permite accesar un dato de la matriz aplanada de forma intuitiva para aplicaciones de la clase.
                 *
                 * Ejemplo de uso:
                 * @code
                 * std::complex<double>> b=A(2,2);
                 * //b={1.0,0.0}
                 * @endcode
                 *
		 * @param i Fila del dato a accesar.
		 * @param j Columna del dato a accesar.
                 * @return Dato en la posición (i,j) de la matriz.
                 */
		std::complex<double>& operator()(int i, int j);
		 /* @brief Permite accesar un dato de la matriz aplanada de forma intuitiva para aplicaciones externas.
                 *
                 * Ejemplo de uso:
                 * @code
                 * std::complex<double>> b=A(2,2);
                 * //b={1.0,0.0}
                 * @endcode
                 *
                 * @param i Fila del dato a accesar.
                 * @param j Columna del dato a accesar.
                 * @return Dato en la posición (i,j) de la matriz.
                 */
   		const std::complex<double>& operator()(int i, int j) const;
		 /* @brief Devuelve el vector complejo que almacena los datos de la matriz.
                 *
                 * Ejemplo de uso:
                 * @code
                 * std::vector<std::complex<double>> b=A.datos()
                 * //b=({1,0},{0,0},{0,1},{0,1})
                 * @endcode
                 *
                 * @param i Fila del dato a accesar.
                 * @param j Columna del dato a accesar.
                 * @return Dato en la posición (i,j) de la matriz.
                 */
		const std::vector<std::complex<double>>& datos() const;
                 /* @brief Devuelve la norma de un vector columna.
                 *
                 * Ejemplo de uso:
                 * @code
                 * double b=A.norma()
                 * //b=sqrt(5)
                 * @endcode
                 *
                 * @return Norma del vector.
                 */
		double norma() const;
};

```


## hamiltoniano.cpp
El archivo hamiltoniano.cpp define una clase Hamiltoniano que construye la matriz Hamiltoniana del modelo de Ising cuántico unidimensional en una red de $$ ( N )$$ espines, empleando productos tensoriales de matrices de Pauli. El término de interacción $$ ( -J \sum_i \sigma^z_i \sigma^z_{i+1} )$$ se implementa mediante ciclos que insertan dos matrices$$ (\sigma^z)$$ en las posiciones correspondientes y matrices identidad en los sitios restantes, considerando condiciones periódicas. El término de campo transversal $$ ( -g \sum_i \sigma^x_i )$$  se construye de forma similar, insertando$$ (\sigma^x)$$  en la posición adecuada. La matriz resultante $$ ( \hat{H} )$$ es de dimensión $$( 2^N \times 2^N )$$ y se guarda internamente para su posterior uso en la evolución dinámica del sistema.
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


## hamiltoniano.hpp


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

## main.cpp

El archivo **main.cpp** realiza la simulación numérica de la dinámica del modelo de Ising  a través de un enfoque de evolución temporal. El usuario ingresa el número de hilos a utilizar mediante OpenMP para aprovechar la paralelización del cálculo. A partir de un estado inicial puro, se construye la matriz Hamiltoniana  \(\hat{H}\)  del sistema utilizando productos tensoriales de matrices de Pauli y se resuelve la ecuación de Schrödinger dependiente del tiempo utilizando el método de Runge-Kutta de cuarto orden, implementado en la clase **rk4**. Además, se evalúa la evolución unitaria exacta mediante una expansión de Taylor truncada, permitiendo comparar y validar el método numérico. Finalmente, se mide el tiempo total de ejecución para evaluar el desempeño computacional según el número de hilos empleados, lo cual es esencial para el análisis de aceleración (**speedup**) en cálculos paralelos.


```cpp
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
```



## pauli.cpp

El archivo $$ \texttt{pauli.cpp}$$ define funciones estáticas de la clase Pauli que retornan las matrices de Pauli $$( \hat{I} )$$, $4( \hat{\sigma}^x )$$  y $$ ( \hat{\sigma}^z )$$  representadas como objetos de la clase Matriz. Estas matrices son fundamentales en la construcción del Hamiltoniano del modelo de Ising cuántico unidimensional, ya que $$( \hat{\sigma}^z )$$  se utiliza para modelar la interacción entre espines vecinos, mientras que $$ ( \hat{\sigma}^x )$$ representa el efecto del campo magnético transversal. La matriz identidad $$( \hat{I} )$$ se emplea para completar los productos tensoriales cuando un operador actúa sobre un sitio específico de la cadena de espines

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



## pauli.hpp

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


## rk4.cpp


El archivo rk4.cpp implementa la clase rk4, que permite resolver numéricamente la ecuación de Schrödinger dependiente del tiempo para un sistema cuántico, utilizando dos métodos: el método de Runge-Kutta de cuarto orden y una aproximación directa de la evolución unitaria basada en el desarrollo en serie de Taylor de la exponencial de matrices. El método aplicar calcula un paso de evolución temporal aplicando el algoritmo de Runge-Kutta a la función de onda $( |\psi(t)\rangle )$, con un paso de tiempo ( h ), usando la matriz Hamiltoniana $( \hat{H} )$. Por otro lado, el método direct evalúa directamente  $( e^{-i\hat{H}t}|\psi(0)\rangle )$ truncando la serie de potencias en $( n ) términos, lo cual es útil para validar la precisión del método numérico. Ambas estrategias son fundamentales para estudiar la dinámica del modelo de Ising a partir de un estado inicial.
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

## rk4.hpp

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

