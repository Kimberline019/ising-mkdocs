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

#endif
