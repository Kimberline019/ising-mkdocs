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
