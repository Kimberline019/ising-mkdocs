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
