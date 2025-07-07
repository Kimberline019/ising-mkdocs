#include "pauli.hpp"

Matriz Pauli::I(){
	return Matriz(2,2,{std::complex<double>(1,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(1,0)});
}
Matriz Pauli::X(){
return Matriz(2,2,{std::complex<double>(0,0),std::complex<double>(1,0),std::complex<double>(1,0),std::complex<double>(0,0)});
}
Matriz Pauli::Z(){
	return Matriz(2,2,{std::complex<double>(1,0),std::complex<double>(0,0),std::complex<double>(0,0),std::complex<double>(-1,0)});
}


