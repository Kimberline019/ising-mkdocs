#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP
#include "Matriz.hpp"
#include <vector>
#include <complex>
/**
 * @brief Esta clase se encarga de realizar las operaciones necesarias para construir el Hamiltoniano del sistema a partir de los valores datos.
 *Almacena la matriz del Hamiltoniano y la cantidad de espines del sistema.
 * 
 */
class Hamiltoniano {
private:
    Matriz H;            
    int N;            
public:
    /**
    * @brief Constructor personalizado que crea una matriz sin inicializar de las dimensiones del Hamiltoniano del sistema.
    *
    * Ejemplo de uso:
    * @code
    * Hamiltoniano H_(espin);
    * @endcode
    */
    Hamiltoniano(int N_);
    /**
    * @brief Realiza los productos tensoriales necesarios para crear la matriz del Hamiltoniano del sistema y los almacena en los datos de la instancia.
    *
    * Ejemplo de uso:
    * @code
    * H_,construir(J,g);
    * @endcode
    *
    * @param J dato de la escala energética de la interacción ferromagnética.
    * @param g dato del parámetro energético del campo transversal.
    */
    void construir(double J, double g); 
        /**
    * @brief Devuelve el hamiltoniano del que la clase almacena.
    *
    * Ejemplo de uso:
    * @code
    * H=H_,matriz();
    * @endcode
    *
    * @return Matriz con el hamiltoniano de la clase.
    */
    const Matriz& matriz() const;
};

#endif

