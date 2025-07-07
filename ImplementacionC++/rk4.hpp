#ifndef RK4_HPP
#define RK4_HPP

#include <vector>
#include <complex>
#include "Matriz.hpp"
/**
 * @brief Almacena los métodos de cálculo para el phi final.
 *
 * Sus operadores permiten calcular limpiamente el estado final de phi. Almacena el hamiltoniano del sistema y el tamaño del paso a aplicar.
 */	
class rk4 {
private:
    Matriz H;
    double h;
public:
    /**
     * @brief Constructor personalizado que almacena los datos necesarios para el método numérico.
     * Ejemplo de uso: 
     * @code 
     * rk4 metodo(H,0.001);
     *  @endcode
     *
     * @param &Hamiltoniano El hamiltoniano del sistema.
     * @param Paso Tamaño del paso que aplicará el método numérico.
     */
    rk4(const Matriz& Hamiltoniano, double paso);
        /**
     * @brief Aplica una iteración del método Runge-Kutta grado 4.
     * Ejemplo de uso:
     * @code
     * Matriz v= metodo.aplicar(phi_0);
     *  @endcode
     *
     * @param &phi0 Matriz con el vector columna que almacena las condiciones de phi antes de la iteración.
     * @return &phi0 Matriz con el vector columna phi en el tiempo t+h.
     */
    Matriz aplicar(const Matriz &phi0) const;
     /**
     * @brief Calcula mediante series de Taylor el valor de exp(-iHt) y la aplica las condiciones iniciales de phi.
     * Ejemplo de uso:
     * @code
     * Matriz v= metodo.direct(phi_0);
     *  @endcode
     *
     *@param &phi0 Matriz con el vector columna que almacena las condiciones de
 phi antes de la iteración.
     *@param &n Entero con el número de términos que se quieren para la expansión de Taylor.
     *@param &t Tiempo final del sistema.
     *@return Matriz con el vector columna con los valores finales de phi.
     *
     * 
     */
    Matriz direct(const Matriz &phi0, int &n,double &t) const;
};

#endif
