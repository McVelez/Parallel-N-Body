#include <iostream>
#include <fstream> 
#include <cmath>
#include <random>
#include <string> 

#include "nbody.h"

// Class methods
Nbody::Nbody(int n_, double dt_, double t_max_) : n(n_), dt(dt_), t_max(t_max_) {
    init_data();
}

// destruct
Nbody::~Nbody() {
    delete[] p; // reiniciar apunt de particulas
    p = 0; 
}

void Nbody::timeIntegration() {
    comp_force();

    for (; t < t_max; t += dt, step += 1) {
        comp_position();
        comp_force();
        comp_velocity();
    }

    print_parameter();
}


// no se necesita paralelizar porque solamente realiza 2 iteraciones
void Nbody::update_velocity(Particle* p) {
    double a = dt * 0.5 / p->m;
    for (int d = 0; d < DIM; d++) {
        p->v[d] += ((p->F[d] + p->F_old[d]) * dt * 0.5) / p->m; // para cada particula calc velocidad
        //bool x = p->F[d] == p->F_old[d]; // para revisar que las fuerzas no fueran iguales
        //printf(x ? "true" : "false");
    }
}

// no se necesita paralelizar porque solamente realiza 2 iteraciones
void Nbody::update_position(Particle* p) {
    for (int d = 0; d < DIM; d++) { 
        p->x[d] += dt * (p->v[d] + ((dt * 0.5 * p->F[d]/ p->m))); // calc posicion = dt * (velocidad + a * fuerza
        p->F_old[d] = p->F[d]; // old force is current force (chain effect)
    }
}

// Esto si puede ser paralelizado
void Nbody::comp_velocity() {
    #pragma omp parallel for schedule(dynamic, 50) num_threads(4)
    for (int i = 0; i < n; i++) {
        update_velocity(&p[i]); // calc velocidad para cada particula
    }
}

// Esto si puede ser paralelizado
void Nbody::comp_position() {
    #pragma omp parallel for schedule(dynamic, 50) num_threads(4)
    for (int i = 0; i < n; i++) { // por cada particula
        update_position(&p[i]); // actualizar posicion
    }
}

// Esto si puede ser paralelizado
void Nbody::comp_force() {
    #pragma omp parallel for schedule(dynamic, 50) num_threads(4)
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < DIM; d++) {
            p[i].F[d] = 0; // reiniciar las fuerzas ? seh
        }
    }

    #pragma omp parallel for schedule(dynamic, 50) num_threads(4)
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            force(&p[i], &p[j]); // calcular fuerza pora cada par de particulas
        }
    }
}

// no se necesita paralelizar porque solamente realiza 2 iteraciones max
void Nbody::force(Particle* i, Particle* j) {
    double r = EPS; // smoothing +/-
    for (int d = 0; d < DIM; d++) {
        r += sqr(j->x[d] - i->x[d]); // distancia entre dos particulas en cada eje
    }

    double f = i->m * j->m / (sqrt(r) * r); // m1 * m2 / sqrt(r) * r
    for (int d = 0; d < DIM; d++) {
        i->F[d] += f * (j->x[d] - i->x[d]); // elian se debe de acordar que el menos es por la direccion de la particula en contrario de la otra
        j->F[d] -= f * (j->x[d] - i->x[d]);
    }
}

void Nbody::print_data() const {
    std::cout.setf(std::ios_base::scientific);
    std::cout.precision(5);
    for (int i = 0; i < n; i++) {
        std::cout << t << " ";
        std::cout << p[i].m << " ";
        for (int d = 0; d < DIM; d++) {
            std::cout << p[i].x[d] << " ";
        }
        for (int d = 0; d < DIM; d++) {
            std::cout << p[i].v[d] << " ";
        }
        for (int d = 0; d < DIM; d++) {
            std::cout << p[i].F[d] << " ";
        }
        std::cout << std::endl;
    }
}

// Esto si puede ser paralelizado
void Nbody::init_data() {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution_x(0.0, 1.0);
    std::uniform_real_distribution<double> distribution_v(-1.0, 1.0);
    
    #pragma omp parallel for schedule(dynamic, 50) num_threads(8)
    for (int i = 0; i < n; i++) {
        p[i].m = 1. / n; // calc masa
        for (int d = 0; d < DIM; d++) {
            p[i].x[d] = i; // calc posicion
            p[i].v[d] = 0.4; // calc velocidad
            p[i].F[d] = 0.0; // init de fuerzas x y y de la particula 
            p[i].F_old[d] = 0.0; // init de fuerzas x y y un paso atras
        }
    }
}

inline void Nbody::print_parameter() const {
    for (int i = 0; i < n; i++)
    {
        std::cout << n << " " << *p[i].x << " " << *p[i].v << " " << *p[i].F << std::endl;
    }
}

// Other Functions

inline double sqr(double x) {
    return x * x;
}