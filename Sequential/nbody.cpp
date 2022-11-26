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

Nbody::~Nbody() {
    delete[] p;
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

void Nbody::update_velocity(Particle* p) {
    double a = dt * 0.5 / p->m;
    for (int d = 0; d < DIM; d++) {
        p->v[d] += a * (p->F[d] + p->F_old[d]);
    }
}

void Nbody::update_position(Particle* p) {
    double a = dt * 0.5 / p->m;
    for (int d = 0; d < DIM; d++) {
        p->x[d] += dt * (p->v[d] + a * p->F[d]);
        p->F_old[d] = p->F[d];
    }
}

void Nbody::comp_velocity() {
    for (int i = 0; i < n; i++) {
        update_velocity(&p[i]);
    }
}

void Nbody::comp_position() {
    for (int i = 0; i < n; i++) {
        update_position(&p[i]);
    }
}

void Nbody::comp_force() {
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < DIM; d++) {
            p[i].F[d] = 0;
        }
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            force(&p[i], &p[j]);
        }
    }
}

void Nbody::force(Particle* i, Particle* j) {
    double r = EPS; // smoothing
    for (int d = 0; d < DIM; d++) {
        r += sqr(j->x[d] - i->x[d]);
    }
    double f = i->m * j->m / (sqrt(r) * r);
    for (int d = 0; d < DIM; d++) {
        i->F[d] += f * (j->x[d] - i->x[d]);
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

void Nbody::init_data() {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution_x(0.0, 1.0);
    std::uniform_real_distribution<double> distribution_v(-1.0, 1.0);

    for (int i = 0; i < n; i++) {
        p[i].m = 1. / n;
        for (int d = 0; d < DIM; d++) {
            p[i].x[d] = i;
            p[i].v[d] = 0.4;
            p[i].F[d] = 0.0;
            p[i].F_old[d] = 0.0;
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