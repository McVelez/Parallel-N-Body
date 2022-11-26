#include "nbody.h"
#include <omp.h>
#include <iostream>

int main(int argc, char* argv[]) {
    // num_particles, step size, max simulation time
    double start = omp_get_wtime();
    Nbody nbody(4500, 0.001, 1);
    nbody.timeIntegration();
    double total = omp_get_wtime() - start;
    printf("Sequential %f", total);
    return 0;
}
