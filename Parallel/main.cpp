#include "nbody.h"
#include <omp.h>
#include <iostream>

int main(int argc, char* argv[]) {
    double start = omp_get_wtime();
    // num_particles, step size, max simulation time
    Nbody nbody(4500, 0.001, 1); // already parallelized
    nbody.timeIntegration();
    double total = omp_get_wtime() - start;
    printf("Sequential: %f", total);
    return 0;
}
