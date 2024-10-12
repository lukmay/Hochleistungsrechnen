#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void calculate_pi(int num_points, int rank, int numProcs) {
    // WTime Initialization
    double starttime, endtime;
    starttime = MPI_Wtime();
    // Define how many points each rank generates and generate them
    int points_for_rank = round(num_points / numProcs);
    int points_in_circle_in_rank = 0;
    int points_in_circle_global = 0;

    for (int count = 0; count < points_for_rank; count++) {
        double random_x = (double)rand() / RAND_MAX;
        double random_y = (double)rand() / RAND_MAX;

        if (random_x * random_x + random_y * random_y <= 1.0) {
            points_in_circle_in_rank++;
        }
    }
    // Gather results from all ranks
    MPI_Reduce(&points_in_circle_in_rank, &points_in_circle_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double pi_approx = 4.0 * points_in_circle_global / num_points;
        endtime = MPI_Wtime();
        printf("WTime: %f seconds\n", endtime - starttime);
        printf("Approximation of Pi: %f\n", pi_approx);
    }    
}

int main(int argc, char *argv[]) {
    // MPI Initialization
    MPI_Init(&argc, &argv);
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    int total_points = atoi(argv[1]);

    if (rank == 0) {
        printf("Total Points used for Approximation: %i \n", total_points);
    }

    calculate_pi(total_points, rank, numProcs);

    MPI_Finalize();
    return 0;
}
