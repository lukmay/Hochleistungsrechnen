#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);
void releaseVector(Vector m);
void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
  // 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 500;

  // MPI Initialization
  MPI_Init(&argc, &argv);
  int rank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  // Time measurement variables
  double starttime, endtime;

  // Calculate chunk size and handle remainder
  int chunkSize = N / numProcs;
  int remainder = N % numProcs;

  // Determine start and end index for each rank
  int startIdx = rank * chunkSize + (rank < remainder ? rank : remainder);
  int endIdx = startIdx + chunkSize + (rank < remainder ? 1 : 0);
  int localSize = endIdx - startIdx;

  // Debug information
  printf("[DEBUG] Rank: %i, localSize: %i, startIdx: %i, endIdx: %i \n", rank, localSize, startIdx, endIdx);

  // Create a buffer for storing temperature fields
  Vector A = createVector(N);

  // Counts and displacements for MPI_Gatherv
  int *counts = NULL;
  int *displs = NULL;
  if (rank == 0) {
    counts = malloc(sizeof(int) * numProcs);
    displs = malloc(sizeof(int) * numProcs);
    for (int i = 0; i < numProcs; i++) {
      int tempStartIdx = i * chunkSize + (i < remainder ? i : remainder);
      int tempEndIdx = tempStartIdx + chunkSize + (i < remainder ? 1 : 0);
      counts[i] = tempEndIdx - tempStartIdx;
      displs[i] = tempStartIdx;
    }
  }

  int source_x = N / 4;

  if (rank == 0) {
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup ---------
    // Set up initial conditions in A
    for (int i = 0; i < N; i++) {
      A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // Heat source in one corner
    A[source_x] = 273 + 60;

    printf("Initial:\t");
    printTemperature(A, N);
    printf("\n");
  }

  // Distribute A to all processes
  MPI_Bcast(A, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  printf("[DEBUG] Broadcasted A to every rank for computation. \n");

  // Start time measurement before the computation begins
  starttime = MPI_Wtime();

  // ---------- compute ----------

  // Create buffers for the computation and initialize with A
  Vector B = createVector(localSize);
  Vector C = createVector(localSize);
  for (int i = startIdx; i < endIdx; i++) {
    B[i - startIdx] = A[i];
  }

  // Create variables for the neighbours
  value_t rightNeighbour = 0;
  value_t leftNeighbour = 0;

  // For each time step
  for (int t = 0; t < T; t++) {
    // Exchange boundary values with neighboring processes
    if (rank > 0) {
      MPI_Sendrecv(&B[0], 1, MPI_DOUBLE, rank - 1, 0, &leftNeighbour, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      leftNeighbour = B[0];
    }
    if (rank < numProcs - 1) {
      MPI_Sendrecv(&B[localSize - 1], 1, MPI_DOUBLE, rank + 1, 1, &rightNeighbour, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
      rightNeighbour = B[localSize - 1];
    }

    // Propagate the temperature
    for (int i = startIdx; i < endIdx; i++) {
      // Heat source remains constant
      if (i == source_x) {
        C[i - startIdx] = B[i - startIdx];
        continue;
      }

      // Current temperature
      value_t tc = B[i - startIdx];

      // Temperatures of adjacent cells
      value_t tl = (i == startIdx) ? leftNeighbour : B[i - 1 - startIdx];
      value_t tr = (i == endIdx - 1) ? rightNeighbour : B[i - startIdx + 1];

      // Compute new temperature
      C[i - startIdx] = tc + 0.2 * (tl + tr + (-2 * tc));
    }

    // Swap buffers
    Vector H = B;
    B = C;
    C = H;

    // Print temperature at intervals
    if (t % 10000 == 0) {
      MPI_Gatherv(B, localSize, MPI_DOUBLE, A, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (rank == 0) {
        printf("Step t=%d:\t", t);
        printTemperature(A, N);
        printf("\n");
      }
    }
  }

  // Gather final results from all processes
  MPI_Gatherv(B, localSize, MPI_DOUBLE, A, counts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  releaseVector(B);
  releaseVector(C);

  // End time measurement after computation
  endtime = MPI_Wtime();

  // ---------- check ----------
  int success = 1;

  if (rank == 0) {
    printf("Final:\t\t");
    printTemperature(A, N);
    printf("\n");

    for (int i = 0; i < N; i++) {
      value_t temp = A[i];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");

    // Print the time
    printf("WTime: %f seconds\n", endtime - starttime);

    // Free counts and displacements arrays
    free(counts);
    free(displs);
  }

  // ---------- cleanup ----------

  releaseVector(A);

  // Finalize MPI
  MPI_Finalize();
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) {
  free(m);
}

void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // Boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // Set the 'render' resolution
  int W = RESOLUTION;

  // Step size in each dimension
  int sW = N / W;

  // Room
  // Left wall
  printf("X");
  // Actual room
  for (int i = 0; i < W; i++) {
    // Get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // Pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // Print the average temperature
    printf("%c", colors[c]);
  }
  // Right wall
  printf("X");
}
