#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string>

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

  // Calculate chunk size, start and end point for each process
  int chunkSize = N / numProcs;
  int remainder = N % numProcs;

  if (rank == 0) {
    if (remainder > 0) {
      printf("The problem size can't be evenly distributed to the number of processes. Please redefine it.");
      return 0;
    }
  }

  int startIdx = rank * chunkSize;
  int endIdx = startIdx + chunkSize;
  int source_x = N / 4;

  // create a buffer for storing temperature fields
  Vector A = createVector(N);

  if (rank == 0) {
    
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

    // ---------- setup --------- 
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
      A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    A[source_x] = 273 + 60;

    printf("Initial:\t");
    printTemperature(A, N);
    printf("\n");
  } 

  // Distribute A to all processes
  MPI_Bcast(A, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // ---------- compute ----------

  // create a second buffer for the computation and initialize with A
  Vector B = createVector(chunkSize);
  std::copy(A + startIdx, A + endIdx, B);
  // create variables for the neighbours
  value_t rightNeighbour;
  value_t leftNeighbour;

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // send and receive boundary elements
    if (rank > 0) {
      MPI_Send(&B[chunkSize - 1], 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < numProcs - 1) {
      MPI_Send(&B[0], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (rank > 0) {
      MPI_Recv(&leftNeighbour, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < numProcs - 1) {
      MPI_Recv(&rightNeighbour, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    // .. we propagate the temperature
    for (long long i = startIdx; i < endIdx; i++) {
      // center stays constant (the heat is still on)
      if (i == source_x) {
        B[i - startIdx] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : tc;
      value_t tr = (i != N - 1) ? A[i + 1] : tc;

      // compute new temperature at current position
      B[i - startIdx] = tc + 0.2 * (tl + tr + (-2 * tc));
    }

    // MPI_Gather(B, chunkSize, MPI_DOUBLE, A, chunkSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      // show intermediate step
      if (!(t % 10000)) {
        MPI_Gather(B, chunkSize, MPI_DOUBLE, A, chunkSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        printf("Step t=%d:\t", t);
        printTemperature(A, N);
        printf("\n");
      }
    }
  }

  releaseVector(B);

  // ---------- check ----------
  int success = 1;

  if (rank == 0) {
    printf("Final:\t\t");
    printTemperature(A, N);
    printf("\n");

    for (long long i = 0; i < N; i++) {
      value_t temp = A[i];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  }

  // ---------- cleanup ----------

  releaseVector(A);

  // done
  MPI_Finalize();
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}