CC=mpicc
CFLAGS=-O2 -std=gnu11 -Wall -Wextra -pedantic

.PHONY: all
all: heat_stencil_1D_mpi

heat_stencil_1D_mpi: heat_stencil_1D_mpi.c
	@$(CC) $(CFLAGS) $< -o heat_stencil_1D_mpi

.PHONY: clean
clean:
	@rm heat_stencil_1D_mpi

.PHONY: run
run: all
	@echo "Parallel:"
	@./heat_stencil_1D_mpi
	@echo