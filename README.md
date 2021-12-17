# ECE1747-Parallel-N-body-simulation

To build the code, use the following command:

for pthread version:

g++ -o test -fopenmp -lpthread test.cpp -O3

for openmp version:

g++ -o test -fopenmp -lpthread test_benchmark_20211125.cpp -O3
