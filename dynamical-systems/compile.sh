g++ -std=c++17 -Wpedantic -Wall -Wextra -fopenmp -c main.cpp -DBE_VEGAS -DBE_CRUDE
g++ main.o -lgsl -lgslcblas -lm -fopenmp -o main.exe
