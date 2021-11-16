g++ -std=c++17 -Wpedantic -Wall -Wextra -DBE_VEGAS -O3 -c main.cpp
g++ main.o -lgsl -lgslcblas -lm -o main

