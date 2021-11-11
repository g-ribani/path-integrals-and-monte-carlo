g++ -std=c++17 -Wpedantic -Wall -Wextra -c main.cpp -DBE_CRUDE -DBE_VEGAS
g++ main.o -lgsl -lgslcblas -lm -o main
