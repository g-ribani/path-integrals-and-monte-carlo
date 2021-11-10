g++ -std=c++17 -Wpedantic -Wall -Wextra -c main.cpp -DBE_VEGAS -DBE_CRUDE
g++ main.o -lgsl -lgslcblas -lm -o main.exe
