echo assuming Boost directory is /usr/local/boost_1_77_0...
g++ anharmonic.cpp -std=c++17 -Wpedantic -Wall -Wextra -I/usr/local/boost_1_77_0 -lm -lgsl -lblas -O3 -o anharmonic