echo assuming Boost is installed at /usr/local/boost_1_77_0...
g++ correlators.cpp -std=c++17 -Wpedantic -Wall -Wextra -I /usr/local/boost_1_77_0 -I ../.. -lm -lgsl -lblas -O3 -o correlators
