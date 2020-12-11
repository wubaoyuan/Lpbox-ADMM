EIGEN_LIB_LOCATION?=/usr/local/include/Eigen/
OPENBLAS_LIB_LOCATION?=/usr/lib/x86_64-linux-gnu/libopenblas.a
CFLAGS?=-O3
CXX = g++

all: cpp_algo demo

demo: demo.cpp
	$(CXX) $(CFLAGS) -o demo demo.cpp

cpp_algo: cpp_algo.cpp solver.h solver.cpp
	$(CXX) $(CFLAGS) -o cpp_algo -I $(EIGEN_LIB_LOCATION) cpp_algo.cpp solver.cpp solver.h -lpthread -lgfortran
