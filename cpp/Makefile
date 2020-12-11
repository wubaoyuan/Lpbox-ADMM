ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
VPATH=./src:./demo
EIGEN_LIB_LOCATION?=/usr/local/include
OPENCV_LIB_LOCATION?=/usr/local/include/opencv4
CLUSTERING_SRC=demo/clustering.cpp src/read_matrix_utils.cpp
IMAGE_SRC=demo/image_segmentation.cpp src/image_segmentation_utils.cpp

SRC_DIR=./src
DEMO_DIR=./demo
LIB_SOLVER_SRC=src/LPboxADMMsolver.cpp
CFLAGS?=-O3
CXX = g++
ifeq ($(PREFIX),)
	PREFIX := /usr/local
endif

all: solver image_segmentation clustering

release: CFLAGS=-O3
release: solver

install: liblpbox_solver.so src/LPboxADMMsolver.h
	install -d $(DESTDIR)$(PREFIX)/lib/
	install -m 644 liblpbox_solver.so $(DESTDIR)$(PREFIX)/lib/
	install -d $(DESTDIR)$(PREFIX)/include/
	install -m 644 src/LPboxADMMsolver.h $(DESTDIR)$(PREFIX)/include/

uninstall:
	rm $(DESTDIR)$(PREFIX)/lib/liblpbox_solver.so 
	rm $(DESTDIR)$(PREFIX)/include/LPboxADMMsolver.h 

clustering: clustering.cpp read_matrix_utils.cpp read_matrix_utils.h LPboxADMMsolver.h LPboxADMMsolver.cpp
	$(CXX) $(CFLAGS) -o demo/clustering -I $(EIGEN_LIB_LOCATION) -I $(ROOT_DIR) -I $(SRC_DIR) -I $(DEMO_DIR) $(CLUSTERING_SRC) liblpbox_solver.so -lpthread

image_segmentation: image_segmentation.cpp image_segmentation_utils.cpp image_segmentation_utils.h LPboxADMMsolver.h LPboxADMMsolver.cpp
	$(CXX) $(CFLAGS) -o demo/image_segmentation -I $(EIGEN_LIB_LOCATION) -I $(OPENCV_LIB_LOCATION) -I $(SRC_DIR) -I $(DEMO_DIR) $(IMAGE_SRC) liblpbox_solver.so -lpthread $(shell pkg-config --libs opencv4) 

solver: LPboxADMMsolver.h LPboxADMMsolver.cpp
	$(CXX) $(CFLAGS) -fPIC -shared $(LIB_SOLVER_SRC) -o liblpbox_solver.so

clean:
	rm -f liblpbox_solver.so demo/image_segmentation_demo demo/clustering_demo
