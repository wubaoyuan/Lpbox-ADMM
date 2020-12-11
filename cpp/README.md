# README

## Environment

This solver is built and passed all the tests on the machine with following specs

* CPU: Ryzen 5 2600
* Graphics Card: GTX 1070Ti
* Memory: 32GB
* Compiler: g++ 9.3.0
* Operating system: Ubuntu 20.04LTS
* Eigen library version: 3.3.8
* Opencv version: 4.4.0

## Build and installation

### Prerequisites

#### Install Eigen library

The whole program depends on the Eigen template library. You can find their main page here http://eigen.tuxfamily.org/index.php?title=Main_Page. 

For Linux and Mac users, after downloading the template library, copy the folder `Eigen` and `unsupported` under the decompressed directory to `/usr/local/include` to finish the installation. 

For Windows users, you can change the makefile attribute `EIGEN_INCLUDE_DIRECTORY` to the parent directory of the Eigen library.

#### Install opencv4

There is a shell script in the demo directory named `install_opencv.sh` that does all the commands below. If you encounter any problem during the running the script, please refer to the detailed commands below. 

If you want to build the image segmentation demo along with the solver shared library, you need to install opencv4 to proceed. First, you need make a directory named `opencv_build`, enter the directory, and clone opencv from the github repository using 

`git clone https://github.com/opencv/opencv.git`

`git clone https://github.com/opencv/opencv_contrib.git`

Then enter the opencv directory by

`cd opencv`

Then checkout to opencv version 4.4.0 using git

`git checkout 4.4.0`

Leave that directory and enter the contrib directory and checkout again

`cd ../opencv_contrib`

`git checkout 4.4.0`

Then enter the opencv directory again and create a directory named `build` and enter that directory by

`cd ../opencv`

`mkdir build`

`cd build`

Then use cmake to generate the makefile by the following command

`cmake -D CMAKE_BUILD_TYPE=RELEASE \
    -D CMAKE_INSTALL_PREFIX=/usr/local \
    -D INSTALL_C_EXAMPLES=ON \
    -D INSTALL_PYTHON_EXAMPLES=ON \
    -D OPENCV_GENERATE_PKGCONFIG=ON \
    -D OPENCV_EXTRA_MODULES_PATH=~/opencv_build/opencv_contrib/modules \
    -D BUILD_EXAMPLES=ON ..`

Notice that in this step, the CMake will download some dependencies. If you fail to download these dependencies (and you can check it by reading the output messages in the terminal), you will fail to build opencv. If you cannot download the required files, you can use the command 

`export https_proxy=${your_proxy_address_here}`

to let CMake download files using your proxy.

After the makefile is generated, use command

`make -j8`

to build opencv where the number after `-j` is the number of threads used in building. After it is successfully build and installed, type

`sudo ldconfig`

to update the symbols inside your operating system.

### Install the solver

If you want to install the solver only, you can type

`make solver && sudo make install && sudo ldconfig`

If you want to build all the demos, change `make solver` to `make all`

## TESTING

You can see the result of the image segmentation demo by entering the demo directory by

`cd demo`

and call the corresponding shell script by

`./segmentation.sh`

You can view the result of the clustering by calling the shell script

`./clustering.sh`

The output of the image segmentation is written to `cameraman.bmp` and the intermediate results are written to `cameraman_out.txt`. The intermediate results of the clustering task is written to `glass_out.txt`, `wine_out.txt` and `iris_out.txt`. For more specific input information, use

`./image_segmentation -h` and `./clustering -h`