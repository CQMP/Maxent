MaxEnt
======

MaxEnt Project

This branch requires the [ALPSCore library](https://github.com/ALPSCore/ALPSCore). 

Table of Contents
=================
  * [MaxEnt](#maxent)
    * [Requirements](#requirements)
      * [Libraries](#libraries)
        * [Boost](#boost)
        * [BLAS/LAPACK](#blaslapack)
        * [ALPSCore](#alpscore)
    * [Installation](#installation)
      * [Tests](#tests)
    * [Usage](#usage)
        * [Input](#input) 
        * [Kernels](#kernels)
        * [Default Models](#default-models)
        * [Grids](#grids)

## Requirements

### Libraries

#### Boost
When compiling both ALPSCore and MaxEnt, be careful to ensure boost was compiled with the same library and stdlib as ALPSCore and MaxEnt.   

We also require `boost/numeric/bindings`. Cmake will automatically try to find it or fetch it using subversion.

#### BLAS/LAPACK
Thanks to Boost this version currently uses uBLAS, but requires LAPACK support. Modify the `CMakeLists.txt` file as needed

#### ALPSCore
Please make sure to install ALPSCore with the cmake flags `-DALPS_HAVE_BLAS=1 -DALPS_HAVE_LAPACK=1`. This will ensure the actual libraries are used. If they are not set, internal ublas will be used instead.

## Installation
To install provide something like:
```
$ git clone https://github.com/CQMP/Maxent  
$ mkdir build  
$ cd build  
$ cmake  ../ -DCMAKE_INSTALL_PREFIX=/path/to/here/Maxent/build -DALPSCore_DIR=/path/to/alpscore/build/share/ALPSCore
$ make -j 8
```
Sometimes it is more convenient to have `CC=gcc CXX=g++` (or clang, etc) before the cmake command.

### Tests
Once compiled please run `./maxent --test`
to ensure everything works.

## Usage
See `./maxent --help` for a list of required and availble parameters. 

#### Input
##### Particle Hole Symmetric Data
The Green's function for PH symmetric data is 0, therefore we only require the imaginary part.
Input file:
```
n imag sigma
//example:
0 imag0 sigma0
1 imag1 sigma1
...
```
##### Non-Particle Hole Symmetric Data
This assumes a non-zero real part of the Green's function. Input data should be:
```
n real sigma_real
n imag sigma_imag
//example:
0 real0 sigma_real0
0 imag0 sigma_imag0
1 real1 sigma_real1
1 imag1 sigma_imag1
```
#### Kernels
![screen shot 2015-07-17 at 3 48 00 pm](https://cloud.githubusercontent.com/assets/7354063/8755755/47c93aba-2c9b-11e5-8743-359ab6271827.png)
![screen shot 2015-07-17 at 3 48 04 pm](https://cloud.githubusercontent.com/assets/7354063/8755770/57c4ab3e-2c9b-11e5-98a3-1a073d67ee34.png)

#### Default Models
* Flat
* Gaussian
  * Shifted Gaussian
  * Double Gaussian
  * Two Gaussians
  * Double General Gaussian
* Lorentzian
  * See Gaussian
* Linear Rise Exponential Decay
* Quadratic Rise Exponential Decay
* Tabulated Default model = "Filename"

#### Grids
Maxent creats a default model on a grid between [0,1]

![grids](https://cloud.githubusercontent.com/assets/7354063/8681331/cb2b0852-2a34-11e5-9485-08c8c6a68274.png)
  
