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
        * [Default Models](#default-models)
        * [Grids](#grids)

## Requirements

### Libraries

#### Boost
When compiling both ALPSCore and MaxEnt, be careful to ensure boost was compiled with the same library and stdlib as ALPSCore and MaxEnt.   

We also require `boost/numeric/bindings`. 
If your boost package does not come with it, simply checkout from [here](https://svn.boost.org/svn/boost/sandbox/numeric_bindings/) and copy the `boost` folder wherever your installation directory is (`/usr/local/Cellar/` for homebrew or `/usr/local/include/`) 

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
  
