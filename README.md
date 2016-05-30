Maxent
======

The Maxent Project: A utility for performing analytic continuation using the method of Maximum Entropy.

Many-body Green's functions calculated on the imaginary axis can be related to a real spectral function, but is an ill-posed problem. One algorithm to solve for the spectral function is the maximum entropy method. This code is an implementation of maximum entropy method as well as useful utilities for dealing with Green's functions and analytic continuation. 

Table of Contents
=================
  * [Maxent](#maxent)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
      * [Libraries](#libraries)
        * [Boost](#boost)
        * [ALPSCore](#alpscore)
        * [Eigen3](#eigen3)
        * [GSL](#gsl)
        * [LAPACK (Optional)](#lapack-optional)
    * [Installation](#installation)
      * [Tests](#tests)
    * [Convention](#convention)
    * [Usage](#usage)
        * [Input](#input)
          * [Particle Hole Symmetric Data](#particle-hole-symmetric-data)
          * [Non-Particle Hole Symmetric Data](#non-particle-hole-symmetric-data)
          * [Time Data](#time-data)
        * [Kernels](#kernels)
        * [Default Models](#default-models)
        * [Grids](#grids)
  * [Utilities](#utilities)
    * [Pade](#pade)
    * [Kramers-Kronig](#kramers-kronig)
    * [Legendre Convert](#legendre-convert)

## Requirements

### Libraries

#### Boost
When compiling both ALPSCore and Maxent, be careful to ensure boost was compiled with the same library and stdlib as ALPSCore and Maxent.   

#### ALPSCore
ALPSCore needs to be properly installed, see [ALPSCore library](https://github.com/ALPSCore/ALPSCore). ALPSCore provides the location of the Boost libraries.

#### Eigen3
For our linear algebra routines we use Eigen3 version >=3.1. If not in your path use `-DEIGEN3_INCLUDE_DIR=`

#### GSL
Maxent requires the GNU Scientific Library (GSL), which can be found [here](https://www.gnu.org/software/gsl/). The choice of BLAS library (the included CBLAS or an external ATLAS/BLAS/etc) does not matter here as the only the integration library is used. If not in your path use `-GSL_ROOT_DIR=` to the path that has `bin/gsl-config`.

#### LAPACK (Optional)
Eigen3 has a good SVD routine, but can be very slow for a large kernel.
Some systems, like OS X or those with Intel MKL, have precompiled BLAS/LAPACK routines that can be faster and as accurate as Eigen3.
To turn on LAPACK support for the SVD, please use the flag `-DUSE_LAPACK=1`. 


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
Once compiled please run `make test`
to ensure everything works.  

## Convention
The Maxent project uses the following conventions:

![convention](https://cloud.githubusercontent.com/assets/7354063/10086355/ef8c8362-62db-11e5-938a-1c24139c72df.png)

![convention_gtau](https://cloud.githubusercontent.com/assets/7354063/10086425/570a68ce-62dc-11e5-8cd3-1e871f89c695.png) 

![convention_A_omega](https://cloud.githubusercontent.com/assets/7354063/10056184/0ce6afd4-6208-11e5-9bdd-556ae958857c.png)

To see more, see [this pdf](examples/conventions_and_kernels.pdf).

## Usage
Upon installation there will be a binary `maxent`. It uses [ALPSCore parameters](https://github.com/ALPSCore/ALPSCore/wiki/Tutorial%3A-parameters), which can take a param file as input (see [examples](./examples)), or command line arguments of the form `--PARAMETER=value`.  
The three required parameters are `BETA` (inverse temperature), `NDAT` (the number of input data points), and either the location of the input data `DATA` or input through the param file using `X_i` (see below).

See `./maxent --help` for a list of required and availble parameters. 

#### Input
##### Particle Hole Symmetric Data
The Green's function for PH symmetric data is 0, therefore we only require the imaginary part.
Input file:
```
omega_n imag sigma
//example:
omega_0  imag0 sigma0
omega_1 imag1 sigma1
...
```
Data can also be stored in the parameter file using:
```
X_0= xxxx
SIGMA_0=xxx
X_1=xxxx
SIGMA_1=xxx
...
X_ndat-1=xxxx
SIGMA_ndat-1=xxx
```
##### Non-Particle Hole Symmetric Data
This assumes a non-zero real part of the Green's function. Input data should be:
```
n real sigma_real
n+1 imag sigma_imag
//example:
omega_0 real0 sigma_real0 imag0 sigma_imag0
omega_1 real1 sigma_real1 imag1 sigma_imag1
```
**_NOTE:_** NDAT=#of points*2 when there is not Particle Hole Symmetric data  

Data can also be stored in the parameter file using:
```
X_0= xxxx
SIGMA_0=xxx
X_1=xxxx
SIGMA_1=xxx
...
```
where `X_0` is the real part and `X_1` is the imaginary part, etc.
##### Time Data
For either symmetric or non-symmetric data, G(tau) is simply input as:
```
tau_n Gtau_n sigma_n
tau_n+1 Gtau_n+1 sigma_n+1
```
You can also include tau points in the parameter file, defined like:
 ```
 TAU_0=xxx...
 TAU_1=xxx
 ...
 TAU_NDAT-1=xxx
 ```
#### Kernels
![Fermionic Kernels](https://cloud.githubusercontent.com/assets/7354063/10101709/42e4cae2-6368-11e5-999b-0483d4f4358f.png)
![Time Kernels](https://cloud.githubusercontent.com/assets/7354063/15372450/754a55fa-1d0e-11e6-8483-e2c827591946.png)

For the `TZero` kernel, supply any `BETA` value.

#### Default Models
[View Examples Here](examples/default_models.pdf)
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

![grids](https://cloud.githubusercontent.com/assets/7354063/14571315/8ac93a8e-0316-11e6-8255-b9756a2710e8.png)
  

# Utilities
## Kramers-Kronig
Requires: [GSL](http://www.gnu.org/software/gsl/), Boost
## Legendre Convert
Requires: Boost
## Optional
### Pade
Requires: [GMP](https://gmplib.org/),[Eigen3.1](http://eigen.tuxfamily.org/index.php?title=Main_Page)
Because Pade requires GMP, it does not build automatically. To include it in your build, either run `cmake` from the pade folder, or in your `maxent` build folder add `-DPADE=1` to the`cmake` command 
