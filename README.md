Maxent
======
[![Build Status](https://travis-ci.org/CQMP/Maxent.svg?branch=master)](https://travis-ci.org/CQMP/Maxent)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.txt)

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
  * [License and citations](#license-and-citations)

## Requirements

### Libraries

#### Boost
When compiling both ALPSCore and Maxent, be careful to ensure boost was compiled with the same library and stdlib as ALPSCore and Maxent.   

#### ALPSCore
ALPSCore needs to be properly installed, see [ALPSCore library](https://github.com/ALPSCore/ALPSCore). ALPSCore provides the location of the Boost libraries.

#### Eigen3
For our linear algebra routines we use Eigen3 version >=3.3. CMake finds it through its `Eigen3Config.cmake`; if it is not in a standard location, add its prefix to `CMAKE_PREFIX_PATH` or set `-DEigen3_DIR=/path/to/share/eigen3/cmake`.

#### GSL
Maxent requires the GNU Scientific Library (GSL), which can be found [here](https://www.gnu.org/software/gsl/). The choice of BLAS library (the included CBLAS or an external ATLAS/BLAS/etc) does not matter here as the only the integration library is used. If it is not in a standard location, use `-DGSL_ROOT_DIR=/path/to/gsl/prefix`.

#### LAPACK (Optional)
Eigen3 has a good SVD routine, but can be very slow for a large kernel.
Some systems, like OS X or those with Intel MKL, have precompiled BLAS/LAPACK routines that can be faster and as accurate as Eigen3.
To turn on LAPACK support for the SVD, use `-DMAXENT_USE_LAPACK=ON`.


## Installation
Maxent needs CMake 3.22 or newer and a C++17 compiler. Boost must be the same
version that ALPSCore was built with; CMake checks this.
```
$ git clone https://github.com/CQMP/Maxent
$ cmake -S Maxent -B build -DALPSCore_DIR=/path/to/alpscore/share/ALPSCore \
        -DBoost_DIR=/path/to/lib/cmake/Boost-<version> \
        -DCMAKE_INSTALL_PREFIX=/path/to/install
$ cmake --build build -j 8
$ cmake --install build
```
The default build type is `Release`. Presets for common configurations are in
`CMakePresets.json` (`cmake --preset release|dev|asan`, pass the paths above
with `-D`). Use `CXX=g++` (or clang++, etc.) before the first `cmake` command to
choose a compiler.

Options:

| Option | Default | Meaning |
|---|---|---|
| `MAXENT_BUILD_TESTS` | ON | unit tests and the regression suite |
| `MAXENT_BUILD_UTILITIES` | ON | `kk` and `legendre_convert` |
| `MAXENT_BUILD_PADE` | OFF | `pade` (needs GMP; currently does not compile) |
| `MAXENT_USE_LAPACK` | OFF | LAPACK instead of Eigen for the SVD of the kernel |
| `MAXENT_WERROR` | OFF | treat warnings as errors |
| `MAXENT_REGRESSION_FULL` | OFF | also run the full-size regression cases (about 1 min) |
| `MAXENT_USE_SYSTEM_GTEST` | OFF | use an installed GoogleTest instead of downloading 1.18.0 |

### Tests
Once compiled, run `ctest --test-dir build` (or `ctest --preset <preset>` for
a preset build) to make sure everything works. This runs the unit tests (label
`unit`) and the regression suite (label `regression-fast`, needs Python 3 with
numpy and h5py); see [test/regression/README.md](test/regression/README.md).

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
Because Pade requires GMP, it does not build automatically. To include it in your build, add `-DMAXENT_BUILD_PADE=ON` to the `cmake` command. (Pade currently does not compile.)

## License and citations

Maxent is distributed under the terms in [LICENSE.txt](LICENSE.txt). See [CITATION.md](CITATION.md) for citation guidance.
