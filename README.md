MaxEnt
======

MaxEnt Project

This branch requires the [ALPSCore library](https://github.com/ALPSCore/ALPSCore). 

Table of Contents
=================
  * [MaxEnt](#maxent)
  * [Table of Contents](#table-of-contents)
    * [Requirements](#requirements)
      * [Libraries](#libraries)
        * [Boost](#boost)
        * [Eigen3](#Eigen3)
        * [ALPSCore](#alpscore)
    * [Installation](#installation)
      * [Tests](#tests)
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
When compiling both ALPSCore and MaxEnt, be careful to ensure boost was compiled with the same library and stdlib as ALPSCore and MaxEnt.   

#### Eigen3
For our linear algebra routines we use Eigen3 version >=3.1. If not in your path use `-DEIGEN3_INCLUDE_DIR=`

#### ALPSCore
ALPSCore needs to be properly installed, see [ALPSCore library](https://github.com/ALPSCore/ALPSCore).

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
n+1 imag sigma_imag
//example:
0 real0 sigma_real0
1 imag0 sigma_imag0
2 real1 sigma_real1
3 imag1 sigma_imag1
```
NOTE: NDAT=#of points*2 when there is not Particle Hole Symmetric data
##### Time Data
For either symmetric or non-symmetric data, G(tau) is simply imput as:
```
n tau_n sigma_n
n+1 tau_n+1 sigma_n+1
```
This assumes an equispaced grid of tau points constructed as the following:
```
 tau[i]= i*beta/ (ndat_-1)
 ```
This produces a grid between [0,beta] inclusive. If your tau points do not follow such a pattern, then in the parameter file they should be defined:
 ```
 TAU_0=xxx...
 TAU_1=xxx
 ...
 TAU_NDAT-1=xxx
 ```
#### Kernels
![screen shot 2015-07-17 at 3 48 00 pm](https://cloud.githubusercontent.com/assets/7354063/8755755/47c93aba-2c9b-11e5-8743-359ab6271827.png)
![screen shot 2015-07-17 at 3 48 04 pm](https://cloud.githubusercontent.com/assets/7354063/8755770/57c4ab3e-2c9b-11e5-98a3-1a073d67ee34.png)

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

![grids](https://cloud.githubusercontent.com/assets/7354063/8681331/cb2b0852-2a34-11e5-9485-08c8c6a68274.png)
  

# Utilities
## Pade
Requires: [GMP](https://gmplib.org/),[Eigen3.1](http://eigen.tuxfamily.org/index.php?title=Main_Page)
To point cmake to the correct eigen directory, use `-DEIGEN3_INCLUDE_DIR=/path/to/eigen`
(This directory should be the root of the downloaded tar.bz2 file) 
## Kramers-Kronig
Requires: [GSL](http://www.gnu.org/software/gsl/), Boost
## Legendre Convert
Requires: Boost
