Legendre Convert
======
A small utility for converting Green's functions in time space to Legendre space - an orthogonal polynomial basis.  
Based on [doi:10.1103/PhysRevB.84.075145](http://dx.doi.org/10.1103/PhysRevB.84.075145) which is also on 
[arxiv](http://arxiv.org/abs/1104.3215).

Table of Contents
=================
  * [Legendre Convert](#legendre-convert)
    * [Requirements](#requirements)
    * [Installation](#installation)
    * [Usage](#usage)
    * [Input Format](#input-format)
    * [Sample Procedure](#sample-procedure)
    * [Integration Trouble](#integration-trouble)
    
## Requirements
All that is required is Boost, which should be found by cmake
## Installation
If building Maxent, this will automatically be built alongside it.
However, to build just this project simply navigate to root, and use  `cmake CMakeLists.txt `  
If you have trouble finding the correct boost the following may help:
```
cmake CMakeLists.txt -DBOOST_ROOT=/path/to/boost \
             -DBOOST_NO_SYSTEM_PATHS=ON
```
## Usage
See `./legendre_convert --help` for a list of required and availble parameters. 
## Input Format
By default, the G(tau) data should be in a file `Gtau.dat` and in the following format:
```
tau0 Gtau0 sigma0
tau1 Gtau1 sigma1
...
```
*Currently only equispaced tau points are supported.*  
It is preferred to have an odd number of points including the end points, so that ndat-1 is even.
At least one end point is also required (tau=0 or tau=beta)

## Sample Procedure
Assuming there is a file in the same directory as legendre_convert labeled `Gtau.dat` with 201 equispaced tau points [0,beta=2] inclusive of Particle-Hole Symmetric Data   
A standard procedure is:
```
$ ./legendre_convert --beta 2 --plot_convergence 100 
```
This creates two files, `Gl.dat` and `Gl_convergence.dat`. Plotting Gl_convergence, it should be clear whata good lmax range is (most likely between 15-40)  
Once a good lmax is known (say, 25)
```
$ ./legendre_convert --beta 2 --lmax 25 --backcontinue
```
This should create 25 G(l) points, a back continued `Gtau_back.dat` file, along with the errors `Gtau_back_errors.dat`  
\> 90% of the points should be within errorbars, otherwise decrease/increase lmax until a maximum is found

## Integration Trouble
One major source of error is the integration of the conversion G(tau)->G(l). To confirm that your data is behaving properly,
the convergence of the tail and low frequency points should be well behaved with and without the tail fix.  
Example:
``` 
$ ./legendre_convert --beta 2 --plot_convergence 100
// generate plot
$ ./legendre_convert --beta 2 --plot_convergence 100 --notail
// generate plot
```
![converg_w_tail](https://cloud.githubusercontent.com/assets/7354063/8853103/22c00bbe-3127-11e5-87c7-3a7674cfa3e9.png)
![converg_wo_tail](https://cloud.githubusercontent.com/assets/7354063/8853105/287b10da-3127-11e5-9d32-f65dbfa6040a.png)

The high frequency points very quickly gain large errors relying on integration alone. On closer inspection:


![converg_wo_tail2](https://cloud.githubusercontent.com/assets/7354063/8853802/9377a8c8-312a-11e5-8f24-b0efd6f5040b.png)



For lmax between 10 and 20 the conversion procedure seems well behaved. 
