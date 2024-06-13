
# VSSBDF:  Variable Stepsize, Semi-implicit, Backward Differentiation Formula Solvers in C++
**VSSBDF** software includes the C++ implementation of Variable-Stepsize IMEX SBDF methods to solve several stiff IVP-ODEs derived from  Advection-Diffusion-Reaction Models



## Installation

Before compiling and link the software, you need to install the libraries LIS (https://www.ssisc.org/lis/index.en.html), and openblas (https://www.openblas.net/).

It is necessary to edit Makefile to indicate the path where the LIS library is installed by setting the variable LISROOT which is set by default to "usr/local".

Once the Makefile has been modified, you have to compile and link  by issuing:

make

This will generate the executable file "IVP_Solver". To execute the program, you have to issue the executable file with the following  arguments:

## Usage

./IVP_Solver <problem id.> <conv. order> <num. grid points (N))>  <tf(final time)> <RK stepsize> <init. SBDF stepsize>  <Num. experiments> <init_tolerance>

problem.id= 0: 1D Advection-Diffusion
problem.id= 1: Combustion
problem.id= 2: 1D Brusselator
problem.id= 3: 2D Heat using a grid with Nx(N/7) points
problem.id= 4: CUSP
problem.id= 5: Simple 1D Advection-Diffusion
problem.id= 6: 2D Brusselator using a grid with NxN points
problem.id= 7: Stiff version of the 1D Brusselator

By default all 1D problems use a grid of N points.

***IVP_ODE-advdiff.h:  it includes the definition of the class IVP_ODE and the implementation of its methods. These methods defines the particularities of the ODEs system (1D Advection-Diffusion problem) including:
  (a) The calculation of Y0 (the initial stage vector).  
  (b) The function F (for the non-stiff part) and G (for the rigid part).
  (c) The function to calculate analytically the Jacobian d and G at (t,Y).

*** IVP_Solver.cpp: Declares an object of the class IVP_ODE to implement the  k-th order Runge-Kutta  method (RK function with k=1,2,3 4) and the  semiimplicit SBDF-order method (SBDF function). The RK function is used to approximate the initial step value using a very small stepsize and  to check the SBDF solution.  In the main function of this file, we run the RK and SBDF-order methods to solve a particular Problem with the parameters (number of grid points, final time, stepsize for each method) provided by the user using the command line.
