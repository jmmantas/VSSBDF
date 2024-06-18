
# VSSBDF:  Variable Stepsize Semi-implicit Backward Differentiation Formula Solvers in C++
**VSSBDF** software package includes the C++ implementation of a family of Variable-Stepsize SBDF methods to solve several stiff IVP-ODEs derived from  Advection-Diffusion-Reaction Models. 

The IVP-ODE  is given by: 

dy/dt=F(t,y(t)) = f1(t,y(t)) + f2(t,y(t)) ,   y(t_0)=y_0,

being f1(t,y(t)) is the nonstiff component in F(t,y(t)) and f2(t,y(t)) is the stiff component.

### Two different types of splitting
The package allows two different types of additive splitting (to separate the stiff component and the nonstiff component in F(t,y(t)) ) to be used:

a) **Physical splitting**: It is the usual splitting, where the terms are split based on their physical properties, e.g. diffusion+reaction could treated implicitly while advection is treated explicitly: dy/dt= f1(t,y(t)) + f2(t,y(t)) = [Adv(t,y)] + [Diff(t,y) + React(t,y)]  


b) **Jacobian Splitting**:  It reflects the numerical properties of the solution. The splitting has not to be provided by the user because it is performed implicitly by the solver.
A splitting Jacobian of a ODE system is written as:     dy/dt= f1(t,y(t)) + f2(t,y(t)) = [J_F*y(t)]  +  [F(t,y(t))-F_F*y(t)]



## Installation

Before compiling and link the software, you need to install the libraries LIS (https://www.ssisc.org/lis/index.en.html), and openblas (https://www.openblas.net/).

It is necessary to edit Makefile to indicate the path where the LIS library is installed by setting the variable LISROOT which is set by default to "usr/local".

Once the Makefile has been modified, you have to compile and link by issuing:

make

This will generate the executable file "IVP_Solver". 

## Usage

To execute the program, you have to issue the executable file with the following  arguments:

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

***IVP_ODE-advdiff.h:  it includes the definition of the class IVP_ODE and the implementation of its methods. These methods define the particularities of the ODEs system (1D Advection-Diffusion problem) including:
  (a) The calculation of Y0 (the initial stage vector).  
  (b) The function F (for the non-stiff part) and G (for the stiff part).
  (c) The function to calculate analytically the Jacobian d and G at (t,Y).

*** IVP_Solver.cpp: Declares an object of the class IVP_ODE to implement the  k-th order Runge-Kutta  method (RK function with k=1,2,3 4) and the  semiimplicit SBDF-order method (SBDF function). The RK function is used to approximate the initial step value using a very small stepsize and  to check the SBDF solution.  In the main function of this file, we run the RK and SBDF-order methods to solve a particular Problem with the parameters (number of grid points, final time, stepsize for each method) provided by the user using the command line.

## License

The code in this repository is released under GPLv3, 2024 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.