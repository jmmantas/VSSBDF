/*This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_stiff_brusselator" which represents  
a stiff variation of the one-dimensional Brusselator model.

VSSBDF Copyright (C) 2024 Jose Miguel Mantas Ruiz (jmmantas@ugr.es) and Raed Ali Mara'Beh (raedmaraabeh@gmail.com)

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
*/




#ifndef STIFF_BRUSSELATOR_H
#define STIFF_BRUSSELATOR_H


#include "IVP_ODE.h"

using namespace std;


//*****************************************************************
// Class for the IVP-ODE  which represents 
// the Brusselator 2D model 
//*****************************************************************

class IVP_ODE_stiff_brusselator : public IVP_ODE {
private:
    int nx; // Number of grid points in the spatial domain
    double dtx; // Spatial step size
    double dtx2; //dtx*dtx
    double dtx_2; //2*dtx
    double DD; //alpha/dtx2
    double AA; //rho/dtx_2
    double alpha=0.01; // Diffusion coefficients for u, v, and w

    //*********************************************************
    // IN THE PAPER rho=0.001!!!!
    double rho=0.01; // Advection coefficients for u, v, and w
    //*********************************************************

    double a = 0.6;  // Reaction parameters
    double b = 2.0; // Reaction parameters
    double eps=0.01; // Parameter to control stiffness of the reaction term for w
    double pi = 3.14159265358979323846;
 

    // Indexation function which maps 2D spatial coordinates (i,j)
    // to a 1D position in a vector
    inline int idx(const int i, const int j) { return i * 3 + j; }


   


    //******************************************************************
    // Set the elements in close to main diagonal in the Jacobian matrix
    //******************************************************************
    void set_tridiag(const int i,
        const double u, const double v, const double w, const double Aii,
        LIS_INT* index, LIS_INT* k_, LIS_SCALAR* value);


public:
    //******************************************************
    // Constructor of the class IVP_ODE_cusp
    //******************************************************
    IVP_ODE_stiff_brusselator(const int nx_points);

    //******************************************************
    // Initialize stage vector Y0 with neqn components
    //******************************************************
    inline void init(double* Y0);

    //******************************************************
    //vector system function for the nonstiff term DY=G(t,Y)
    //******************************************************
    void G(const double t, const double* Y, double* DY);

    //******************************************************
    //vector system function for the nonstiff term DY=F(t,Y)
    //******************************************************
    void F(const double t, const double* Y, double* DY);

    //******************************************************
    //ODE vector system function feval: DY=feval(t,Y)
    //******************************************************
    void feval(const double t, const double* Y, double* DY);



    //******************************************************
    //Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
    // in the IVP_ODE
    //****************************************************** 
    void Compute_G_Matrix_Exact(const double t, double a, double* Y, LIS_MATRIX As);


    //******************************************************
    //Compute matrix Jf= Jacobian of the function 
    // feval defining the IVP_ODE
    //******************************************************       
    void Compute_Feval_Jacobian_exact(const double t, double* Y, LIS_MATRIX As);

    //***************************************************
    // Function to obtain the exact solution Y=U(t,X)
    //***************************************************
    void exact_solution(const double t, double* Y);
};

//*****************************************************
// Constructor of the class IVP_ODE_brusselator2D
//*****************************************************
IVP_ODE_stiff_brusselator::IVP_ODE_stiff_brusselator(const int nx_points) {
     IVP_name = "stiff_Brusselator";
    nx = nx_points;
    // Number of ODEs
    neqn = 3 * nx;

    // Number of non-zero elements in the Jacobian matrix
    nnz_G = 9 * nx-12;
    nnz_FEVAL = 14 * nx - 12;
    dtx = 1.0 / (nx - 1); // Adjusted for boundary conditions
    dtx2 = dtx * dtx;
    dtx_2 = 2*dtx ;
    DD = alpha / dtx2;
    AA = -rho / dtx_2;
}

    // Initialize stage vector Y0 with neqn components
//******************************************************
    inline void IVP_ODE_stiff_brusselator::init(double *Y0)
    {
      for (int i = 0; i < nx; i++) {
        double x = i * dtx;
        Y0[idx(i, 0)] = a + 0.1 * sin(pi * x); // Initial condition for u
        Y0[idx(i, 1)] = b / a + 0.1 * sin(pi * x); // Initial condition for v
        Y0[idx(i, 2)] = b + 0.1 * sin(pi * x); // Initial condition for w
      }
    }

    inline void IVP_ODE_stiff_brusselator::feval(const double t, const double* Y, double* DY)
    {
        double Ytmp[neqn];
        G(t, Y, DY);
        F(t, Y, Ytmp);
        cblas_daxpy(neqn, 1.0, Ytmp, 1, DY, 1);
    }
//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//***************************************************
void IVP_ODE_stiff_brusselator::G(const double t, const double *Y, double *DY) 
{
   // First, set the boundary conditions to zero derivative (stationary)
   DY[idx(0, 0)] = 0; // u at the left boundary
   DY[idx(nx - 1, 0)] = 0; // u at the right boundary
   DY[idx(0, 1)] = 0; // v at the left boundary
   DY[idx(nx - 1, 1)] = 0; // v at the right boundary
   DY[idx(0, 2)] = 0; // w at the left boundary
   DY[idx(nx - 1, 2)] = 0; // w at the right boundary

   // Apply the diffusion term for interior points
   for (int i = 1; i < nx - 1; i++) {
    // Indexes for u, v, and w at grid point i
    const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
    // Apply the diffusion term using the Laplacian discretization for interior points
    DY[i0] = DD * (Y[idx(i - 1, 0)] - 2.0 * Y[i0] + Y[idx(i + 1, 0)]);
    DY[i1] = DD * (Y[idx(i - 1, 1)] - 2.0 * Y[i1] + Y[idx(i + 1, 1)]);
    DY[i2] = DD * (Y[idx(i - 1, 2)] - 2.0 * Y[i2] + Y[idx(i + 1, 2)]);
   } 

}

//*******************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//*******************************************************
void IVP_ODE_stiff_brusselator::F(const double t, const double *Y, double *DY) {
    // Apply stationary boundary conditions at the first and last points explicitly for advection terms
  // Assuming no change due to advection at the boundaries
    DY[idx(0, 0)] = 0;
    DY[idx(0, 1)] = 0;
    DY[idx(0, 2)] = 0;
    DY[idx(nx - 1, 0)] = 0;
    DY[idx(nx - 1, 1)] = 0;
    DY[idx(nx - 1, 2)] = 0;

    for (int i = 1; i < nx - 1; i++) {
      const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
      const double u = Y[i0];
      const double v = Y[i1];
      const double w = Y[i2];

      // Calculate central differences for advection terms
      double du_dx = AA*(Y[idx(i + 1, 0)] - Y[idx(i - 1, 0)]) ;
      double dv_dx = AA*(Y[idx(i + 1, 1)] - Y[idx(i - 1, 1)]) ;
      double dw_dx = AA*(Y[idx(i + 1, 2)] - Y[idx(i - 1, 2)]) ;

      // Reaction terms
      DY[i0] = du_dx+a - (w + 1) * u + u * u * v ; // Adding advection term for u
      DY[i1] = dv_dx + w * u - u * u * v ; // Adding advection term for v
      DY[i2] = dw_dx + (b - w) / eps - w * u; // Adding advection term for w
    }

   //  For boundary points, calculate the reaction terms without advection contributions
    // This maintains the stationary boundary condition by not allowing flux due to advection
    // at the boundaries but still accounts for reactions happening at all points
    
    for (int i : {0, nx - 1}) {
   
        const int i0 = idx(i, 0), i1 = i0 + 1, i2 = i0 + 2;
        const double u = Y[i0];
        const double v = Y[i1];
        const double w = Y[i2];

        // Only reaction terms are considered at the boundaries
        DY[i0] = a - (w + 1) * u + u * u * v; // u reaction term
        DY[i1] = w * u - u * u * v; // v reaction term
        DY[i2] = (b - w) / eps - w * u; // w reaction term
    }
}


//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the stiff Brusselator model
//******************************************************       
void IVP_ODE_stiff_brusselator::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX As)
{
  LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  
  const double Aii = 1 + 2.0 * A * DD; // Coefficient for the diagonal elements
  const double other = -A * DD; // Coefficient for the off-diagonal elements

  init_row_insertion(ptr, k, row);

  // row  (0,j) 
  for (int j = 0; j <= 2; j++) {
    new_entry(index, value, k, idx(0, j), 1.0);
    next_row(ptr, k, row);
  }

  // row(i,j), i=1,...,nx-2
  for (int i = 1; i <= nx - 2; i++) {
    for (int j = 0; j <= 2; j++) {
      new_entry(index, value, k, idx(i-1, j), other);
      new_entry(index, value, k, idx(i, j)  , Aii);
      new_entry(index, value, k, idx(i+1, j), other);
      next_row(ptr, k, row);
    }
  }
  // row (nx-1,j)
  for (int j = 0; j <= 2; j++) {
      new_entry(index, value, k, idx(nx-1,j), 1.0);
      next_row(ptr, k, row);
  }
 
}



//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_stiff_brusselator::Compute_Feval_Jacobian_exact(const double t, double *Y, LIS_MATRIX As)
{
  LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  
  // Coefficient for the diagonal elements
  const double Aii = - 2.0 * DD; 

  const double otherm1 = DD-AA; 
  const double otherp1 = DD+AA;   

  init_row_insertion(ptr, k, row);

  int i;
  double u,v,w;
  
  // row  (0,0) 
  i=0; u = Y[idx(i, 0)];v = Y[idx(i, 1)];w = Y[idx(i, 2)];
  new_entry(index, value, k, idx(i, 0), - (w + 1)+2.0*u*v);
  new_entry(index, value, k, idx(i, 1),  u*u );
  new_entry(index, value, k, idx(i, 2), -u    );
  next_row(ptr, k, row);
  // row  (0,1) 
  new_entry(index, value, k, idx(i, 0) , w-2.0*u*v );
  new_entry(index, value, k, idx(i, 1), -u*u);
  new_entry(index, value, k, idx(i, 2), u);
  next_row(ptr, k, row);
  // row  (0,2)
  new_entry(index, value, k, idx(i, 0), -w);
  new_entry(index, value, k, idx(i, 2), -1.0/eps-u);
  next_row(ptr, k, row);
  

  // row(i,j), i=1,...,nx-2
  for ( i = 1; i <= nx - 2; i++) {
    u = Y[idx(i, 0)]; v = Y[idx(i, 1)]; w = Y[idx(i, 2)];
    // row  (i,0)
    new_entry(index, value, k, idx(i-1, 0), otherm1);
    new_entry(index, value, k, idx(i, 0), Aii- (w + 1)+2.0*u*v);
    new_entry(index, value, k, idx(i, 1),  u*u );
    new_entry(index, value, k, idx(i, 2), -u    );
    new_entry(index, value, k, idx(i+1, 0), otherp1);
    next_row(ptr, k, row);
    // row  (i,1)
    new_entry(index, value, k, idx(i-1, 1), otherm1); 
    new_entry(index, value, k, idx(i, 0) , w-2.0*u*v );
    new_entry(index, value, k, idx(i, 1), Aii-u*u);
    new_entry(index, value, k, idx(i, 2), u);
    new_entry(index, value, k, idx(i+1, 1), otherp1);
    next_row(ptr, k, row);
    // row  (i,2)
    new_entry(index, value, k, idx(i-1, 2), otherm1);
    new_entry(index, value, k, idx(i, 0), -w);
    new_entry(index, value, k, idx(i, 2), Aii-1.0/eps-u);
    new_entry(index, value, k, idx(i+1, 2), otherp1);
    next_row(ptr, k, row);    
  }


  // row (nx-1,0)
  i=nx-1; u = Y[idx(i, 0)];v = Y[idx(i, 1)];w = Y[idx(i, 2)];
  new_entry(index, value, k, idx(i, 0), - (w + 1)+2.0*u*v);
  new_entry(index, value, k, idx(i, 1),  u*u );
  new_entry(index, value, k, idx(i, 2), -u    );
  next_row(ptr, k, row);
  // row  (nx-1,1) 
  new_entry(index, value, k, idx(i, 0) , w-2.0*u*v );
  new_entry(index, value, k, idx(i, 1), -u*u);
  new_entry(index, value, k, idx(i, 2), u);
  next_row(ptr, k, row);
  // row  (nx-1,2)
  new_entry(index, value, k, idx(i, 0), -w);
  new_entry(index, value, k, idx(i, 2), -1.0/eps-u);
  next_row(ptr, k, row);
 
}







//***************************************************
// Function to obtain the exact solution Y=U(t,X)
//***************************************************
void IVP_ODE_stiff_brusselator::exact_solution(const double t, double *Y) {
    // Implementation of exact_solution

}

#endif