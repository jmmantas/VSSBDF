/*This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_brusselator" which represents  
a one-dimensional Brusselator model.

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



#ifndef BRUSSELATOR_H
#define BRUSSELATOR_H

#include "IVP_ODE.h"

using namespace std;

//*****************************************************************
// Class for the IVP-ODE representing the 1D Brusselator model 
//*****************************************************************
class IVP_ODE_brusselator:public IVP_ODE{

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  double DD;
  const double alpha=1.0/50;
  const double  A=1.0;   
  const double  B=3.0;   
  const double pi=3.14159265358979;

  double f(const double y);
 
  // Indexation function which maps 2D spatial coordinates (i,j)
  // to a 1D position in a vector
  //inline int idx(const int i, const int j) { return j * nx + i; }
  inline int idx(const int i, const int j) { return i * 2 + j; }

public:
  //******************************************************
  // Constructor of the class IVP_ODE_brusselator
  //******************************************************
  IVP_ODE_brusselator (const int nx_points);

 
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0) 
  { 
    for (int i=0;i<nx;i++)
    {  double x_i=(double)(i+1)*dtx;
      Y0[idx(i,0)]=A+sin(2*pi*x_i);
      Y0[idx(i,1)]=B;
    }  
  }

  //******************************************************
  //vector system function for the nonstiff term DY=G(t,Y)
  //******************************************************
  void G (const double t, const double *Y, double *DY);

  //******************************************************
  //vector system function for the nonstiff term DY=F(t,Y)
  //******************************************************
  void F (const double t, const double *Y, double *DY);
  

  //******************************************************
  //Compute sparse matrix As=I-a*Jf, where Jf= Forward Diff 
  //Approx of the Jacobian of the function feval in (t,Y)
  //******************************************************  
  void Compute_Matrix_FD(const double t, double a, double *Y, LIS_MATRIX * As);

  //******************************************************
  //Compute matrix I-a*Jf where Jf= Jacobian of the function 
  // feval defining the Advection-Diffusion-1D IVP_ODE_brusselator
  //******************************************************       
  void Compute_Matrix_Exact(const double t, double a, double *Y, LIS_MATRIX As);


  //******************************************************
  //Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
  // in the IVP_ODE
  //****************************************************** 
  void Compute_G_Matrix_Exact(const double t, double a,double *Y, LIS_MATRIX As);

  //******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
  void Compute_Feval_Jacobian_exact(const double t,double *Y, 
                                    LIS_MATRIX As);

};


// Constructor of the class IVP_ODE_brusselator
IVP_ODE_brusselator::IVP_ODE_brusselator(const int nx_points)
//*****************************************************
{ IVP_name="Brusselator_1D";
  nx=nx_points;
  // Number of ODEs
  neqn=2*nx;
  // Compute Spatial step
  dtx=1.0/(nx+1);
  dtx2=dtx*dtx;
  DD=alpha/dtx2; 
  nnz_G = 3 * neqn - 4;
  nnz_FEVAL = nnz_G + neqn;
 }



double IVP_ODE_brusselator::f(const double y)
{ 
    return(((y-0.7)*(y-1.3))/((y-0.7)*(y-1.3)+0.1));
}

//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//**************************************************
void IVP_ODE_brusselator::G (const double t, const double *Y, double *DY)
{
  const double C[2]={A,B};
  for (int j = 0; j < 2; j++) {
    const int first   = idx(0,j), last=idx(nx-1,j); 
    DY[first] = DD * (Y[first+2] - 2.0 * Y[first] + C[j]);
    DY[last]  = DD * ( C[j]- 2.0 * Y[last] + Y[last-2] );
  }
  for (int i = 1; i < nx-1; i++) {
    for (int j = 0; j < 2; j++) {
      const int ij   = idx(i,j); 
      DY[ij] = DD * (Y[ij+2] - 2.0 * Y[ij] + Y[ij-2]);
    }
  }      
}

//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_brusselator::F(const double t, const double* Y, double* DY)
{
  for (int i = 0; i < nx; i++) {
    const int i0=idx(i,0), i1=idx(i,1);
    const double ui = Y[i0];
    const double vi = Y[i1];
    const double u2v=ui*ui*vi;
    DY[i0] = A+u2v-(B+1)*ui;
    DY[i1] = B*ui-u2v; 
    //DY[i0] = A+ ui*vi-(B+1)*ui;
    //DY[i1] = B*ui+vi; 
  }
}


//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the Brusselator 1D model
//******************************************************       
void IVP_ODE_brusselator::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX As)
{
  
  LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  
  const double Aii = 1 + 2.0 * A * DD;
  const double other = -A * DD;

  init_row_insertion(ptr, k, row);
  // row  (0,j) 
  for (int j = 0; j <= 1; j++) {
    new_entry(index, value, k, row, Aii);
    new_entry(index, value, k, row+2, other);
    next_row(ptr, k, row);
  }
    // row(i,j), i=1,...,nx-2
  for (int i = 1; i <= nx - 2; i++) {
    for (int j = 0; j <= 1; j++) {
      new_entry(index, value, k, row-2, other);
      new_entry(index, value, k, row, Aii);
      new_entry(index, value, k, row+2, other);
      next_row(ptr, k, row);
    }
  }
  // row (nx-1,j)
  for (int j = 0; j <= 1; j++) {
      new_entry(index, value, k, row-2, other);
      new_entry(index, value, k, row, Aii);
      next_row(ptr, k, row);
  }
 
}





//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_brusselator::Compute_Feval_Jacobian_exact(const double t, double* Y, LIS_MATRIX As)
{
  // Number of non-zero entries in As
  LIS_INT k,row;

  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;

  
  int i;
  double Aii = -2.0*DD;
  const double other = DD;
  const double B1 = B + 1.0;
  double ui,vi;

  init_row_insertion(ptr, k, row);
 
  //***  row  (0,j)***
  ui = Y[row];  vi = Y[row+1];
  //i=0, j=0
  new_entry(index, value, k, row, Aii+2* ui * vi - B1);
  new_entry(index, value, k, row+1, ui*ui);
  new_entry(index, value, k, row+2, other);

  next_row(ptr, k, row);

  //i=0, j=1
  new_entry(index, value, k, row-1, B-2*ui*vi);
  new_entry(index, value, k, row, Aii- ui* ui);
  new_entry(index, value, k, row+2, other);

  next_row(ptr, k, row);

  // row(i,j), i=1,...,nx-2
  for (i = 1; i <= (nx - 2); i++) {
    ui = Y[row];  vi = Y[row+1];
    //j=0
    new_entry(index, value, k, row-2, other);
    new_entry(index, value, k, row, Aii+2* ui * vi - B1);
    new_entry(index, value, k, row+1, ui*ui);    
    new_entry(index, value, k, row+2, other);

    next_row(ptr, k, row);

    //j=1
    new_entry(index, value, k, row-2, other);
    new_entry(index, value, k, row-1, B-2*ui*vi);
    new_entry(index, value, k, row, Aii- ui* ui);
    new_entry(index, value, k, row+2, other);

    next_row(ptr, k, row);
  }
  i=nx-1;
  // row (nx-1,j)
  ui = Y[row];  vi = Y[row+1];
  //i=nx-1, j=0
  new_entry(index, value, k, row-2, other);
  new_entry(index, value, k, row, Aii+2* ui * vi - B1);
  new_entry(index, value, k, row+1, ui*ui);   

  next_row(ptr, k, row);

  // i=nx-1, j=1
  new_entry(index, value, k, row-2, other);
  new_entry(index, value, k, row-1, B-2*ui*vi);
  new_entry(index, value, k, row, Aii- ui* ui);

  next_row(ptr, k, row);
    
}

  
#endif 
  
  
