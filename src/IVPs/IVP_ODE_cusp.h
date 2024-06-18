
/*This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_cusp" which represents  
the CUSP model combining the Zeeman’s CUSP catastrophe model and the Van der
Pol oscillator model (see  Hundsdorfer, Verwer 2003: "Numerical Solution of time-dependent advection-diffusion-reaction equations",
Vol. 33 of Springer series in computational mathematics, Springer-Verlag).

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



#ifndef CUSP_H
#define CUSP_H

#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing the CUSP model 
//*****************************************************************
class IVP_ODE_cusp:public IVP_ODE{

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  double DD; //sigma/dtx2
  const double eps=0.0001;
  const double minveps=-(1.0 /eps );
  const double  sigma=1.0/144;   
  const double pi=3.14159265358979;


  // Indexation function which maps 2D spatial coordinates (i,j)
  // to a 1D position in a vector
  //inline int idx(const int i, const int j) { return j * nx + i; }
  inline int idx(const int i, const int j) { return i * 3 + j; }


  //******************************************************************
  // Set values to buid the Jacobian matrix
  //******************************************************************
  void set_values(const int i, double * Y, double *yi_,double *ai_, double *bi_, 
                             double * v_prime); 


  //******************************************************************
  // Set the elements in close to main diagonal in the Jacobian matrix
  //******************************************************************
  void set_tridiag(const int i, const int j,
                         const double  yi, const double ai, const double  bi,
                         const double v_prime, const double Aii, 
                        LIS_INT *index, LIS_INT *k_, LIS_SCALAR * value);  


public:
  //******************************************************
  // Constructor of the class IVP_ODE_cusp
  //******************************************************
  IVP_ODE_cusp (const int nx_points);

  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0); 

  //******************************************************
  //vector system function for the nonstiff term DY=G(t,Y)
  //******************************************************
  void G (const double t, const double *Y, double *DY);

  //******************************************************
  //vector system function for the nonstiff term DY=F(t,Y)
  //******************************************************
  void F (const double t, const double *Y, double *DY);
  

  //******************************************************
  //Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
  // in the IVP_ODE
  //****************************************************** 
  void Compute_G_Matrix_Exact(const double t, double a,double *Y, LIS_MATRIX As);


  //******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
  void Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX As);

};



// Constructor of the class IVP_ODE_cusp
IVP_ODE_cusp::IVP_ODE_cusp(const int nx_points)
//*****************************************************
{ IVP_name="CUSP";
  nx=nx_points;
  // Number of ODEs
  neqn=3*nx;

  // Number of non-zero elements in the Jacobian matrix
  nnz_G=3*neqn;
  nnz_FEVAL=5*neqn;  

  // Compute Spatial step
  dtx=1.0/nx;
  dtx2=dtx*dtx;
  DD=sigma/dtx2; 

 }


  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void IVP_ODE_cusp::init(double *Y0) 
  { 
    for (int i=0;i<nx;i++){ 
      //const double pi_2_xi=2.0*pi*(double)(i+1)*dtx;
        // I made this change to fit the periodic conditions
        const double pi_2_xi = 2.0 * pi * (double)(i) * dtx;
      Y0[idx(i,0)]=0.0;
      Y0[idx(i,1)]=-2*cos(pi_2_xi);
      Y0[idx(i,2)]= 2*sin(pi_2_xi);
    }  
  }

//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//**************************************************
void IVP_ODE_cusp::G (const double t, const double *Y, double *DY)
{
  for (int i = 0; i < nx; i++) {
    const int i0=idx(i,0), i1=i0+1, i2=i0+2;
    const int im1=(i > 0) ? i-1 : nx-1;
    const int ip1=(i < nx-1) ? i+1 : 0;
    DY[i0] = DD * (Y[idx(ip1,0)] - 2.0 * Y[i0] + Y[idx(im1,0)]);
    DY[i1] = DD * (Y[idx(ip1,1)] - 2.0 * Y[i1] + Y[idx(im1,1)]);
    DY[i2] = DD * (Y[idx(ip1,2)] - 2.0 * Y[i2] + Y[idx(im1,2)]);
  }
}

//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_cusp::F (const double t, const double *Y, double *DY)
{   // Compute partially DY in inner points
 for (int i = 0; i < nx; i++) {
    const int i0=idx(i,0), i1=i0+1, i2=i0+2;
    const double yi = Y[i0];
    const double ai = Y[i1];
    const double bi = Y[i2];  
    const double u = (yi - 0.7) * (yi - 1.3);
    const double v = u / (u + 0.1);
    DY[i0] = minveps * (yi*yi*yi + ai*yi + bi);
    DY[i1] = 0.07*v + bi;
    DY[i2] = (1- ai*ai)*bi - ai - 0.4*yi + 0.035*v;
  }

}

//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_cusp
//******************************************************       
void IVP_ODE_cusp::Compute_G_Matrix_Exact(const double t, double A, double *Y,  LIS_MATRIX As)
  {  
    LIS_INT k;

    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
 
    const double Aii=1.0+2*A*DD;
    const double other=-A*DD;
   
    int row=0; ptr[row] = 0; k=0;

    for(int j=0;j<=2;j++){ 
      index[k] = row;        value[k] = Aii;   k++;
      index[k] = row+3;      value[k] = other; k++; 
      index[k] = idx(nx-1,j);value[k] = other; k++; 
      row++;
      ptr[row]=k;
    }

    for(int i=1;i<=nx-2;i++){ 
      for(int j=0;j<=2;j++){ 
        index[k] = row-3;  value[k] = other; k++;
        index[k] = row;    value[k] = Aii;   k++; 
        index[k] = row+3;  value[k] = other; k++; 
        row++;
        ptr[row]=k;
      }   
    }

    for(int j=0;j<=2;j++){ 
      index[k] = idx(0,j); value[k] = other; k++;
      index[k] = row-3;    value[k] = other; k++; 
      index[k] = row;      value[k] = Aii;   k++; 
      row++;
      ptr[row]=k;
    }
    
}



//******************************************************************
// Set values to buid the Jacobian matrix
//******************************************************************
void IVP_ODE_cusp::set_values(const int i, double * Y, double *yi_,double *ai_, double *bi_, 
                             double * v_prime) 
{

  double u, u_prime, denom;
  const int row=idx(i,0);   
   
    *yi_ = Y[row]; *ai_ = Y[row+1]; *bi_ = Y[row+2];
    const double yi=*yi_;
    u = (yi - 0.7) * (yi - 1.3);
    denom=u+0.1;
    u_prime=2*yi-2;
    *v_prime=0.1*u_prime/(denom*denom);
}            




//******************************************************************
// Set the elements in close to main diagonal in the Jacobian matrix
//******************************************************************
void IVP_ODE_cusp::set_tridiag(const int i, const int j,
                         const double  yi, const double ai, const double  bi,
                         const double v_prime, const double Aii, 
                        LIS_INT *index, LIS_INT *k_, LIS_SCALAR * value)  
{
  LIS_INT k=*k_;

  if (j==0){
    index[k] = idx(i,0); value[k] = Aii+minveps* (3*yi*yi + ai);k++;
    index[k] = idx(i,1); value[k] = minveps*yi;k++;    
    index[k] = idx(i,2); value[k] = minveps;k++;          
  }   
  if (j==1){
    index[k] = idx(i,0); value[k] = 0.07*v_prime;k++;  
    index[k] = idx(i,1); value[k] = Aii;k++;
    index[k] = idx(i,2); value[k] = 1.0;k++;
  }  
  if (j==2){
    index[k] = idx(i,0); value[k] = 0.035*v_prime-0.4;k++; 
    index[k] = idx(i,1); value[k] = -2*ai*bi-1;k++;           
    index[k] = idx(i,2); value[k] = Aii+1.0-ai*ai;k++;
  }

  *k_=k;
}            




//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_cusp::Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX As)  
{ 
    LIS_INT k;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;

    const double Aii=-2*DD;
    const double other=DD;
    double yi,ai,bi;
    double v_prime;
  
    int row=0, i=0;
    ptr[row] = 0; k=0;
    set_values(i, Y, &yi,&ai, &bi,&v_prime); 
    for(int j=0;j<=2;j++){ 
      set_tridiag(i,j,yi,ai,bi,v_prime, Aii, index, &k, value);  
      index[k] = row+3;       value[k] = other;        k++; 
      index[k] = idx(nx-1,j); value[k] = other;        k++;       
      row++;
      ptr[row]=k;
    }

    for(i=1;i<=nx-2;i++){ 
      set_values(i, Y, &yi,&ai, &bi,&v_prime); 
      for(int j=0;j<=2;j++){ 
        index[k] = row-3; value[k] = other;         k++;
        set_tridiag(i,j,yi,ai,bi,v_prime, Aii, index, &k, value);  
        index[k] = row+3; value[k] = other;         k++; 
        row++;
        ptr[row]=k;
      }   
    }

    i=nx-1;
    set_values(i, Y, &yi,&ai, &bi,&v_prime); 
    for(int j=0;j<=2;j++){ 
      index[k] = idx(0,j);value[k] = other;        k++;
      index[k] = row-3;   value[k] = other;        k++; 
      set_tridiag(i,j,yi,ai,bi,v_prime, Aii, index, &k, value);  
      row++;
      ptr[row]=k;
    }
    
}

#endif
