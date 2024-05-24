
#ifndef HEAT_H
#define HEAT_H

#include "IVP_ODE.h"
using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing a 1D Advection-Diffusion model 
//*****************************************************************
class IVP_ODE_heat:public IVP_ODE{

private:
  int nx; // number of grid points at x direction
  int ny; // number of grid points at y direction
  double * coeff; // Coefficient vectors for the RHS
  const double l=0.1;
  const double h=0.7;
  const double d=1.38e-7;
  const double Uavg = 1e-6 * 4.806 * (h) * (h) / (12 * 8.9e-4);

  //const double Uavg=4.5e-6;
  const double T0=298.15;
  const double Tin=303.15;

  double dtx; // Spatial step
  double d_dtx2; //d/(dtx*dtx)
  double dty; // Spatial step
  double d_dty2; //d/(dty*dty)
  

  const double pi=3.14159265358979;



  // Indexation function which maps 2D spatial coordinates (i,j)
  // to a 1D position in a vector 
  inline int idx(const int i, const int j)
      { return i*ny+j;}

  // Function which maps a 1D position index to 2D spatial coordinates (i,j) 
  inline void idx2d(const int index, int &i, int &j)
      { i=index/ny;j=index-i*ny;}


public:
  //******************************************************
  // Constructor of the class IVP_ODE-advdiff1d
  //******************************************************
  IVP_ODE_heat(const int nx_points,const int ny_points);


  // Destructor of the class IVP_ODE_heat
  ~IVP_ODE_heat();
  
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0) 
  {
    for (int i=0;i<neqn;i++)
    {      
      Y0[i]=T0;
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


// Constructor of the class IVP_ODE_heat
//IVP_ODE_heat::IVP_ODE_heat(const int nx_points,const int ny_points)
IVP_ODE_heat::IVP_ODE_heat(const int nx_points,const int ny_points)
//*****************************************************
{
  IVP_name="Heat 2D";
  nx=nx_points;
  ny=ny_points;
  // Number of ODEs
  neqn=nx*ny;

  nnz_G = 5 * neqn - 2 * nx - ny ;
  nnz_FEVAL = nnz_G;

  // Compute Spatial step
  dtx=l/nx;
  d_dtx2=d/(dtx*dtx);
  dty=h/ny;
  d_dty2=d/(dty*dty);
  const double h_2=0.5*h;
  coeff =new double[ny];
  for (int j=0;j<ny;j++){
    const double term=(-h_2+(j+0.5)*dty)/h_2;
    coeff[j]=-1.5*Uavg*(1.0-term*term)/dtx;
   }

 }



// Destructor of the class IVP_ODE_heat
IVP_ODE_heat::~IVP_ODE_heat()
//*****************************************************
{ 
  delete[] coeff;
}


//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//***************************************************
void IVP_ODE_heat::G (const double t, const double *Y, double *DY)
{ double ij, im1j,  ijm1, ij1; 
  for(int i=0;i<nx;i++){
	  for (int j=0;j<ny;j++){
      ij   = Y[idx(i, j)];
      im1j = (i == 0)    ? Tin        : Y[idx(i - 1, j)];
      ijm1 = (j == 0)    ? ij : Y[idx(i    , j-1)];
      ij1  = (j == ny-1) ? ij : Y[idx(i    , j+1)];
      const double y_derivative2=d_dty2 * (ijm1 - 2 * ij + ij1);
      const double x_derivative2=( i < nx-1 ) ? d_dtx2 * (im1j             - 2.0 * ij   + Y[idx(i +1, j)]):
                                                d_dtx2 * (-Y[idx(i - 2, j)]+ 2.0 * Y[idx(i - 1, j)] - ij);
      DY[idx(i, j)] = x_derivative2 + y_derivative2;
    }  
  } 
}

//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_heat::F (const double t, const double *Y, double *DY)
{
  double ij, im1j;
  for(int i=0;i<nx;i++){
    for (int j=0;j<ny;j++){
      ij = Y[idx(i, j)];
      im1j = (i == 0) ? Tin: Y[idx(i - 1, j)];
      DY[idx(i, j)] = coeff[j] * (ij -  im1j);
    }
  }
}

//***************************************************



//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_heat
//******************************************************       
void IVP_ODE_heat::Compute_G_Matrix_Exact(const double t, double A,  
                                                double *Y,  LIS_MATRIX As)
{ LIS_INT k,row;
  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  const double Aii_inner=1.0+2*A*(d_dty2+d_dtx2);
  const double Ai0=Aii_inner-A*d_dty2;
  const double Aii_last=1.0+A*(2*d_dty2+d_dtx2);
  const double A0i_last=Aii_last-A*d_dty2;
  const double Aij1=-A*d_dty2;
  const double Ai1j=-A*d_dtx2;
  init_row_insertion(ptr, k, row);
 
  // loop over sparse Jacobian rows
 for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++) {
      // ROW for (i,j)
      if (i<nx-1){ // i=0,...,nx-2 
        if (i>0)    new_entry(index,value,k,idx(i-1,j), Ai1j);        
        if (j>0)    new_entry(index,value,k,idx(i,j-1), Aij1);                          
        if (j<ny-1) {
          new_entry(index,value,k,idx(i,j), (j>0) ? Aii_inner: Ai0); 
          new_entry(index,value,k,idx(i,j+1), Aij1);
        }
      else  new_entry(index,value,k,idx(i,j), Ai0);       
      new_entry(index,value,k,idx(i+1,j), Ai1j);
      }
      else { // i=nx-1  
        new_entry(index,value,k,idx(i-2,j), -Ai1j); 
        new_entry(index,value,k,idx(i-1,j), 2*Ai1j);
        if (j>0) new_entry(index,value,k,idx(i,j-1), Aij1);             
        if (j<ny-1) {
          new_entry(index,value,k,idx(i,j), (j>0) ? Aii_last: A0i_last);
          new_entry(index,value,k,idx(i,j+1), Aij1);
        }
        else new_entry(index,value,k,idx(i,j), A0i_last);
    } 
    next_row(ptr, k, row); 
  }
}
//***************************************************


///******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
void IVP_ODE_heat::Compute_Feval_Jacobian_exact(const double t,double *Y,                                     
                                                     LIS_MATRIX As)
  { LIS_INT k,row;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;

    const double Aii_inner=-2*(d_dty2+d_dtx2);
    const double Ai0=Aii_inner+d_dty2;
    const double Aii_last=-(2*d_dty2+d_dtx2);
    const double A0i_last=Aii_last+d_dty2;
    const double Aij1=d_dty2;
    const double Ai1j=d_dtx2;

    init_row_insertion(ptr, k, row);
 
    // loop over sparse Jacobian rows
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++) {
        // ROW for (i,j)
        if (i<nx-1){ // i=0,...,nx-2 
          if (i>0)    new_entry(index,value,k,idx(i-1,j), Ai1j-coeff[j]);        
          if (j>0)    new_entry(index,value,k,idx(i,j-1), Aij1);                          
          if (j<ny-1) {
            new_entry(index,value,k,idx(i,j),             (j>0) ? Aii_inner+coeff[j]: Ai0+coeff[j]); 
            new_entry(index,value,k,idx(i,j+1),           Aij1);
          }
          else  new_entry(index,value,k,idx(i,j),         Ai0+coeff[j]);       
          new_entry(index,value,k,idx(i+1,j),             Ai1j);
        }
        else { // i=nx-1  
          new_entry(index,value,k,idx(i-2,j),          -Ai1j); 
          new_entry(index,value,k,idx(i-1,j),          2*Ai1j-coeff[j]);
          if (j>0) new_entry(index,value,k,idx(i,j-1), Aij1);             
          if (j<ny-1) {
            new_entry(index,value,k,idx(i,j),         (j>0) ? Aii_last+coeff[j]: A0i_last+coeff[j]);
            new_entry(index,value,k,idx(i,j+1),        Aij1);
          }
          else new_entry(index,value,k,idx(i,j),       A0i_last+coeff[j]);
        } 
        next_row(ptr, k, row); 
    }
  }

#endif