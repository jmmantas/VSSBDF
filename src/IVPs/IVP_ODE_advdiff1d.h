#ifndef ADVDIFF1D_H
#define ADVDIFF1D_H


#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing a 1D Advection-Diffusion model 
//*****************************************************************
class IVP_ODE_advdiff1d:public IVP_ODE {

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  //constant scalars representing the strength of advection and diffusion
  const double a=10,d=1; 
  const double pi=3.14159265358979;
  double f(const double x, const double t);


public:
  //*****************************************************
  // Constructor of the class IVP_ODE_advdiff1d
  IVP_ODE_advdiff1d(const int nx_points);
  //*****************************************************

  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0); 

  //***************************************************
  //vector system function for the stiff term DY=G(t,Y)
  //***************************************************
  void G (const double t, const double *Y, double *DY);

  //***************************************************
  //vector system function for the nonstiff term DY=F(t,Y)
  //***************************************************
  void F (const double t, const double *Y, double *DY);

  //******************************************************
  //Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
  // in the IVP_ODE_advdiff1d
  //******************************************************       
  void Compute_G_Matrix_Exact(const double t, double A,             
                              double *Y,  LIS_MATRIX As);

  //******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
  void Compute_Feval_Jacobian_exact(const double t,double *Y, 
                                    LIS_MATRIX As);
};







//*****************************************************
// Constructor of the class IVP_ODE_advdiff1d
IVP_ODE_advdiff1d::IVP_ODE_advdiff1d(const int nx_points)
//*****************************************************
{ IVP_name="1D_Advection-Diffusion";
  nx=nx_points;
  // Number of ODEs
  neqn=nx;  
  nnz_G=3*neqn;
  nnz_FEVAL=3*neqn;
  // Compute Spatial step
  dtx=1.0/nx;
  dtx2=dtx*dtx;
  
}
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void IVP_ODE_advdiff1d::init(double *Y0) 
  {
    for (int i=0;i<neqn;i++) { 
      double x_i=(double)(i+1)*dtx;
      Y0[i]=sin (2.0*pi*x_i);
    }
  }







double IVP_ODE_advdiff1d::f(const double x, const double t)
{ 
  const double pi2xpt=2.0*pi*x + t;
  const double c=cos(pi2xpt);
  const double s=sin(pi2xpt); 
  return( c + 2*a*pi* s*c + 4*d*pi*pi*s - s);
}

//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//***************************************************
void IVP_ODE_advdiff1d::G (const double t, const double *Y, double *DY)
{
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)
    {   
      DY[i]=d*(Y[i+1] -2*Y[i]+ Y[i-1])/dtx2;
    }

  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =d*(Y[1] -2*Y[0]   + Y[nx-1])/dtx2;
  DY[nx-1]=d*(Y[0] -2*Y[nx-1]+ Y[nx-2])/dtx2;

}

//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_advdiff1d::F (const double t, const double *Y, double *DY)
{
  const double dtx_4=4.0*dtx;
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)
    {   
      DY[i]=-a*(Y[i+1]*Y[i+1]-Y[i-1]*Y[i-1])/dtx_4;
    }

  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =-a*(Y[1]*Y[1]-Y[nx-1]*Y[nx-1])/dtx_4;
  DY[nx-1]=-a*(Y[0]*Y[0]-Y[nx-2]*Y[nx-2])/dtx_4;
  
  // Complete the computation of DY
  for(int i=0;i<nx;i++)
    {   
      DY[i]+=Y[i]+f((i+1)*dtx,t);
    }
}
//***************************************************

//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_advdiff1d
//******************************************************       
void IVP_ODE_advdiff1d::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX As)
  
  { LIS_INT i,k;
    
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    const double t2=1.0/dtx2;

    const double Aii=1.0+2*A*t2*d;
    const double other=-A*t2*d;

    ptr[0] = 0;
    index[0] = 0;      value[0] = Aii; //Jn(0,0)
    index[1] = 1;      value[1] = other; //Jn(0,1)
    index[2] = neqn-1; value[2] = other; //Jn(0,nx-1)

    ptr[1]=3;
    k=3;
    for(i=1;i<neqn-1;i++)
    {  
      index[k] = i-1; value[k] = other;//Jn(i,i-1) 
      k++;
      index[k] = i;   value[k] = Aii; //Jn(i,i)
      k++;
      index[k] = i+1; value[k] = other; //Jn(i,i)
      k++;
      ptr[i+1] = k;
    }
    index[k] = 0;      value[k] = other; //Jn(nx-1,0)
    k++;
    index[k] = neqn-2; value[k] = other; //Jn(nx-1,nx-2)
    k++;
    index[k] = neqn-1;      value[k] = Aii; //Jn(nx-1,nx-1)
    k++;
    ptr[neqn] = k;
  
}
//***************************************************


//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_advdiff1d::Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX As)  
{ 
  LIS_INT i,k;
 
  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  

  const double dtx_2=dtx*2.0;
  const double t2=d/dtx2;
  const double Aii=1.0-2.0*t2;
  ptr[0] = 0;
  index[0] = 0;      value[0] = Aii; //Jn(0,0)
  index[1] = 1;      value[1] = -a*Y[1]/dtx_2+t2; //Jn(0,1)
  index[2] = neqn-1; value[2] = a*Y[neqn-1]/dtx_2 + t2; //Jn(0,nx-1)
  ptr[1]=3;
  k=3;
  for(i=1;i<neqn-1;i++){  
    index[k] = i-1; value[k] = a*Y[i-1]/dtx_2 +t2;//Jn(i,i-1) 
    k++;
    index[k] = i;   value[k] = Aii; //Jn(i,i)
    k++;
    index[k] = i+1; value[k] = -a*Y[i+1]/dtx_2 +t2; //Jn(i,i)
    k++;
    ptr[i+1] = k;
  }
  index[k] = 0;      value[k] = -a*Y[0]/dtx_2    + t2; //Jn(nx-1,0)
  k++;
  index[k] = neqn-2; value[k] = a*Y[nx-2]/dtx_2 + t2; //Jn(nx-1,nx-2)
  k++;
  index[k] = neqn-1;      value[k] = Aii; //Jn(nx-1,nx-1)
  k++;
  ptr[neqn] = k;


}
//******************************************************

#endif
