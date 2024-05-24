#ifndef COMBUSTION_H
#define COMBUSTION_H

#include "IVP_ODE.h"

using namespace std;
#include <algorithm>


//*****************************************************************
// Class for the IVP-ODE  which represents 
// combustion model 
//*****************************************************************
class IVP_ODE_combustion:public IVP_ODE {

private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  double pi_dtx_div_L;
  const double a=1.0;
  const double d=1.0; //constant scalars representing the strength of advection and diffusion 
  const double pi=3.14159265358979;
  const double gamma=0.1;
  const double beta=1.0;
  const double GoB=gamma/beta;
  const double alpha1 = 0.1;
  const double cs=0.6;
  const double alpha2 = alpha1/((1-cs)*(1-cs));
  const double m = 10;
  const double alpha3 = alpha1 / (4 * (pow(m / (m + 1), m) * (1 - m / (m + 1))));
  const double x0=20.5;
  const double xi=10.0;
  const double xf=50.0;
  const double sigma=10;
  const double U0=0.00;
  int Type = 2; // Defining Type of the reaction 1 for FKPP, 2 for Ignition, and 3 for Fisher

  const double L=1;
  //double f_ignition(const double y);
  //double f_Fisher(const double y);
  double reaction(const double y);
  double reaction_derivative(const double y);

public:


    

  //******************************************************
  // Constructor of the class IVP_ODE
  //******************************************************
  IVP_ODE_combustion (const int nx_points);

  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0)
  {
    for (int i=0;i<neqn;i++){ 
     // for (int i = 0; i <= neqn; i++) {

      double x_i=xi+(double)i*dtx;
      double dif=x_i-x0;
     Y0[i]=exp(-dif*dif/sigma);
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
  void Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX As);  


 
};


//*****************************************************
// Constructor of the class IVP_ODE_combustion
IVP_ODE_combustion::IVP_ODE_combustion(const int nx_points)
//*****************************************************
{ IVP_name="Combustion";
  nx=nx_points;
  // Number of ODEs
  neqn=nx;
  // Number of non-zero elements in the Jacobian matrix
  nnz_G=3*nx;
  nnz_FEVAL=3*nx; 
  // Compute Spatial step
  dtx=(xf-xi)/nx;
  dtx2=dtx*dtx;
  pi_dtx_div_L=pi*dtx/L;
  
 }
//***************************************************

double IVP_ODE_combustion::reaction(const double y)
{ 
    if (Type == 1) { 
        return (alpha1 * (y) * (1 - y));
        
    }
    else if (Type == 2) {
        return (alpha2 * (1 - y) * std::max(y - cs, 0.0));
    }
    else if (Type == 3) {
        return(alpha3 * pow(y, m) * (1 - y));
    }
    else {
        std::cerr << "Invalid Reaction Type." << std::endl;
        return 0.0;
        
    }
}



 double IVP_ODE_combustion::reaction_derivative(const double y)
{

    if (Type == 1) {
        return (alpha1 *  (1 - 2*y));

    }
    else if (Type == 2) {
         double r= cs > y ? 0 : alpha2 * (cs - 2 * y + 1);
          //  alpha2 * (1 - y) * std::max(y - cs, 0.0));
         return r;
    }
    else if (Type == 3) {
        return(alpha3 *(m* pow(y, m-1) - (m+1) * pow(y, m)));
    }
    else {
        std::cerr << "Invalid Reaction Type." << std::endl;
        return 0.0;

    }

}

//double IVP_ODE_combustion::f_ignition(const double y)
//{
   
 //   return(alpha2*(1-y)* std::max(y-cs,0.0));
//}

//double IVP_ODE_combustion::f_Fisher(const double y)
//{

 //   return(alpha3 * pow(y,m) * (1-y));
//}



//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//***************************************************
void IVP_ODE_combustion::G (const double t, const double *Y, double *DY)
{
  // Compute partially DY in inner points
  for(int i=1;i<nx-1;i++)  {   
    DY[i]=(1+U0*sin(i*pi_dtx_div_L))*GoB*(Y[i+1] -2*Y[i]+ Y[i-1])/dtx2;
  }
  // Compute partially DY in boundary points (i=0 and i=nx-1)
  DY[0]   =GoB*(Y[1] -2*Y[0]+ Y[nx-1])/dtx2;
  DY[nx-1]=(1+U0*sin((nx-1)*pi_dtx_div_L))*GoB*(Y[0] -2*Y[nx-1]+ Y[nx-2])/dtx2;
}

//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_combustion::F (const double t, const double *Y, double *DY)
{
  
  // Compute partially DY in inner points
    for (int i = 1; i < nx - 1; i++) {
        // DY[i]=-a*(1+U0*sin(i*pi_dtx_div_L))*(Y[i+1]-Y[i-1])/(2*dtx)      +     alpha1*(Y[i])*(1-Y[i]);
        DY[i] = -a * (1 + U0 * sin(i * pi_dtx_div_L)) * (Y[i + 1] - Y[i - 1]) / (2*dtx) + reaction(Y[i]);
    };
  // Compute partially DY in boundary points (i=0 and i=nx-1)
 //DY[0]   =-a*(Y[1]-Y[nx-1])/ 2 * dtx +  alpha1*(Y[0])*(1-Y[0]);
  DY[0] = -a * (Y[1] - Y[nx - 1]) / 2 * dtx + reaction(Y[0]);
 // DY[nx-1]=-a*(1+U0*sin((nx-1)*pi_dtx_div_L))*(Y[0]-Y[nx-2])/ 2 * dtx
                       //         +  alpha1*(Y[nx-1])*(1-Y[nx-1]);
   DY[nx-1]=-a*(1+U0*sin((nx-1)*pi_dtx_div_L))*(Y[0]-Y[nx-2])/ 2 * dtx
                                  + reaction(Y[nx-1]);
}
 
//***************************************************



//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_combustion
//******************************************************       
 void IVP_ODE_combustion::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX As)
  
  { LIS_INT k, row;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;


    const double At2=A/dtx2;
    const double mAt2_GoB=-At2*GoB;
    double other=mAt2_GoB;
    init_row_insertion(ptr, k, row);

    ptr[0] = 0;
    new_entry(index, value, k, 0   , 1.0-2*other);//Jn(0,0)
    new_entry(index, value, k, 1   , other);//Jn(0,1)
    new_entry(index, value, k, nx-1, other);//Jn(0,nx-1)
    next_row(ptr, k, row);

    for(int i=1;i<nx-1;i++) { 
      other=mAt2_GoB*(1+U0*sin(i*pi_dtx_div_L));
      new_entry(index, value, k, i-1   , other);//Jn(i,i-1)
      new_entry(index, value, k, i     , 1.0-2*other);//Jn(i,i)
      new_entry(index, value, k, i+1   , other);//Jn(i,i+1)
      next_row(ptr, k, row);
    }
    other=mAt2_GoB*(1+U0*sin((nx-1)*pi_dtx_div_L));
    new_entry(index, value, k, 0   , other);//Jn(nx-1,0)
    new_entry(index, value, k, nx-2 , other);//Jn(nx-1,nx-2)
    new_entry(index, value, k, nx-1 , 1.0-2*other);//Jn(nx-1,nx-1)
    next_row(ptr, k, row);

} 


//***************************************************

//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_combustion::Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX As)  
  {   
    LIS_INT k;
    
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;



    const double t2=1.0/dtx2;
    const double GoB_t2=GoB*t2;
    const double t1 = 1.0 / (2*dtx);
    const double a_t2=a*t1;
    double other1, other2, F2_i, factor_U0;

    other1=GoB_t2;
    other2=a_t2;
    //F2_i=alpha1*(1.0-2.0*Y[0]); 
    F2_i = reaction_derivative(Y[0]);
    ptr[0] = 0;
    index[0] = 0;    value[0] =-2*other1+F2_i;//Jn(0,0)
    index[1] = 1;    value[1] = other1-other2;//Jn(0,1)
    index[2] = nx-1; value[2] = other1+other2;//Jn(0,nx-1)
    ptr[1]=3; k=3;
    for(int i=1;i<nx-1;i++) {
      factor_U0=(1+U0*sin(i*pi_dtx_div_L));
      other1=GoB_t2*factor_U0;
      other2=a_t2*factor_U0;
      //F2_i=alpha1*(1.0-2.0*Y[i]);  
      F2_i = reaction_derivative(Y[i]);
      index[k] = i-1; value[k] = other1+other2;k++;//Jn(i,i-1) 
      index[k] = i;   value[k] =-2*other1+F2_i;k++;//Jn(i,i)
      index[k] = i+1; value[k] = other1-other2;k++;//Jn(i,i)
      ptr[i+1] = k;
    }
    factor_U0=(1+U0*sin((nx-1)*pi_dtx_div_L));
    other1=GoB_t2*factor_U0;
    other2=a_t2*factor_U0;
    //F2_i=alpha1*(1.0-2.0*Y[nx-1]);  
    F2_i = reaction_derivative(Y[nx - 1]);
    index[k] = 0;    value[k] = other1-other2; k++; //Jn(nx-1,0)
    index[k] = nx-2; value[k] = other1+other2; k++; //Jn(nx-1,nx-2)
    index[k] = nx-1; value[k] = -2*other1+F2_i;k++; //Jn(nx-1,nx-1)
    ptr[neqn] = k;
  
  }


#endif
