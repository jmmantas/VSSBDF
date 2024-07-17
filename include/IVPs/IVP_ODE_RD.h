
/* 
This file contains a C++ subclass of "IVP_ODE" called "IVP_ODE_RD" which represents  
a one-dimensional Reaction-Diffusion model.

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


#ifndef RD_H
#define RD_H
#include "IVP_ODE.h"

using namespace std;



//*****************************************************************
// Class for the IVP-ODE representing the CUSP model 
//*****************************************************************
class IVP_ODE_RD:public IVP_ODE{



double PiGammaInc(double E)
{
    if (E < 0) {
        return sqrt(M_PI) / 2.0;
    } else {
        return 0.5 * sqrt(M_PI) * tgamma(0.5 * E + 1);
    }
}


private:
  int nx; // number of grid points at each dimension
  double dtx; // Spatial step
  double dtx2; //dtx*dtx
  double DD1;
  double DD2;
  //double DD;
  const double  mA=2;
  const double  mC=4;
  const double 	mB=3.5;
  const double 	mY=mA+mC-mB;
  const double 	mZ=2*mY-mB;
  const double  nu1=0.03; //nuA1A2
  const double  nu2=0.04; //nuZ2Z1
  const double  nu3=0.001; //nuB1AC
  const double  nu4=100; //nu11ZB
  const double  nu1A=15000;
  const double  nu1B=15000;
  const double  nu1C=15000;
  const double  nu2A=2;
  const double  nu2B=2;
  const double  nu2C=2;
  const double  nA=1;
  const double  nB=1;
  const double  nC=1;
  const double  EA=4.5;
  const double  EC=4;
  const double  EB=3.6;

//////////////////////
// I changes the E1 and E2 based on the paper parameters
  const double  E1=7.2;
  const double  E2=E1-2;
  const double  EZ=8.2;
////////////////////

  const double  a=pow(((mB*mY)/(mA*mC)),(3/2))*exp(EA+EC-EB-E1)*((nA*nC)/nB);
  const double  b=(nu1*PiGammaInc(E2-E1)*nA)/(nu3*PiGammaInc(EA+EC-EB-E1)*nB);
  const double  d=(nu2*PiGammaInc(E1-E2))/(nu3*PiGammaInc(EA+EC-EB-E1)*nB*nB*pow(((mY*mY)/(mZ*mB)),(3/2))*exp(EZ+EB-E1-E1));
  const double  alfaAY=mA/(mA+mY);
  const double  alfaBY=mB/(mB+mY);
  const double  alfaCY=mC/(mC+mY);
  const double  D1=    1/(mY*(nu1A*nA*alfaAY+nu1B*nB*alfaBY+nu1C*nC*alfaCY)*(nu3*PiGammaInc(EA+EC-EB-E1)*nB));
  const double  D2=  1/(mY*(nu2A*nA*alfaAY+nu2B*nB*alfaBY+nu2C*nC*alfaCY)*(nu3*PiGammaInc(EA+EC-EB-E1)*nB));
  const double  ddelta=(nu2A*nA*alfaAY+nu2B*nB*alfaBY+nu2C*nC*alfaCY)/(nu1A*nA*alfaAY+nu1B*nB*alfaBY+nu1C*nC*alfaCY);


////////////////////// Domain Measure //////////////////////////////
  const double  L=30;
  
 
public:
  //******************************************************
  // Constructor of the class IVP_ODE_RD
  //******************************************************
  IVP_ODE_RD (const int nx_points);


  // Indexation function which maps 2D spatial coordinates (i,j)
  // to a 1D position in a vector
  inline int idx(const int i, const int j) { return j * nx + i; }


  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  inline void init(double *Y0) 
  { srand( (unsigned)time( NULL ) );
         vector<double> rand_vec(2*nx-4);

  for (int i = 0; i < 2*nx-4; i++)
    {
        rand_vec[i] = (static_cast<double>(rand()) / RAND_MAX);
    }

    vector<double> n1_0(nx, 0.0);
    vector<double> n2_0(nx, 0.0);
    const double pert=0.1;

  for (int i = 2; i < nx-2; i++)
    {
        n1_0[i] = a + pert * rand_vec[i-2];
        n2_0[i] = (b / (a * d)) + pert * rand_vec[i-2];

    }

        cout << n1_0[5] << endl;

 	
  for (int i=0;i<nx;i++)
    { //double v=(double) (rand())/RAND_MAX ;
      //double v=0;
      Y0[idx(i,0)]=n1_0[i];
      Y0[idx(i,1)]=n2_0[i];
    }


/*
    for (int i=0;i<neqn;i++)
    { 
      double x_i=(double)(i+1)*dtx;
      if (i<=nx) 
        Y0[i]=a;
      else if(nx<i && i<=2*nx)
        Y0[i]=(b/(a*d));
         }
*/        
    
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
  void Compute_G_Matrix_Exact(const double t, double a,double *Y, LIS_MATRIX * As);


  //******************************************************
  //Compute matrix Jf= Jacobian of the function 
  // feval defining the IVP_ODE
  //******************************************************       
  void Compute_Feval_Jacobian_exact(const double t,double *Y, 
                                    LIS_MATRIX * As);

};


// Constructor of the class IVP_ODE_RD
IVP_ODE_RD::IVP_ODE_RD(const int nx_points)
//*****************************************************
{ IVP_name="Reaction-Diffusion 1D";
  nx=nx_points;
  // Number of ODEs
  neqn=2*nx;
  // Compute Spatial step
  dtx=1.0/nx;
  dtx2=dtx*dtx;
  DD1=D1/dtx2;
  DD2=D2/dtx2; 
 }




//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//**************************************************
void IVP_ODE_RD::G (const double t, const double *Y, double *DY)
{

  for (int j = 0; j <= 1; j++) {
    double DD=( j == 0 ) ? DD1:  DD2;
    for (long i = 0; i < nx; i++) {
      const double ij = Y[idx(i,j)];
      const double im1j = (i == 0) ? Y[idx(i+1, j)] : Y[idx(i - 1, j)];
      const double ip1j = (i == nx-1) ? Y[idx(i-1, j)] : Y[idx(i + 1, j)];
      DY[idx(i,j)] = DD * (ip1j - 2.0 * ij + im1j);
    }
  }
      
}

//***************************************************


//***************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//***************************************************
void IVP_ODE_RD::F(const double t, const double* Y, double* DY)
{
  // Compute partially DY in inner points
//for (long i = 0; i < 2*nx; i++) DY[i]=0.0;

 for (long i = 0; i < nx; i++) {
    const double n1i = Y[idx(i,0)];
    const double n2i = Y[idx(i,1)];
    DY[idx(i,1)] = b*n1i-d*n1i*n1i*n2i ;
    DY[idx(i,0)] = a - n1i - DY[idx(i,1)] ;
  }

}

//***************************************************


//******************************************************
//Compute matrix I-a*Jf where Jf= Jacobian of the function 
// feval defining the Advection-Diffusion-1D IVP_ODE_RD
//******************************************************       
void IVP_ODE_RD::Compute_Matrix_Exact(const double t, double A,double *Y,  
      LIS_MATRIX * As)
  { // Number of non-zero entries in As
    LIS_INT nnz=3*neqn-2;
    //Arrays associated to the CSR format of the matrix As  
    LIS_INT *ptr, *index;
    LIS_SCALAR * value; 
    LIS_INT i,k;

    lis_matrix_malloc_csr(neqn,nnz,&ptr,&index,&value);
    
    const double dtx_2= dtx*2.0;
    const double t2=1.0/dtx2;
    const double Aii=1.0-A*(1.0-2.0*t2);

    ptr[0] = 0;
    index[0] = 0;      value[0] = Aii; //Jn(0,0)
    index[1] = 1;      value[1] = -A*(-Y[1]/dtx_2+t2); //Jn(0,1)
    index[2] = neqn-1; value[2] = -A*(Y[neqn-1]/dtx_2 + t2); //Jn(0,nx-1)
    ptr[1]=3;
    k=3;
    for(i=1;i<neqn-1;i++)
    {  
      index[k] = i-1; value[k] = -A*(Y[i-1]/dtx_2 +t2);//Jn(i,i-1) 
      k++;
      index[k] = i;   value[k] = Aii; //Jn(i,i)
      k++;
      index[k] = i+1; value[k] = -A*(-Y[i+1]/dtx_2 +t2); //Jn(i,i)
      k++;
      ptr[i+1] = k;
    }
    index[k] = 0;      value[k] = -A*(-Y[0]/dtx_2    + t2); //Jn(nx-1,0)
    k++;
    index[k] = neqn-2; value[k] = -A*(Y[nx-2]/dtx_2 + t2); //Jn(nx-1,nx-2)
    k++;
    index[k] = neqn-1;      value[k] = Aii; //Jn(nx-1,nx-1)
    k++;
    ptr[neqn] = k;


    auto start = chrono::high_resolution_clock::now(); 
    lis_matrix_set_csr(nnz,ptr,index,value,*As);

    auto end = chrono::high_resolution_clock::now();
    double runtime_set = chrono::duration_cast<chrono::nanoseconds>(end - start).count()* 1e-9;   
    
    start = chrono::high_resolution_clock::now(); 
    lis_matrix_assemble(*As);
    end = chrono::high_resolution_clock::now();
    double runtime_assemble = chrono::duration_cast<chrono::nanoseconds>(end - start).count()* 1e-9;   
    
    cout<<"****************SET="<<runtime_set<<"****************ASSEMBLE="<<runtime_assemble<<endl;
    
}

//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_RD
//******************************************************       
void IVP_ODE_RD::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX * As)
  { // Number of non-zero entries in As
    LIS_INT nnz=3*neqn-2; 
    LIS_INT *ptr, *index;
    LIS_INT k;
    LIS_SCALAR * value;
    lis_matrix_malloc_csr(neqn,nnz,&ptr,&index,&value);
 
    ptr[0] = 0; k=0;

    for(int j=0;j<=1;j++){
      double DD=( j == 0 ) ? DD1:  DD2;
      const double Aii=1.0+2*A*DD;
      const double otherA1=-2*A*DD ;
      const double other1= -A*DD;

      index[k] = idx(0,j);      value[k] =Aii;  k++; //Jn(j*nx,j*nx)
      index[k] = idx(1,j);      value[k] =otherA1; k++; //Jn(j*nx,j*nx+1)
      ptr[idx(1,j)]=k;

      for(int i=idx(1,j);i<=idx(nx-2,j);i++){  
        index[k] = i-1; value[k] =other1; k++; //Jn(j*nx+i,i-1)
        index[k] = i;   value[k] =Aii;   k++; //Jn(j*nx+i,i)
        index[k] = i+1; value[k] =other1; k++; //Jn(j*nx+i,i+1)
        ptr[i+1] = k;
      }
      
      index[k] = idx(nx-2,j);   value[k] = otherA1; k++;//Jn(j*nx+nx-1,nx-1)
      index[k] = idx(nx-1,j);   value[k] =Aii;      k++; //Jn(j*nx+nx-1,nx-1)
      ptr[idx(nx,j)] = k;
    }



    lis_matrix_set_csr(nnz,ptr,index,value,*As);
      
    lis_matrix_assemble(*As);
}

//***************************************************

//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************       
void IVP_ODE_RD::Compute_Feval_Jacobian_exact(const double t,double *Y, LIS_MATRIX * As)  
{ 

  
}


#endif

