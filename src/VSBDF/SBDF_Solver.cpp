/*
This is the implementation of the C++ class "SBDF_Solver" which implements the VSSBDF solvers with support for
orders 1 to 4. The solvers apply adaptive time-stepping techniques to optimize the computational
process where the accuracy is ensured by error checking. 

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

#include "IVP_ODE.h"
#include "SBDF_Solver.h"

using namespace std;
#include <assert.h>

//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PRIVATE FUNCTIONS
//*******************************************************************

//**************************************************************
// Check important errors between sparse matrices A1 and A2
void SBDF_Solver::compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2)
//**************************************************************
{
  LIS_INT	n,nnz;
	LIS_INT	i,j,jj1,jj2;
	n   = A1->n;
	nnz = A1->nnz;
    LIS_SCALAR value1, value2;
    int errors=0;
    const double tol=1.0e-2;
    double max_value=A1->value[0];
    max_value=0;
    for(i=1;i<nnz;i++) 
       max_value=max(max_value, fabs(A1->value[i]));

    for(i=0;i<n;i++)
	{
		for(j=A1->ptr[i];j<A1->ptr[i+1];j++)
		{
			jj1 = A1->index[j];
            jj2 = A2->index[j];
			value1 = A1->value[j];
            value2 = A2->value[j];
            double diff=fabs((value2-value1)/max_value);
            if (jj1!=jj2) {
                cout<<"ROW  "<<i<<" ....  COL "<<jj1<<" in A1 is COL "<<jj2<<" in A2"<<endl;
                errors++;
            }  
            else if (diff>tol) {
                cout<<"ROW  "<<i<<" ....  COL "<<jj1<<" Value in A1 is  "<<value1
                                    <<"    but is   "<<value2<<"  in A2"<<endl;
                errors++;
            } 
          
		}
	}
    if (errors>0) {
        cout<<endl<<endl<<"*****************  "<<errors<< 
                       " ERRORS  ***********************" <<endl; 
       // exit(-1);
    } 
}


  //**************************************************************************
  // Compute R=R+a*F(t,Y) 
  //**************************************************************************
  void  SBDF_Solver::Add_aFY(const double t, const double a, double* Y) {   
    double* DY = new double[neqn];
    // Physical splitting
    if (splitting_type==0) {
      IVP->F(t, Y, DY);
      cblas_daxpy(neqn, a, DY, 1, R, 1);
    }
    else { // Jacobian splitting
      double DY2[neqn];
      IVP->feval(t, Y, DY);
      IVP->Jacobian_based_G(t, Y, DY2, Jf);
      cblas_daxpy(neqn, -1.0, DY2, 1, DY, 1);
      cblas_daxpy(neqn, a, DY, 1, R, 1);
    }   
    delete [] DY;
  }

//***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  void SBDF_Solver::init_CSR_LIS_matrix(const int splitting_type, LIS_MATRIX *A) {
    LIS_INT * ptr_A, * index_A;
    LIS_SCALAR* value_A;
    lis_matrix_create(0, A);
    lis_matrix_set_size(*A, 0, neqn);
    lis_matrix_malloc_csr(neqn, nnz, &ptr_A, &index_A, &value_A);
    lis_matrix_set_csr(nnz,ptr_A,index_A,value_A,*A);
    lis_matrix_assemble(*A);
  }
//***************************************************


  //***************************************************
  // Init CSR Sparse LIS matrices As, Jf, Jf2 and As2
  void SBDF_Solver::init_matrices() {
    init_CSR_LIS_matrix(splitting_type, &As);
    init_CSR_LIS_matrix(splitting_type, &Jf); 
    init_CSR_LIS_matrix(splitting_type, &Jf2); 
    init_CSR_LIS_matrix(splitting_type, &As2); 
  }
//***************************************************


//******************************************************
// Compute the values for every CSR Sparse LIS matrix
//******************************************************
  void SBDF_Solver::compute_matrices(const double t, const double h, double ** Y) {
    if (splitting_type==0) {
      //Compute matrix As=Id-a*Jg
      //IVP->Compute_G_Matrix_FD(t, coef_idx, Y[idx], As);

      //IVP->Compute_G_Matrix_FD(t, coef_idx, Y[idx], As2);
      IVP->Compute_G_Matrix_Exact(t, coef_idx, Y[idx], As);
      //compare_matrix_csr(As, As2);
    }
    else {
      //IVP->Compute_Feval_Jacobian_FD(t,Y[idx], Jf);

      //IVP->Compute_Feval_Jacobian_FD(t,Y[idx], Jf2);
      IVP->Compute_Feval_Jacobian_exact(t,Y[idx], Jf);
      //compare_matrix_csr(Jf, Jf2);
      IVP->Compute_Feval_Matrix (t, coef_idx, Jf, As);
    }
  }
//***************************************************


///**************************************************************************
// Compute the variable term (in a time step) of Right Hand Side at each Newton Iteration
//**************************************************************************
void  SBDF_Solver::compute_RHS_SBDF(const double t, const double h, const double* R0, double* Y1, double * R){ 
  double* DY = new double[neqn];
  //R=Y1-R0 
  cblas_dcopy(neqn, Y1, 1, R, 1);
  cblas_daxpy(neqn, -1.0, R0, 1, R, 1);

  if (splitting_type==0) 
     IVP->G(t + h, Y1, DY);
  else 
     IVP->Jacobian_based_G(t+h, Y1, DY, Jf);     

  cblas_daxpy(neqn, -coef_G[idx]*h, DY, 1, R, 1);
  delete[] DY;
}


//**************************************************************************

//**************************************************************************
// Compute the constant term (in a time step) of Right Hand Side (R) at each Newton Iteration
//**************************************************************************
void  SBDF_Solver::compute_RHS0_SBDF(const double t, const double * h_vector, double** Y)
{
    const double h=h_vector[idx];
    double ti[4];
    ti[idx]=t;
    for (int i = idx-1; i >= 0; i--){
      ti[i]=ti[i+1]-h_vector[i];
    }

    //R=Y[0] 
    cblas_dcopy(neqn, Y[0], 1, R, 1);

    if (order==1){
      // SBDF1: Compute R(Y1)=h*F(t,Y[0])+Y[0]
      Add_aFY(t, h, Y[0]);
    }
    else {
      // R=alpha_{-1}*Y[0]
      cblas_dscal(neqn, coef_Y[idx][0],      R, 1);

      for (int i = 1; i <= idx; i++) {
        cblas_daxpy(neqn, coef_Y[idx][i],      Y[i], 1, R, 1);
      }
      for (int i = 0; i <= idx; i++) {
        Add_aFY(ti[i],      coef_F[idx][i]*h,    Y[i]);
      }
            
    }  
} 


//******************************************************
// Update coefficients of VSSBDF-order scheme
//******************************************************
void SBDF_Solver::Update_coefs(const int order, const double * h_vector) {
/*
double coef2_G[4]={1.0, 2.0/3, 6.0/11, 12.0/25};
double coef2_F[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-2.0/3, 4.0/3, 0.0, 0.0 },
    { 6.0/11,  -18.0/11, 18.0/11, 0.0 },
    {-12.0/25, 48.0/25, -72.0/25, 48.0/25 }
  };
double coef2_Y[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-1.0/3, 4.0/3, 0.0, 0.0 },
    {2.0/11, -9.0/11, 18.0/11, 0.0 },    
    {-3.0/25, 16.0/25, -36.0/25, 48.0/25 }
  };
*/


 double w[idx];
 for (int i=0;i<=idx;i++){
   w[i]=h_vector[i+1]/h_vector[i];
 } 

 if (order==2){
    const double denom=1.0+2.0*w[0];
    const double wp1=1.0+w[0];
    coef_G[idx]=wp1/denom;
    coef_Y[idx][0]=-w[0]*w[0]/denom;
    coef_Y[idx][1]=wp1*wp1/denom;
    coef_F[idx][0]=-w[0]*wp1/denom;
    coef_F[idx][1]=coef_Y[idx][1];    
  } 
  else if (order==3){ 
    const double w0=w[0];
    const double w1=w[1];
    const double w0_sq=w0*w0;
    const double w1_sq=w1*w1;
    const double w1p1=1 + w1;
    const double w0p1=1 + w0;
    const double w0w1p1=w0 * w1p1;
      
    const double alpha_0 = -(w0_sq * w1_sq * w0w1p1) / (w0p1 * (w0p1 + w0 * w1));
    const double alpha_1 = w1_sq * (w0 + 1.0 / w1p1);
    const double alpha_2 = -w1p1 - w1 * w0w1p1 / w0p1;
    const double alpha_3 = 1.0 + (w1 / w1p1) + (w0 * w1) / (1.0 + w0w1p1);
    const double beta_0 = (w0_sq * w1 * w1p1) / w0p1;
    const double beta_1 = -w1 * (1.0 + w0w1p1);
    const double beta_2 = (w1p1 * (1.0 + w0w1p1)) / w0p1;
  
    coef_G[idx]=      1.0/alpha_3;
    coef_Y[idx][0]=  - alpha_0/alpha_3;
    coef_Y[idx][1]=  - alpha_1/alpha_3;
    coef_Y[idx][2]=  - alpha_2/alpha_3;
    coef_F[idx][0]=  beta_0/alpha_3;
    coef_F[idx][1]=  beta_1/alpha_3;
    coef_F[idx][2]=  beta_2/alpha_3;
 
  }
  else if (order==4){ 
    const double w0=w[0];
    const double w1=w[1];
    const double w2=w[2];
    const double w0_sq=w0*w0;
    const double w1_sq=w1*w1;
    const double w2_sq=w2*w2;
    const double w0p1=1 + w0;
    const double w1p1=1 + w1;
    const double w2p1=1 + w2;
    
    const double A1 = 1 + w0 * w1p1;
    const double A2 = 1 + w1 * w2p1;
    const double A3 = 1 + w0 * A2;  


    const double alpha_0 = -(w2p1 /w0p1) * (A2 / A1) * (w0_sq*w0_sq * w1_sq*w1 * w2_sq) / A3;
    const double alpha_1 = w1_sq*w1 * w2_sq * (w2p1 /w1p1) * (A3 / A2);
    const double alpha_2 = -w2 * ((w2 / w2p1) + w1 * w2 * (A3 + w0) /w0p1);    
    const double alpha_3 = 1 + w2 * (1 + ((w1 * w2p1) /w1p1) * (1 + (w0 * A2) / A1));
    const double alpha_4 = 1 + (w2 / w2p1) + (w1 * w2 / A2) + (w0 * w1 * w2 / A3);
    const double beta_0 = -pow(w0, 3) * w1_sq * w2 * (w2p1 /w0p1) * (A2 / A1);
    const double beta_1 = w1_sq * w2 * (w2p1 /w1p1) * A3;
    const double beta_2 = -A2 * A3 * (w2 /w0p1);
    const double beta_3 = (w1 * w2p1 /w1p1) * (w2p1 * (A3 + w0) + w0p1 / w1) / A1;

    coef_G[idx]=1.0/alpha_4;
    coef_Y[idx][0]=alpha_0/alpha_4;
    coef_Y[idx][1]=alpha_1/alpha_4;
    coef_Y[idx][2]=alpha_2/alpha_4;
    coef_Y[idx][3]=alpha_3/alpha_4;
    coef_F[idx][0]=beta_0/alpha_4;
    coef_F[idx][1]=beta_1/alpha_4;
    coef_F[idx][2]=beta_2/alpha_4;
    coef_F[idx][3]=beta_3/alpha_4; 
 
  }

 }

  //******************************************************
  // Compute LTE estimate from:
  // the coarse solution Yc_sol and fine solution Yf_sol
  // and the step size coarse vector h_vector
  //******************************************************
  double SBDF_Solver::compute_LTE(const double * h_vector, 
                      double * Yc_sol, 
                      double * Yf_sol, double * LTE)
  //******************************************************
  { 
    double factor;
    const double h_np1=h_vector[idx]; // h_{n+1}
    const double h_n  =h_vector[idx-1]; //h_n
    const double h_n_sq  =h_n*h_n; 
    const double h_np1_sq=h_np1*h_np1;
 
    switch(order) {
      case 1:
        factor=2.0;  break;

      case 2:
      { 
        factor=h_np1_sq*(h_np1 + h_n); //h_{n+1}^2 *(h_{n+1}+h_n)
        const double denom=(5.0*h_np1+7.0*h_n)*h_np1_sq/8.0; 
        factor=factor/denom;    
        break;
      }

      case 3:
      { /*const double h_nm1=h_vector[idx-2];
        const double denom=15.0*h_n_sq +14.0*h_n*h_nm1 
                            + 30.0* h_n*h_np1 + 10.0*h_np1*h_nm1 + 11.0*h_np1_sq;
        factor=16.0*(h_n+h_np1)*(h_np1 + h_n + h_nm1)/denom;
        break;*/

        const double h_nm1=h_vector[idx-2];
        const double denom=14.0*h_n_sq +16.0*h_n*h_nm1 
                            + 27.0* h_n*h_np1 + 16.0*h_np1*h_nm1 + 11.0*h_np1_sq;
        factor=16.0*(h_n+h_np1)*(h_np1 + h_n + h_nm1)/denom;
        break;
      }  

      case 4:  
      {  
        const double h_nm1  =h_vector[idx-2]; //h_{n-1}
        const double h_nm2  =h_vector[idx-3]; //h_{n-2}
        const double h_nm1_sq=h_nm1*h_nm1;

        const double C_1 = h_n_sq*(28.0*h_n+32.0*h_nm2+62.0*h_nm1+84.0*h_np1) 
                           + 32.0*h_n *(h_nm2*h_nm1+2.0*h_nm2*h_np1+h_nm1_sq);
                           
        const double C_2 = 125.0 * h_n * h_nm1 * h_np1 + 79.0 *h_n*h_np1_sq + 32.0*h_nm2*h_nm1*h_np1
            + 32.0 * h_nm2 * h_np1_sq + 32.0 * h_nm1_sq * h_np1 + 63.0 * h_nm1 * h_np1_sq
            + 23.0 * h_np1_sq * h_np1;

        const double denom = C_1 + C_2;

        double term = (h_np1 + h_n);
        term= term*(term + h_nm1) * (term + h_nm1 + h_nm2);
        factor= (32.0 * term)/denom;    
        break;




      }
      default: 
        cout <<"order should belong to {1,2,3,4} "<<endl; exit(-1);
    }  
    // Compute difference LTE=factor*(Yc-Yf)
    cblas_dcopy(neqn, Yc_sol, 1, LTE, 1);
    cblas_daxpy(neqn, -1.0, Yf_sol, 1, LTE, 1);
    cblas_dscal(neqn, factor,      LTE, 1);
    double error = cblas_dnrm2(neqn, LTE, 1);// / pow(neqn, 0.5);
    return(error);
  }
  //******************************************************


  //**************************************************************************
  // Update intermediate vectors Yf and Y using the current vector solution Y1 
  // for the adaptive time stepper before the next integration step
  //**************************************************************************
  void SBDF_Solver::Update_adaptive_intermediate_vectors(double *h_vector_half, 
                                                         double *h_vector_half2, 
                                                         double *h_vector, 
                                                         double ** Yf, double ** Y, 
                                                         double *Y1)
  {  
    //Update h_vectors for the coarse and fine approximation  
    for (int i=0;i<idx;i++) {
      h_vector_half[i]=h_vector_half2[i+1];
      h_vector[i]=h_vector[i+1];
    }

    // Update coarse grain approximation vectors Y
    Update_intermediate_vectors(Y, Y1);

    // Update fine grain approximation vectors Yf from coarse vectors Y
    switch(order) {
      case 1:
        cblas_dcopy (neqn, Y[idx],1,    Yf[0],1);   break;
      case 2:
      { cblas_dcopy (neqn, Yf[order],1, Yf[0],1);
        cblas_dcopy (neqn, Y[idx],1,    Yf[1],1);   break;
      }
      case 3:
      { cblas_dcopy (neqn, Y[1],1,      Yf[0],1);
        cblas_dcopy (neqn, Yf[order],1, Yf[1],1);
        cblas_dcopy (neqn, Y[idx],1,    Yf[2],1);   break;
      }
      case 4:
      { cblas_dcopy (neqn, Yf[2],1, Yf[0],1);
        cblas_dcopy (neqn, Y[2],1,        Yf[1],1);
        cblas_dcopy (neqn, Yf[4],1,   Yf[2],1);
        cblas_dcopy (neqn, Y[idx],1,      Yf[3],1); break;
      }
    }
  }    
  
  
   
 //**************************************************************************
// Update intermediate vectors before the next integration step
//**************************************************************************
/*void Update_intermediate_vectors(const int order, const int neqn, double** Y, double* Y1)
{   double * Ytmp=Y[0];
    for (int i = 1; i < order; i++) {
        Y[i-1]=Y[i];
    }    
    Y[order-1]=Ytmp;
    cblas_dcopy(neqn, Y1, 1, Y[order-1], 1);
}*/

  //**************************************************************************
  // Update the coarse vector Y1_c from the fine vector Y1_f 
  // by using Richardson extrapolation 
  //**************************************************************************
  void SBDF_Solver::Extrapolate(const double * h_vector, 
                                const double *Y1_f, double * Y1_c){
    double alpha, beta;
    if (order==1) {
      alpha=-1.0;
      beta=2.0;
    } 
    else if (order==2) {
      const double denom=7.0*h_vector[0]+5.0*h_vector[1];
      alpha=-(h_vector[0]+3.0*h_vector[1])/denom;
      beta=8.0*(h_vector[0]+h_vector[1])/denom;
    }  
    else if (order==3) {
      const double h1_sq=h_vector[1]*h_vector[1];
      const double h2_sq=h_vector[2]*h_vector[2];
      const double h1h0=h_vector[1]*h_vector[0];
      const double h2h1=h_vector[2]*h_vector[1];
      const double h2h0=h_vector[2]*h_vector[0];
      const double h2ph1=h_vector[2]+h_vector[1];
      const double denom=15.0*h1_sq +14.0*h1h0 + 30.0*h2h1+10.0*h2h0 + 11.0*h2_sq;
      alpha=-(h1_sq +2.0*(h1h0+h2h1+3.0*h2h0)+5.0*h2_sq)/denom;
      beta=16.0*(h2ph1*(h2ph1+h_vector[0]))  /denom;
    }
    else //(order==4) 
    {

      const double h_np1=h_vector[idx]; // h_{n+1}
      const double h_n  =h_vector[idx-1]; //h_n
      const double h_nm1  =h_vector[idx-2]; //h_{n-1}
      const double h_nm2  =h_vector[idx-3]; //h_{n-2}

      const double C_1 = 31.0 * pow(h_n, 3) + 30.0 * pow(h_n, 2) * h_nm2 + 60.0 * pow(h_n, 2) * h_nm1
            + 93.0 * pow(h_n, 2) * h_np1 + 28.0 * h_n * h_nm2 * h_nm1 + 60.0 * h_n * h_nm2 * h_np1
            + 28.0 * h_n * pow(h_nm1, 2);
      const double C_2 = 120.0 * h_n * h_nm1 * h_np1 + 93.0 * h_n * pow(h_np1, 2) + 20.0 * h_nm2 * h_nm1 * h_np1
            + 22.0 * h_nm2 * pow(h_np1, 2) + 20.0 * pow(h_nm1, 2) * h_np1 + 44.0 * h_nm1 * pow(h_np1, 2)
            + 23.0 * pow(h_np1, 3);
      const double term = (h_np1 + h_n) * (h_np1 + h_n + h_nm1) * (h_np1 + h_n + h_nm1 + h_nm2);
      const double denom = C_1 + C_2;
      alpha = ((denom-32)*term) / denom;
      beta = (32 *term) / denom;
    }  



    cblas_dscal (neqn, alpha, Y1_c , 1);
    cblas_daxpy (neqn, beta, Y1_f , 1, Y1_c, 1);    

  }



//*******************************************************************
//*******************************************************************
// IMPLEMENTATION OF THE PUBLIC FUNCTIONS
//*******************************************************************

//******************************************************
// Constructor of the class SBDF
//******************************************************
SBDF_Solver::SBDF_Solver (const int _order, IVP_ODE* _IVP, 
                                const int _splitting_type)
//*****************************************************
{ 
  order=_order;
  idx=order-1; 
  IVP=_IVP;
  splitting_type=_splitting_type;
  neqn = IVP->get_num_ODEs();
  nnz=IVP->get_nnz(splitting_type);
  
  // Create LIS vectors
  lis_vector_create(0, &b);
  lis_vector_create(0, &x);
  lis_vector_set_size(b, 0, neqn);
  lis_vector_set_size(x, 0, neqn);

  R = new double[neqn];
  
  // Create and initialize LIS solver with the suitable switches
  lis_solver_create(&solver);
  //char opts[] = "-i bicgstab -p ilut  -tol 1.0e-12";
  //char opts[] = "-i icr  -p jacobi -tol 1.0e-12";
  char opts[] = "-i fgmres -p ilut  -tol 1.0e-15";
  lis_solver_set_option(opts, solver);
 }

 //******************************************************
// Destructor of the class SBDF
//******************************************************
SBDF_Solver::~SBDF_Solver (){

  // Destroy LIS vectors b ad x
  lis_vector_destroy(x);
  lis_vector_destroy(b);

  // Destroy LIS matrices As and Jf
  lis_matrix_destroy(As);
  lis_matrix_destroy(Jf);
  lis_matrix_destroy(Jf2);
  lis_matrix_destroy(As2); 

  // Destroy the LIS solver
  lis_solver_destroy(solver);
  delete[] R;
  
}

//**************************************************************************
// Update intermediate vectors before the next integration step
//**************************************************************************
void SBDF_Solver::Update_intermediate_vectors(double** Y, double* Y1){   
  double * Ytmp=Y[0];
  for (unsigned i = 1; i < order; i++) {
      Y[i-1]=Y[i];
  }    
  Y[idx]=Ytmp;
  cblas_dcopy(neqn, Y1, 1, Y[idx], 1);
}




//***************************************************
// Store the results of an experiment in a data file
//***************************************************
void SBDF_Solver::Store_result(const double tf, 
                              const double* Y1, 
                              const double tol)
//***************************************************

{
  ostringstream sstr;
  string splitting=(splitting_type==0)?"PS":"JS";
  sstr << IVP->get_name()<<"_" << neqn<<"_"<<splitting 
          <<"_tf-"<<tf<<"_tol-"<<tol << ".dat";
  ofstream str;


  str.open( sstr.str().c_str(),ios_base::out );

  if (!str.is_open()) {
        std::cerr << "Failed to open: " << sstr.str().c_str()<< std::endl;
        return;
  }      
  str << scientific << setprecision(16);
  str << IVP->get_name() << endl;
  str << neqn << endl;
  for( int i=0; i<neqn; ++i )
    {
	    str << Y1[i] << endl;
     }

     str.close();

     return;
}


//*******************************************************************
// Function implementing one time step using the 
// SBDF-order Time Integrator (order=1,..,4)
// splitting  = 0 --> Physical splitting  (Fimp(Y) + Fexp(Y))
// splitting  = 1 --> Jacobian-based splitting (Jf*Y  + (f(Y)-Jf*Y) )
//*******************************************************************
void  SBDF_Solver::Time_Step(const double t, const double *h_vector, double** Y, double* Y1, 
const bool variable_tstep){

  if (variable_tstep) {
   Update_coefs(order, h_vector); 
  }
  const double h=h_vector[idx];
  coef_idx=coef_G[idx]*h;
  compute_matrices(t,h,Y);
  // Get initial approximation of Y1=Yidx+h*Yidx
  cblas_dcopy(neqn, Y[idx], 1, Y1, 1); 
  cblas_daxpy(neqn, h, Y[idx], 1, Y1, 1);
  // Modified Newton Iteration to approximate Y1 at the current time step
  double norm = 1.0e23;
  int it = 0;
  bool convergence = false;
  const double EPSTOL = 1.0e-12;
  const int max_iterations = 200;
  // Compute the constant term in a time step of vector R
  compute_RHS0_SBDF(t, h_vector, Y);
  while (!convergence) {
    it++;
    // Complete the  RHS vector R, adding the variable terms  
    compute_RHS_SBDF(t, h, R, Y1, b->value);
    // Initialize vector x of the linear solver and b=R
    // to solve iteratively the system A*x=b to approximate the residual vector x
    lis_vector_set_all(0, x);
    int error = lis_solve(As, b, x, solver);
    if (error != 0) {
      cout << "######## ############# ERROR IN LINEAR SYSTEM SOLUTION" << endl; exit(0);
    }
    //Update Y1=Y1-x
    cblas_daxpy(neqn, -1.0, x->value, 1, Y1, 1);
    // Compute the norm of the residual vector x
    norm = cblas_dnrm2(neqn, x->value, 1) / sqrt((double)neqn);
    convergence = (norm < EPSTOL) || (it > max_iterations);
  }  // End of modified Newton Iteration
  //cout<< "TIME="<<t<< " ............#NEWTON_ITERATIONS="<<it<<endl;
  if (it >= max_iterations) {
    cerr << ".................. ##################### Max iterations achieved!!!!" << endl;
    return;
  }

}


//***************************************************
// Function implementing the order 1-4 SBDF Time Integrator
// splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
// splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
// It assumes a constant time step h
//***************************************************
void SBDF_Solver::Const_dt_SBDF_Integrate(const double t0, const double tf,
    const double h, double** Y_, double* Y1)  {        
    //***************************************************
    const double EPSTOL=1.5e-14;
    const int neqn = IVP->get_num_ODEs();
    // Initialize intermediate stage vectors Y
    double *Y[order];
    for (unsigned i = 0; i < order; i++) {
        Y[i] = new double[neqn];
        cblas_dcopy(neqn, Y_[i], 1, Y[i], 1); 
    }
       
    double t = t0;  
    const int idx=order-1;
    //int Num_steps = (int)((tf - t0) / h)-idx;
    t = t + idx * h;
    double h_vector[4]={h,h,h,h};
    init_matrices();
    bool variable_tstep=false;
    int N_steps=0;
    double h_now=h;
    bool end=false;
    //*********  TIME LOOP  ****************
    while (!end) {
    //**************************************
      const double distance = tf - t;  
      if (distance<h) {    
        h_now=distance;
        h_vector[idx]=h_now;
        variable_tstep=true;                                 
      } 
 
      Time_Step(t, h_vector, Y, Y1, variable_tstep);

      Update_intermediate_vectors(Y, Y1);
      t+=h_now;  
      N_steps++;
      if (fabs(tf-t)<EPSTOL) {end=true;}
    } // End of time stepping

    for (unsigned i = 0; i < order; i++) 
        delete[] Y[i];
}


//********************************************************************************************
// Compute Coarse and Fine Time steps using Variable SBDF scheme as well as the LTE
//********************************************************************************************
void SBDF_Solver::Variable_Time_Step(const double t, double * h_vector, 
                                     double * h_vector_half, double * h_vector_half2, double ** Y, 
                                     double ** Yf, double * Y1, double * LTE, double *epsilon_c)
{ const bool  variable_tstep=true;

  // COARSE STEP
  Time_Step(t, h_vector, Y, Y1, variable_tstep);

  h_vector_half[idx]=h_vector[idx]/2.0;  
  for (int i=0;i<idx;i++) {
    h_vector_half2[i]=h_vector_half[i+1];
  }
  h_vector_half2[idx]=h_vector[idx]/2.0;

  // FINE STEP 1
  Time_Step(t,h_vector_half, Yf, Yf[order], variable_tstep);

  //FINE STEP 2
  Time_Step(t+h_vector_half[idx], h_vector_half2, &(Yf[1]), Yf[order+1], variable_tstep);
  
  // Computation of the LTE vector and the scalar error
  *epsilon_c=compute_LTE(h_vector, Y1, Yf[order+1],LTE);
}



//***************************************************
// Function implementing the order 1-4 SBDF Time Integrator
// splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
// splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
// It assumes a adaptive time step
//***************************************************
void SBDF_Solver::Adaptive_dt_SBDF_Integrate0(const double t0, const double tf,
    const double h, double** Y_init, double ** Yf_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps)  {        
//***************************************************
    // Several constants used in the adaptive time stepping
    const double  p_inv=1.0/(order+1), range = tol/3, 
                  alpha=0.8, eta_min = 0.2, eta_max = 5,
                  EPSTOL=1.0e-20;
    const unsigned i_max = 5;

    // index:  indicates the position of the last solution vector (in the last integration step)
    // in the intermediate step vectors Y and Yf 
    const int idx=order-1;

    //neqn: number of ODEs
    const int neqn = IVP->get_num_ODEs();

    // Initialize intermediate step vectors Y, Yf and the LTE vector
    double *Y[order], *Yf[order+2];
    double *LTE=new double[neqn];
    for (unsigned int i = 0; i < (order+2); i++) {
        Yf[i] = new double[neqn]; 
    }
    for (unsigned int i = 0; i < order; i++) {
        Y[i] = new double[neqn];
        cblas_dcopy(neqn, Y_init[i], 1, Y[i], 1);
        cblas_dcopy(neqn, Yf_init[i], 1, Yf[i], 1);    
    }
    // Initialize next time instant t in terms of the initial time t0, 
    // the initial stepsize  h  and the order (idx=order-1) of the method
    double t = t0 + idx * h;

    // Vectors of intermediate time steps h_vector, h_vector_half and h_vector_half2  
    double h_vector[4]={h,h,h,h}; // for coarse grain approximation
    double h_vector_half[4]={h/2,h/2,h/2,h/2}; // for the 1st step in the fine grain approximation
    double h_vector_half2[4];// for the 2nd step in the fine grain approximation
    
    // Init CSR Sparse LIS matrices: As, Jf, Jf2 and As2 before compute them
    init_matrices();

    // Counter for the number of time steps
    int N_steps=0;
    // Boolean to indicate the end of the time loop  
    bool end=false;
    // Total number of inner iterations
    unsigned total_isteps=0;

    //*********  TIME LOOP  ****************
    while (!end) {
    //**************************************
      // Distance to the final time
      const double distance = tf - t;
      // accepted: indicates if the current time step is accepted
      bool accepted=false;
      // last_step: indicates if the current time step is probably the last one
      bool last_step=false;
      // imposed_stability: indicates if the stability factor has been applied in the current time step
      bool imposed_stability=false;
      // Local trucation error epsilon_c
      double epsilon_c;
      // number of iteration for a single adaptive time step
      unsigned i_step=0;
      // factor: value used to update the time step by multiplying the current time step
      double factor;   

      // INNER LOOP FOR A TIME STEP
      while (! accepted)  {
    
        // If the current step size is greater than tf-t, 
        // then new step size= tf-t and it is probably the last step 
        if (h_vector[idx]>distance) {
          h_vector[idx]=distance; 
          last_step=true;
          //cout<<"** T= "<<t<<"   ***********COMPUTING FINAL DT="
          //            <<h_vector[idx]<<"  *****  h_const= "<<h<<endl<<flush;        
        } 

        // Compute a coarse and fine time steps in Y1, using the VSSBDF scheme
        // and the error estimate epsilon_c
        Variable_Time_Step(t, h_vector, h_vector_half, h_vector_half2, Y, Yf, Y1, LTE, &epsilon_c);
       
        // Check if the current step size is accepted
        accepted=(fabs(epsilon_c-tol)<=range);


        // If the stability factor has been applied 
        // or it is the last time step, 
        // or the number of inner iterations is greater than i_max
        // then the current step size is accepted when the error is sufficiently small  
        if ( i_step>i_max || imposed_stability || last_step){
          accepted=(epsilon_c <= (tol+range)); // Weaker acceptance condition in special cases
          imposed_stability=false;
        }

        // when the step size is not accepted and should be updated
        if (!accepted) {      
          // If the number of inner iterations gets its limit 
          // and the error is greater than the tolerance
          // then update the factor without any restriction
          if (i_step>i_max && epsilon_c > tol) {
            // Use more agressive update when we are stuck
            factor=alpha*pow(tol/epsilon_c,p_inv); 
          }
          else {
            // Update the factor using the known formula
            factor=min(max(alpha*pow(tol/epsilon_c,p_inv),eta_min  ),eta_max);  
          }  
          // Update the step size using the current factor
          h_vector[idx]*=factor;  

          // If hnew/hold is greater than the stability limit, update the step size 
          if (order>1  && ( (h_vector[idx]/h_vector[idx-1])>= stability_factor[idx])){
            imposed_stability=true;
            h_vector[idx]=h_vector[idx-1]*stability_factor[idx];       
          }
        }
        // Update the number of inner iterations for this time step         
        i_step++;     
      } // END OF INNER LOOP OF A TIME STEP


      // Update the total number of inner iterations
      total_isteps+=i_step;
      // Update t using the old step size stored in h_vector[idx]
      t+=h_vector[idx];
      // Update the number of time steps
      N_steps++;
      // Check if the end of the time loop is reached
      end=(fabs(tf-t)<EPSTOL); 

      // Update Y1 using Richardson Extrapolation using LTE: Y1=Y1-LTE  
      cblas_daxpy(neqn, -1.0, LTE, 1, Y1, 1);

      if (!end) { // If the end of the time loop is not reached
        // Update intermediate vectors and stepsizes for the next time step
        Update_adaptive_intermediate_vectors(h_vector_half, h_vector_half2, h_vector, Yf, Y, Y1);
      }
    } // End of time stepping


    //cout << endl << "FINAL TIME= " << t << "    tf-t= "<<tf-t<<
    //"   Actual Number of time steps=" << N_steps <<
    //"   TOTAL COMPUTED TIME STEPS= "<<3*total_isteps<< endl;

    // Set the final values for the output parameters 
    // denoting the number of time steps and inner iterations 
    *nsteps=N_steps; 
    *n_isteps=total_isteps; 

    // Free memory allocated for intermediate vectors Y, Yf and LTE
    for (unsigned i = 0; i < order; i++)
        delete[] Y[i];
    for (unsigned i = 0; i < (order+2); i++)
        delete[] Yf[i];
    delete[] LTE;    
}


//***************************************************
// Function implementing the order 1-4 SBDF Time Integrator
// splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
// splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
// It assumes a adaptive time step
//***************************************************
void SBDF_Solver::Adaptive_dt_SBDF_Integrate1(const double t0, const double tf,
    const double h, double** Y_init, double ** Yf_init, double* Y1, 
    const double tol, int *nsteps, int * n_isteps)  {        
//***************************************************
  const double  p_inv=1.0/(order+1);
  const double range=tol/3.0;
  const double alpha=0.8;
  const double eta_min= 0.2;
 // const double eta_max= stability_factor[idx];
  const double eta_max = 5;
  const double EPSTOL=1.0e-20;
  const int neqn = IVP->get_num_ODEs();
  // Initialize intermediate stage vectors Y, Tf and the LTE vector
  double *Y[order], *Yf[order+2];
  double *LTE=new double[neqn];
  for (unsigned int i = 0; i < (order+2); i++) {
    Yf[i] = new double[neqn]; 
  }
  for (unsigned int i = 0; i < order; i++) {
    Y[i] = new double[neqn];
    cblas_dcopy(neqn, Y_init[i], 1, Y[i], 1);
    cblas_dcopy(neqn, Yf_init[i], 1, Yf[i], 1);    
  }

  // Init the time counter
  double t = t0;  
  const int idx=order-1;
  t = t + idx * h;
  double h_vector[4]={h,h,h,h}; double h_vector_half[4]={h/2,h/2,h/2,h/2};
  double h_vector_half2[4]; double next_dt=h;

  init_matrices();

  int N_steps=0, total_isteps=0;
  bool end=false;
  //*********  TIME LOOP  ****************
  while (!end) {
  //**************************************
    const double distance = tf - t;  
    bool found_h_now=false;
    // Local trucation error epsilon_c and factor to modify dt_new
    double epsilon_c,factor;
    // number of iteration for a single adaptive time step
    unsigned i_step=0;

    h_vector[idx]=next_dt;
    //**************************
    // INNER LOOP OF A TIME STEP
    //**************************
    while (! found_h_now)  {

      // If the hnew is greater than distance=tf-t, then new h=tf-t and last step is possible 
      if (h_vector[idx]>distance) {
        h_vector[idx]=distance; 
        cout<<"** T= "<<t<<"   ***********COMPUTING FINAL DT="
            <<h_vector[idx]<<"  *****  h_const= "<<h<<endl<<flush;        
      }     

      // New coarse and fine integration step
      Variable_Time_Step(t, h_vector, h_vector_half, h_vector_half2, Y, Yf, Y1, LTE, &epsilon_c);

      // Check the error condition
      double error_diff=epsilon_c-tol;
      found_h_now=(error_diff<=range);
      factor=min(max(alpha*pow(tol/epsilon_c,p_inv),eta_min  ),eta_max);

      if (found_h_now){
        if (fabs(error_diff)<=range) next_dt=h_vector[idx];
        else                         next_dt=h_vector[idx]*factor;
      } 
      else { 
        cout<< "TOL="<<tol<<" ..... STEP NOT ACCEPTED    NEXT DT="<< next_dt
            <<"     ERROR= "<<epsilon_c<<"   FACTOR=  "<<factor<<endl<<flush;
        h_vector[idx]*=factor;
      }  
      i_step++;     
    //*************************************  
    }   // END OF INNER LOOP OF A TIME STEP
    //*************************************

    total_isteps+=i_step;
    // Update t using the old time step stored in h_vector[idx]
    t+=h_vector[idx]; N_steps++;

    // Check the end ot time integration 
    end=(fabs(tf-t)<EPSTOL); 

    // Richardson Extrapolation using LTE:   Y1=Y1-LTE
    cblas_daxpy(neqn, -1.0, LTE, 1, Y1, 1);

    if (!end) {
      // Update intermediate vectors and stepsizes
      Update_adaptive_intermediate_vectors(h_vector_half, h_vector_half2, h_vector, Yf, Y, Y1);
    }      
  //*********************
  } // End of time stepping
  //*********************

  cout << endl << "FINAL TIME= " << t << "    tf-t= "<<tf-t<<
    "   Actual Number of time steps=" << N_steps <<
    "   TOTAL COMPUTED TIME STEPS= "<<3*total_isteps<< endl;

  *nsteps=N_steps;  *n_isteps=total_isteps; 
  for (unsigned i = 0; i < order; i++)
    delete[] Y[i];
  for (unsigned i = 0; i < (order+2); i++)
    delete[] Yf[i];
  delete[] LTE;    
}
