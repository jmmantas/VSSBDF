#ifndef SBDF_SOLVER_H
#define SBDF_SOLVER_H


#include "IVP_ODE.h"

using namespace std;
#include <assert.h>



//*****************************************************************
// Class for the IVP Solver SBDF-order 
//*****************************************************************
class SBDF_Solver {

//*****************************************************
// PRIVATE VARIABLES AND FUNCTIONS
private:
//*****************************************************

  unsigned order; // order of the SBDF IVP solver
  IVP_ODE* IVP; // pointer to the particular IVP Object
  int splitting_type; // splitting type (0=physical splitting, 1=Jacobian-based splitting)
  // Coefficients of SBDF-order method 
  double coef_G[4]={1.0, 2.0/3, 6.0/11, 12.0/25};
  double coef_F[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-2.0/3, 4.0/3, 0.0, 0.0 },
    { 6.0/11,  -18.0/11, 18.0/11, 0.0 },
    {-12.0/25, 48.0/25, -72.0/25, 48.0/25 }
  };
  double coef_Y[4][4]={    
    {1.0, 0.0, 0.0, 0.0 },
    {-1.0/3, 4.0/3, 0.0, 0.0 },
    {2.0/11, -9.0/11, 18.0/11, 0.0 },    
    {-3.0/25, 16.0/25, -36.0/25, 48.0/25 }
  };

  double stability_factor[4]= {5.0 ,2.414,1.501,1.101};
 // double stability_factor[4] = { 1.101 ,1.101,1.101,1.101 };


  double coef_idx;
  int neqn; // number of equations of the IVP
  LIS_MATRIX As,Jf,As2,Jf2; // Auxiliary sparse matrices to store Jacobian matrices
  LIS_SOLVER solver; // LIS Iterative solver
  LIS_VECTOR b, x; // Auxiliary LIS vectors
  double * R; //Auxiliary vector to store the residual of the modified NR iteration
  int idx;
  int nnz; // Number of non-zero entries in Jacobian matrix

  //**************************************************************************
  // Compute the variable term (in a time step) of Right Hand Side at each Newton Iteration
  //**************************************************************************
  void  compute_RHS_SBDF(const double t, const double h, const double* R0, double* Y1,double* R);

  //**************************************************************
  // Check important errors between sparse matrices A1 and A2
  void compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2);
  //**************************************************************

  //**************************************************************************
  // Compute R=R+a*F(t,Y) 
  //**************************************************************************
  void  Add_aFY(const double t, const double a,double* Y);

  //***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  void init_CSR_LIS_matrix(const int splitting_type, LIS_MATRIX *A);
  //***************************************************

  //**************************************************************************
  // Compute the constant term (in a time step) of Right Hand Side 
  // at each Newton Iteration
  //**************************************************************************
  void  compute_RHS0_SBDF(const double t, const double * h_vector, double** Y);

  //***************************************************
  // Init every CSR Sparse LIS matrix
  void init_matrices();
  //***************************************************

  //******************************************************
  // Compute the values for every CSR Sparse LIS matrix
  //******************************************************
  void compute_matrices(const double t, const double h, double ** Y);
  //***************************************************

  //******************************************************
  // Update coefficients of VSSBDF-order scheme
  //******************************************************
  void Update_coefs(const int order,const double * h_vector);
  //******************************************************

  //******************************************************
  // Compute LTE estimate from:
  // the coarse solution Yc_sol and fine solution Yf_sol
  // and the step size coarse vector h_vector
  //******************************************************
  double compute_LTE(const double * h_vector, 
                      double * Yc_sol, double * Yf_sol, double * LTE);
  //******************************************************

  //**************************************************************************
  // Update intermediate vectors Yf and Y using the current vector solution Y1 
  // for the adaptive time stepper before the next integration step
  //**************************************************************************
  void Update_adaptive_intermediate_vectors(double *h_vector_half, double *h_vector_half2, 
                                            double *h_vector, double ** Yf, double ** Y, 
                                            double *Y1);

  //**************************************************************************
  // Update the coarse vector Y1_c from the fine vector Y1_f 
  // by using Richardson extrapolation 
  //**************************************************************************
  void Extrapolate(const double * h_vector, const double *Y1_f, double *Y1_c);
   



//******************************************************
// PUBLIC METHODS
public:
//******************************************************

  //******************************************************
  // Constructor of the class SBDF
  //******************************************************
  SBDF_Solver (const int _order, IVP_ODE* _IVP, const int _splitting_type);

  //*******************************************************************
  // Function implementing one time step using the 
  // SBDF-order Time Integrator (order=1,..,4)
  // splitting  = 0 --> Physical splitting  (Fimp(Y) + Fexp(Y))
  // splitting  = 1 --> Jacobian-based splitting (Jf*Y  + (f(Y)-Jf*Y) )
  //*******************************************************************
  void  Time_Step(const double t, const double *h_vector, double** Y, 
                  double* Y1, const bool variable_tstep);

  //**************************************************************************
  // Update intermediate vectors before the next integration step
  //**************************************************************************
  void Update_intermediate_vectors(double** Y, double* Y1);
  //******************************************************

  //******************************************************
  // Destructor of the class SBDF
  //******************************************************
  ~SBDF_Solver();

  //***************************************************
  // Function implementing the order 1-4 SBDF Time Integrator
  // splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
  // splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
  // It assumes a constant time step h
  //***************************************************
  void Const_dt_SBDF_Integrate(const double t0, const double tf,
                               const double h, double** Y_, 
                               double* Y1);

  //***************************************************
  // Coarse and Fine Time steps using Variable SBDF scheme
  //***************************************************
  void Variable_Time_Step(const double t, double * h_vector, 
                  double * h_vector_half, double * h_vector_half2, double ** Y, 
                  double ** Yf, double * Y1, double * LTE, double *epsilon_c);
  //***************************************************

  //***************************************************
  // Function implementing the order 1-4 SBDF Time Integrator
  // splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
  // splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
  // It assumes a adaptive time step taking into account a stability factor
  //***************************************************
  void Adaptive_dt_SBDF_Integrate0(const double t0, const double tf,
                                  const double h, double** Y_init, 
                                  double ** Yf_init, double* Y1, 
                                  const double tol, int *nsteps, 
                                  int * n_isteps);        
  //***************************************************


  //***************************************************
  // Function implementing the order 1-4 SBDF Time Integrator
  // splitting  =0 -----> Physical splitting is assumed (Fimp(Y) + Fexp(Y))
  // splitting =1 -----> Jacobian-based splitting is assumed (Jf*Y  + (f(Y)-Jf*Y) )
  // It assumes a adaptive time step without stability factor
  //***************************************************
  void Adaptive_dt_SBDF_Integrate1(const double t0, const double tf,
                                   const double h, double** Y_init, 
                                   double ** Yf_init, double* Y1, 
                                   const double tol, int *nsteps, 
                                   int * n_isteps);        
  //***************************************************


  //***************************************************
  // Store the results of an experiment in a data file
  //***************************************************
  void Store_result(const double tf, const double* Y1, 
                    const double tol);
  //***************************************************
   
};


#endif