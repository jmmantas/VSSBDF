#ifndef Included_IVP_ODE_H

#define Included_IVP_ODE_H


#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <chrono>
#include "lis.h"
#include "lis_config.h"
#include <vector>
#include <cblas.h>
#include <sstream>
#include <iomanip>

using namespace std;

//*****************************************************************
// Class for the IVP-ODE representing a model 
//*****************************************************************
class IVP_ODE {

protected:
  int neqn; // number of equations
  string IVP_name; // Name of the IVP model
  int nnz_G; //number of non-zero elements in the Jacobian for the stiff term G
  int nnz_FEVAL; //number of non-zero elements in the Jacobian for F+G
  

  //***************************************************
  // Introduce a new element "entry" in a CSR sparse matrix, 
  // in the k-th position of the index vector on the column "col"
  //*************************************************** 
  inline void new_entry(LIS_INT* index, LIS_SCALAR* value,  
                      LIS_INT & k, const int col, const double entry){
    index[k] = col;
    value[k] = entry;
    k++;
  }  
  //***************************************************

  //***************************************************
  // Insert a new row in a CSR sparse matrix with pointer vector ptr
  // this new row start in position k of the index vector
  //*************************************************** 
  inline void next_row(LIS_INT* ptr, const LIS_INT k, int & row){
    row++;
    ptr[row] = k;
  }  
  //***************************************************

  //***************************************************
  // Init row 0 in CSR sparse matrix with pointer vector ptr
  // this new row start in position k=0 of the index vector
  //*************************************************** 
  inline void init_row_insertion(LIS_INT* ptr, LIS_INT & k, int & row){
    row = 0;
    ptr[row] = 0;
    k = 0;
  }  
  //***************************************************



//***************************************************
// PUBLIC FUNCTIONS
public:
//***************************************************
  
  //******************************************************
  // Get number of ODEs
  //******************************************************
  inline int get_num_ODEs() {return neqn;};
  //******************************************************


  //***************************************************
  // Get name of the IVP
  inline string get_name() {return IVP_name;}
  //***************************************************                                              

  //***************************************************
  // Get number of nonzero elements in Jacobian
  // depending on the splitting type
  inline int get_nnz(const int splitting_type ) 
  {
    int nnz=(splitting_type==0)? nnz_G: nnz_FEVAL;     
    return nnz;
  }
  //***************************************************                                              


  //***************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  //******************************************************
  inline void init_CSR_LIS_matrix(const int splitting_type, LIS_MATRIX *A) 
  {
    const int nnz=(splitting_type==0)? nnz_G: nnz_FEVAL;
    LIS_INT * ptr_A, * index_A;
    LIS_SCALAR* value_A;
    lis_matrix_create(0, A);
    lis_matrix_set_size(*A, 0, neqn);
    lis_matrix_malloc_csr(neqn, nnz, &ptr_A, &index_A, &value_A);
    lis_matrix_set_csr(nnz,ptr_A,index_A,value_A,*A);
    lis_matrix_assemble(*A);
  }
  //***************************************************                                              

 
  //******************************************************
  // Initialize stage vector Y0 with neqn components
  //******************************************************
  virtual void init(double *Y0)=0; 
 
  //******************************************************
  //vector system function for the stiff term DY=G(t,Y)
  //******************************************************
  virtual void G (const double t, const double *Y, 
                                               double *DY)=0;

  //******************************************************
  //vector system function for the nonstiff term DY=F(t,Y)
  //******************************************************
  virtual void F (const double t, const double *Y, double *DY)=0;


  //*************************************************************
  //ODE vector system function feval: DY=feval(t,Y)=F(t,Y)+G(t,Y)
  //*************************************************************
  //virtual void feval (const double t, const double *Y, double *DY)=0;
  

  //*************************************************************
  //ODE vector system function feval: DY=feval(t,Y)=F(t,Y)+G(t,Y)
  //*************************************************************
  inline void feval (const double t, const double *Y, double *DY)
  { double Ytmp[neqn];
    G (t, Y, DY);
    F (t, Y, Ytmp);
    cblas_daxpy(neqn,1.0,Ytmp,1,DY,1);
  }

  //******************************************************
  //Compute matrix I-a*Jg where Jg= Jacobian of the stiff Term 
  // in the IVP_ODE
  //****************************************************** 
  virtual void Compute_G_Matrix_Exact(const double t, double a,
                                    double *Y, LIS_MATRIX As)=0;
  //******************************************************
  // Compute matrix Jf= Jacobian of the function feval defining the  IVP_ODE
  //******************************************************     
  virtual void Compute_Feval_Jacobian_exact(const double t,double *Y,LIS_MATRIX As)=0;


  ///******************************************************
  // Compute matrix  A=I-a*Jg, where Jf= Forward Difference Approximation of
  // Jacobian of the function G defining the stiff term of the IVP_ODE
  //******************************************************       
  void Compute_G_Matrix_FD(const double t, double a,double *Y,LIS_MATRIX As)
  { 
    double * DY0=new double[neqn];
    double * DYF=new double[neqn];
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    LIS_INT  row=0, k;
    init_row_insertion(ptr, k, row);
    double Aij;
    const double uround=1.0e-15;
    double ysafe[neqn];
    //int nonzeros=0;
    G(t, Y, DY0);
    //for(int j=0;j<neqn;j++) ysafe=Y[j];
    cblas_dcopy(neqn, Y,1, ysafe, 1);
    for(int i=0;i<neqn;i++){
      for(int j=0;j<neqn;j++){
        const double delta=sqrt(uround*max(1.e-5,abs(ysafe[j])));
        Y[j]=ysafe[j]+delta;
        G(t, Y, DYF);
        Y[j]=ysafe[j];
        Aij=-a*(DYF[i]-DY0[i])/delta;
        if (i==j) Aij+=1.0;
        if (fabs(Aij)>uround){  
          new_entry(index, value, k, j, Aij);
          //nonzeros++;
        } 
      } 
      next_row(ptr, k, row);
    }  
    //cout<<"nonzeros="<<nonzeros<<endl;  
    delete[] DY0;
    delete[] DYF;    
  }

//******************************************************


///******************************************************
//Compute matrix  A=I-a*Jfeval, from the sparse matrix Jf
//******************************************************     
void Compute_Feval_Matrix(const double t, double a, LIS_MATRIX Jf,  LIS_MATRIX As)
  { 
    LIS_INT nnz=As->nnz;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    memcpy(ptr,Jf->ptr, (neqn+1)*sizeof(LIS_INT));
    memcpy(index,Jf->index, nnz*sizeof(LIS_INT));    
    cblas_dcopy(nnz, Jf->value,1, value, 1);
    cblas_dscal(nnz, -a, value, 1);
    for(int i=0;i<neqn;i++){   
      int j=ptr[i];
      while (index[j]!=i) {j++;}
      value[j]+=1.0;
    }
}
//******************************************************


///******************************************************
// Compute matrix Jfeval= Forward Difference Approximation of
// Jacobian of the function Feval=F+G defining the stiff term of the IVP_ODE
//******************************************************       
void Compute_Feval_Jacobian_FD(const double t,double *Y, LIS_MATRIX As){ 
  double * DY0=new double[neqn];
  double * DYFminus=new double[neqn];
  double * DYFplus=new double[neqn];
  LIS_INT *ptr=As->ptr;
  LIS_INT *index=As->index;
  LIS_SCALAR *value=As->value;
  LIS_INT  row=0, k;
  init_row_insertion(ptr, k, row);
  double Aij;
  const double uround=1.0e-15;
  double ysafe[neqn];
  feval(t, Y, DY0);
  cblas_dcopy(neqn, Y,1, ysafe, 1);
  for(int i=0;i<neqn;i++){
    for(int j=0;j<neqn;j++){
      const double delta=sqrt(uround*max(1.e-5,abs(ysafe[j])));
      Y[j]=ysafe[j]+delta;
      feval(t, Y, DYFplus);
      Y[j]=ysafe[j]-delta;
      feval(t, Y, DYFminus);
      Y[j]=ysafe[j];      
      Aij=(DYFplus[i]-DYFminus[i])/(2.0*delta);
      if (fabs(Aij)>uround){  
        new_entry(index, value, k, j, Aij);
      } 
    } 
    next_row(ptr, k, row);
  }  
  delete[] DY0;
  delete[] DYFplus;    
  delete[] DYFminus; 
}
 


///******************************************************
// Compute DY= Jfeval*Y where Jfeval= Forward Difference Approximation of
// Jacobian of the function Feval=F+G defining the stiff term of the IVP_ODE
//******************************************************       
void Jacobian_based_G(const double t, double * Y, double * DY,  LIS_MATRIX As) 
{  
  LIS_INT * jj0 = As->index;
	LIS_SCALAR * vv0 = As->value;
	for(int i=0;i<neqn;i++)
		{
			double t0  = 0;
			LIS_INT is = As->ptr[i];
			LIS_INT ie = As->ptr[i+1];
			for(int j=is;j<ie-0;j+=1)
			{
				const int j0 = jj0[j+0];
				t0 += vv0[j+ 0]*Y[j0];
			}
			DY[i] = t0;
		}
}
//******************************************************



};


#endif
