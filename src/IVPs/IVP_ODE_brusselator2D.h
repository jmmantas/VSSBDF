#ifndef BRUSSELATOR2D_H
#define BRUSSELATOR2D_H


#include "IVP_ODE.h"

using namespace std;


//*****************************************************************
// Class for the IVP-ODE  which represents 
// the Brusselator 2D model 
//*****************************************************************

class IVP_ODE_brusselator2D : public IVP_ODE {
private:
    int nx; // number of grid points at x dimension
    int ny; // number of grid points at y dimension
    double dtx; // Spatial step for both dimensions
    double dtx2; // dtx*dtx
    double DD;
    const double alpha = 0.002;
    const double A = 1.0;
    const double B = 3.4;


    // Auxiliary function f
    inline double f(const int i, const int j, const double t) {
        const double x = (i + 1) * dtx;    const double y = (j + 1) * dtx;
        const double xmxc = x - 0.3; const double ymyc = y - 0.5;
        const double r = 0.1;
        double result = ((xmxc * xmxc + ymyc * ymyc) <= r * r && t >= 1.1) ? 5.0 : 0.0;
        return result;
    }

    // Indexation function which maps 2D spatial coordinates (i,j)
    // to a 1D position in a vector
    inline int idx(const int i, const int j, const int k) { return 2 * (i * ny + j) + k; }

public:
    IVP_ODE_brusselator2D(const int nx_points);

    inline void init(double* Y) {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double x_i = (double)(i + 1) * dtx;
                double y_i = (double)(j + 1) * dtx;
                Y[idx(i, j, 0)] = 22 * y_i * pow(1 - y_i, 1.5);
                Y[idx(i, j, 1)] = 27 * x_i * pow(1 - x_i, 1.5);
            }
        }
    }
    //******************************************************
    //vector system function for the nonstiff term DY=G(t,Y)
    //******************************************************
    void G(const double t, const double* Y, double* DY);


    //******************************************************
    //ODE vector system function feval: DY=feval(t,Y)
    //******************************************************
    void F(const double t, const double* Y, double* DY);

    //******************************************************
    //Compute sparse matrix As=I-a*Jf, where Jf= Forward Diff 
    //Approx of the Jacobian of the function feval in (t,Y)
    //******************************************************  
    void Compute_Matrix_FD(const double t, double a, double* Y, LIS_MATRIX* As);

    //******************************************************
    //Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
    // in the IVP_ODE
    //****************************************************** 
    void Compute_G_Matrix_Exact(const double t, double A, double *Y,  LIS_MATRIX As);

    //******************************************************
    //Compute matrix Jf= Jacobian of the function 
    // feval defining the IVP_ODE
    //******************************************************       
    void Compute_Feval_Jacobian_exact(const double t, double* Y,
        LIS_MATRIX As);

};
//*****************************************************
// Constructor of the class IVP_ODE_brusselator2D
//*****************************************************
IVP_ODE_brusselator2D::IVP_ODE_brusselator2D(const int nx_points) {
    IVP_name = "Brusselator_2D";
    nx = nx_points;
    ny = nx;
    neqn = 2 * nx * ny;
    dtx = 1.0 / (nx + 1);
    dtx2 = dtx * dtx;
    DD = alpha / dtx2;
    nnz_G = 10 * nx * ny;
    nnz_FEVAL = nnz_G + neqn;
}

//***************************************************
//vector system function for the stiff term DY=G(t,Y)
//***************************************************
void IVP_ODE_brusselator2D::G(const double t, const double* Y, double* DY) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const double u_ij   = Y[idx(i, j, 0)];
            const double v_ij   = Y[idx(i, j, 1)];
            const double u_im1j = (i == 0)      ? Y[idx(nx - 1, j, 0)]      : Y[idx(i - 1, j, 0)];
            const double u_ip1j = (i == nx - 1) ? Y[idx(0, j     , 0)]      : Y[idx(i + 1, j, 0)];
            const double u_ijm1 = (j == 0)      ? Y[idx(i, ny - 1, 0)]      : Y[idx(i, j - 1, 0)];
            const double u_ijp1 = (j == ny - 1) ? Y[idx(i, 0     , 0)]      : Y[idx(i, j + 1, 0)];

            const double v_im1j = (i == 0)      ? Y[idx(nx - 1, j, 1)]      : Y[idx(i - 1, j, 1)];
            const double v_ip1j = (i == nx - 1) ? Y[idx(0, j     , 1)]      : Y[idx(i + 1, j, 1)];
            const double v_ijm1 = (j == 0)      ? Y[idx(i, ny - 1, 1)]      : Y[idx(i, j - 1, 1)];
            const double v_ijp1 = (j == ny - 1) ? Y[idx(i, 0     , 1)]      : Y[idx(i, j + 1, 1)];

            DY[idx(i, j, 0)] = DD * (u_im1j + u_ijm1 - 4.0 * u_ij + u_ip1j + u_ijp1);
            DY[idx(i, j, 1)] = DD * (v_im1j + v_ijm1 - 4.0 * v_ij + v_ip1j + v_ijp1);
        }
    }
}

//*******************************************************
//vector system function for the nonstiff term DY=F(t,Y)
//*******************************************************
void IVP_ODE_brusselator2D::F(const double t, const double* Y, double* DY) {
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            const double ui = Y[idx(i, j, 0)];
            const double vi = Y[idx(i, j, 1)];
            const double term = B * ui - ui * ui * vi;
            DY[idx(i, j, 0)] = A - term - ui + f(i, j, t);
            DY[idx(i, j, 1)] = term;
        }
    }
}
//******************************************************
//Compute matrix I-a*Jg where Jg= Jacobian of the Diffusive Term 
// in the IVP_ODE_brusselator2D
//******************************************************      
void IVP_ODE_brusselator2D::Compute_G_Matrix_Exact(const double t, double A, 
                                       double *Y,  LIS_MATRIX As)
{  
    LIS_INT k, row;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    const double Aii = 1.0+4.0 * A * DD;
    const double other = -A * DD;

    init_row_insertion(ptr, k, row);    
    for (int i = 0; i <= nx-1; i++) {    
        for (int j = 0; j <= ny-1; j++) {   
            for (int var = 0; var <= 1; var++) {
                if ( i==ny-1) {new_entry(index, value, k, idx(0,j,var)   , other);} // 2ny steps
                if ( i>0)     {new_entry(index, value, k, idx(i-1,j,var) , other);} // 2(nx-1)ny steps
                if (j==ny-1)  {new_entry(index, value, k, idx(i,0,var)   , other);} // 2 nx steps
                if (j>0)      {new_entry(index, value, k, idx(i,j-1,var) , other);} // 2nx (ny-1) steps
                              {new_entry(index, value, k, idx(i,j,var)   , Aii);  } // 2 nx ny steps
                if (j<ny-1)   {new_entry(index, value, k, idx(i,j+1,var) , other);} // 2nx (ny-1) steps
                if (j==0)     {new_entry(index, value, k, idx(i,ny-1,var), other);} // 2 nx steps 
                if (i<nx-1)   {new_entry(index, value, k, idx(i+1,j,var) , other);} // 2(nx-1)ny steps 
                if (i==0)     {new_entry(index, value, k, idx(nx-1,j,var), other);} // 2ny steps
                next_row(ptr, k, row);
            }
        }
    }

}


//******************************************************
//Compute matrix Jf= Jacobian of the function 
// feval defining the IVP_ODE
//******************************************************  
//******************************************************
// Compute the Jacobian matrix of the full right-hand side 
// of the Brusselator 2D system (diffusion + reaction)
//****************************************************** 
void IVP_ODE_brusselator2D::Compute_Feval_Jacobian_exact(const double t, double* Y, LIS_MATRIX As) {
   
    LIS_INT k, row;
    LIS_INT *ptr=As->ptr;
    LIS_INT *index=As->index;
    LIS_SCALAR *value=As->value;
    const double Aii_term = -4.0 * DD;
    const double other = DD;
    const double B1 = B + 1.0;
    double Aii[2];

    init_row_insertion(ptr, k, row);    
    for (int i = 0; i <= nx-1; i++) {    
        for (int j = 0; j <= ny-1; j++) {
            const double ui = Y[idx(i, j, 0)];
            const double vi = Y[idx(i, j, 1)];
            Aii[0]=Aii_term +2*ui*vi-B1;
            Aii[1]=Aii_term -  ui*ui;
            for (int var = 0; var <= 1; var++) {
                if ( i==ny-1) {new_entry(index, value, k, idx(0,j,var)   , other);}
                if ( i>0)     {new_entry(index, value, k, idx(i-1,j,var) , other);}
                if (j==ny-1)  {new_entry(index, value, k, idx(i,0,var)   , other);} 
                if (j>0)      {new_entry(index, value, k, idx(i,j-1,var) , other);}
                if (var==1)   {new_entry(index, value, k, idx(i,j,var)-1 ,B-2*ui*vi);}          
                              {new_entry(index, value, k, idx(i,j,var)   , Aii[var]);}
                if (var==0)   {new_entry(index, value, k, idx(i,j,var)+1 ,ui*ui);}           
                if (j<ny-1)   {new_entry(index, value, k, idx(i,j+1,var) , other);}
                if (j==0)     {new_entry(index, value, k, idx(i,ny-1,var), other);}
                if (i<nx-1)   {new_entry(index, value, k, idx(i+1,j,var) , other);} 
                if (i==0)     {new_entry(index, value, k, idx(nx-1,j,var), other);}
                next_row(ptr, k, row);
            }
        }
    }
}


#endif