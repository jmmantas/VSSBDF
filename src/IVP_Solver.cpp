/*
This C++ driver program run the experiments associated with a given IVP model using a particular VSSBDF method  
and generate the numerical errors for both splitting approaches, the execution times, and the 
solutions at the final simulation times. 

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

#include "IVP_ODE_advdiff1d.h"
#include "IVP_ODE_combustion.h"
#include "IVP_ODE_brusselator.h"
#include "IVP_ODE_brusselator2D.h"
#include "IVP_ODE_stiff_brusselator.h"
#include "IVP_ODE_heat.h"
#include "IVP_ODE_cusp.h"
#include "IVP_ODE_simpleadvdiff1d.h"

#include "SBDF_Solver.h"
using namespace std;
using namespace std::chrono;

typedef std::chrono::time_point<high_resolution_clock, nanoseconds> time_ns;



//**************************************************************************
//AUXILIARY  FUNCTIONS
//**************************************************************************


//***************************************************
// Print vector values on standard output
void print_vector(const int neqn, const double* Y, string str)
//***************************************************                                              
{
    for (int i = 0; i < neqn; i++)
    {
        cout << str<<"[" << i << "]=" << Y[i] << endl;
    }
}

//**************************************************************
// Check important errors between sparse matrices A1 and A2
void compare_matrix_csr(LIS_MATRIX A1, LIS_MATRIX A2)
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
//**************************************************************



//******************************************************************
  // Init CSR Sparse neqn x neqn LIS matrix
  inline void init_CSR_LIS_matrix(const int neqn, const int nnz, 
                            const int splitting_type, LIS_MATRIX *A) 
//******************************************************************
  {
    LIS_INT * ptr_A, * index_A;
    LIS_SCALAR* value_A;
    lis_matrix_create(0, A);
    lis_matrix_set_size(*A, 0, neqn);
    lis_matrix_malloc_csr(neqn, nnz, &ptr_A, &index_A, &value_A);
    lis_matrix_set_csr(nnz,ptr_A,index_A,value_A,*A);
    lis_matrix_assemble(*A);
  }
//******************************************************************




//***************************************************
//1st, 2nd and 3rd order Runge-Kutta Time Integrator
//***************************************************
void RK(const int nstages, IVP_ODE* IVP, const double t0,
    const double tf, const double h0, const double* Y0, double* Y1)
//***************************************************
{
    const int neqn = IVP->get_num_ODEs();

    double* Y0_tmp = new double[neqn];
    double* DY1 = new double[neqn];
    double* DY2 = new double[neqn];
    double* DY3 = new double[neqn];
    double* DY4 = new double[neqn];

    double t = t0;
    bool end = false;
    double h = h0;
    int N = 0;
    cblas_dcopy(neqn, Y0, 1, Y0_tmp, 1);

    while (!end) {
        if ((tf - t) <= h) { h = tf - t; end = true; }
        // DY1=F(t,Y0)  
        IVP->feval(t, Y0_tmp, DY1);
        switch (nstages)
        {
        case 1:
            cblas_daxpy(neqn, h, DY1, 1, Y0_tmp, 1);
            break;
        case 2:
            cblas_dscal(neqn, 0.5 * h, DY1, 1);
            cblas_daxpy(neqn, 1, Y0_tmp, 1, DY1, 1);
            IVP->feval(t + 0.5 * h, DY1, Y1);
            cblas_daxpy(neqn, h, Y1, 1, Y0_tmp, 1);
            break;
        case 3:
            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5 * h, DY1, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY2);

            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, -h, DY1, 1, Y1, 1);
            cblas_daxpy(neqn, 2 * h, DY2, 1, Y1, 1);
            IVP->feval(t + h, Y1, DY3);

            cblas_daxpy(neqn, h / 6    , DY1, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, 4 * h / 6, DY2, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 6    , DY3, 1, Y0_tmp, 1);
            break;
        // RK4 
        case 4:
            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5*h, DY1, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY2);

            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, 0.5*h, DY2, 1, Y1, 1);
            IVP->feval(t + 0.5 * h, Y1, DY3);

            cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);
            cblas_daxpy(neqn, h, DY3, 1, Y1, 1);
            IVP->feval(t + h, Y1, DY4);

            cblas_daxpy(neqn, h / 6, DY1, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 3, DY2, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 3, DY3, 1, Y0_tmp, 1);
            cblas_daxpy(neqn, h / 6, DY4, 1, Y0_tmp, 1);
            break;
        }
        N++;
        t = t + h;
    }
    cblas_dcopy(neqn, Y0_tmp, 1, Y1, 1);

    delete[] DY1;
    delete[] DY2;
    delete[] DY3;
    delete[] DY4;
    delete[] Y0_tmp;
}
//***************************************************


//**************************************************************************
// Initialize  intermediate vectors starting the multistep integration
//**************************************************************************
void Init_Vectors(const int order, const int neqn,  const double t0, 
                  const double h,  IVP_ODE* IVP, double ** Y_init, double ** Yf_init) 
{
    const int RK_order=3, idx=order-1;
    IVP->init(Y_init[0]); //First  step initial value

    // Approximation of the intermediate step initial values 
    // using RK4 and very small stepsize
	const double h_RK=1.0e-8;
    // Initialize initial condition vectors for coarse time integration
    double t2=t0;           
    for (int i = 1; i < order; i++){  
        RK(RK_order, IVP, t2, t2+h, h_RK, Y_init[i-1], Y_init[i]);
        t2+=h;
    } 
    // Set up initial condition vectors for 
    // fine time integration (half of the time step size)
    const double h_over2=h/2;
    if (order==2) {
        t2=t0;
        RK(RK_order, IVP, t2, t2+h_over2, h_RK, Y_init[0], Yf_init[0]);
    }
    else if (order==3) {
        cblas_dcopy (neqn, Y_init[1],1, Yf_init[0] , 1);
        t2=t0+h;
        RK(RK_order, IVP, t2, t2+h_over2, h_RK, Yf_init[0], Yf_init[1]);            
    } else if (order==4) {
        t2=t0+h;
        RK(RK_order, IVP, t2, t2+h_over2, h_RK, Y_init[1], Yf_init[0]);
        cblas_dcopy (neqn, Y_init[2],1, Yf_init[1] , 1);
        t2=t0+2*h;
        RK(RK_order, IVP, t2, t2+h_over2, h_RK, Yf_init[1], Yf_init[2]);
    }
    cblas_dcopy (neqn, Y_init[idx],1, Yf_init[idx] , 1);
}

//**************************************************************************
// Free vectors in dynamic memory	
//**************************************************************************
void Free_vectors(const int order, double **Y_init, double ** Yf_init,
                   double * Y1_RK,double *Y1_variable)
{

    for (int i = 0; i < order; i++) 
        delete[] Y_init[i];
    for (int i = 0; i < (order+2); i++) 
        delete[] Yf_init[i];  
    delete[] Y1_variable;
    delete[] Y1_RK;
}
//**************************************************************************


//***************************************************
// MAIN PROGRAM
int main(int argc, char** argv)
//***************************************************
{

    //const double EPSTOL=1.0e-12;
    // Initializa LIS environment
    lis_initialize(&argc, &argv);
    const int num_args=9;
    // Check the number of parameters
    if (argc < num_args) {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] 
        << " <problem id.> <conv. order> <Neqn(N)>  <tf(final time)> <RK stepsize> "
        << " init. SBDF stepsize>  <Num. experiments> <init_tolerance>" << endl<<endl;

        cerr << "problem.id= 0: 1D Advection-Diffusion" <<endl;
        cerr << "problem.id= 1: Combustion" <<endl;
        cerr << "problem.id= 2: 1D Brusselator" <<endl;
        cerr << "problem.id= 3: 2D Heat" <<endl;
        cerr << "problem.id= 4: CUSP" <<endl;
        cerr << "problem.id= 5: Simple 1D Advection-Diffusion" <<endl;
        cerr << "problem.id= 6: 2D Brusselator" << endl;
        cerr << "problem.id= 7: Stiff Brusselator" << endl;

        return 1;
    }
    int problem_id=atoi(argv[1]);
    const int RK_order=4;
    const unsigned order = atoi(argv[2]);
    const int Num_points = atoi(argv[3]);
    const double tf = atof(argv[4]);
    //Time integration step for RK3 and SBDF1 methods
    double h_RK = atof(argv[5]); //Time integration stepsize for RK4 method
    double h_imex = atof(argv[6]); //Time integration stepsize for SBDF method
    const int num_tests=atoi(argv[7]); // Number of experiments
    double start_tol = atof(argv[8]); //Tolerance for the first experiment


    // Declaration of the C++ objects representing 
    // the IVP-ODE problems to be solved 
    IVP_ODE_advdiff1d   IVP_advdiff1d(Num_points);
    IVP_ODE_combustion  IVP_combustion(Num_points);
    IVP_ODE_brusselator IVP_brusselator(Num_points);
    IVP_ODE_heat        IVP_heat(Num_points,Num_points/7); 
    IVP_ODE_cusp        IVP_cusp(Num_points);
    IVP_ODE_simpleadvdiff1d IVP_simpleadvdiff1d (Num_points);
    IVP_ODE_brusselator2D IVP_brusselator2D(Num_points);
    IVP_ODE_stiff_brusselator IVP_stiff_brusselator(Num_points);

    const int num_IVPs = 8;

    IVP_ODE* IVP_list[num_IVPs];

    IVP_list[0]= &IVP_advdiff1d;
    IVP_list[1]= &IVP_combustion;
    IVP_list[2]= &IVP_brusselator;
    IVP_list[3]= &IVP_heat;
    IVP_list[4]= &IVP_cusp;
    IVP_list[5]= &IVP_simpleadvdiff1d;
    IVP_list[6]= &IVP_brusselator2D;
    IVP_list[7]= &IVP_stiff_brusselator;
    
    //Time-step conditions
    const double t0 = 0.0;

    ///////////////////////

    /// This part to save some results in txt file

    ofstream file1("TTime_P_S0_order" + to_string(order) + ".txt", ios_base::app);
    ofstream file2("Error_P_S0_order" + to_string(order) + ".txt");
    ofstream file3("TTime_J_S0_order" + to_string(order) + ".txt", ios_base::app);
    ofstream file4("Error_J_S0_order" + to_string(order) + ".txt");


    //ofstream file5("Time_P_S1_order" + to_string(order) + ".txt", ios_base::app);
    //ofstream file6("Error_P_S1_order" + to_string(order) + ".txt");
    //ofstream file7("Time_J_S1_order" + to_string(order) + ".txt", ios_base::app);
    //ofstream file8("Error_J_S1_order" + to_string(order) + ".txt");
    // Check if the file streams are open before attempting to write to files
    if (!file1.is_open() || !file2.is_open() || !file3.is_open() || !file4.is_open()) {
        std::cerr << "Error: One or more files could not be opened." << std::endl;
        return 1; // Exit program, error opening one of the files
    }


    ////////////
    
    // Get the number of ODEs in the test problem
    const int neqn = IVP_list[problem_id]->get_num_ODEs();
    cout <<"NEQN="<<neqn<<endl;
    cout << "................................................" << endl;

    cout<<endl<<    IVP_list[problem_id]->get_name()<<" IVP with Neqn = "<<neqn<<endl;
    cout << "................................................" << endl;

    //Set up Initial condition Vectors: Y[0] in t=t0,
    // Y_init[1] in t=t0+h, ... Y_init[order-1] in t=t0+(order-1)*h
    double *Y_init[order];
    for (unsigned i = 0; i < order; i++) {
            Y_init[i] = new double[neqn];
    }
    //Set up Initial condition Vectors for the fine approximation: 
    // Yf_init[0] in t=t0, Yf_init[1] in t=t0+h/2, Yf_init[2] in t+h, ...
    // Setup initial condition vectors and final vectors, Yf_init, for the fine time steps
    double *Yf_init[order+2];
    for (unsigned int i = 0; i < (order+2); i++) {
        Yf_init[i] = new double[neqn];
    }
        
    //Declare vector Y1_variable to store the numerical approximation 
    // to the solution	 
    double* Y1_variable = new double[neqn];

    //Declare vector Y1 to store the numerical approximation 
    // to the solution with RK solver	 
    double* Y1_RK = new double[neqn];

    cout.precision(3);

    //**********************************************
    double runtime_RK;
    cout << "................................................" << endl;
    // Numerical solution with the Explicit Runge-Kutta method 
    cout << "<<<< Explicit Runge-Kutta order "<<RK_order
         <<" ...   Time Stepsize h=" << h_RK << " >>>> " << endl;
    cout << "<<<< N= " << neqn << "...t0=" << t0 << "...tf= " << tf << endl;

    //First  step initial value    
    IVP_list[problem_id]->init(Y_init[0]);

    //print_vector(neqn, Y[0], "Y0");

    // Execution of the RK Solver to approximate accurately the solution
    time_ns start,end;
    start = high_resolution_clock::now();
    RK(RK_order, IVP_list[problem_id], t0, tf, h_RK, Y_init[0], Y1_RK);
    //IVP_list[problem_id]->exact_solution (tf , Y1_RK);    
    end = high_resolution_clock::now();
    runtime_RK = duration_cast<nanoseconds>(end - start).count() * 1e-9;
    cout << "................................................" << endl;
    cout << "Runtime RK=  " << runtime_RK << endl;
    cout << "................................................" << endl;
    //print_vector(neqn, Y1_RK, "Y1_RK");
        
    // Numerical solution with the semiimplicit SBDF-order Method
    cout << "<<<< SBDF" << order << "...   Time Step h=" << h_imex << " >>>>" << endl;

    // Executions of SBDF solver
    double error[2][num_tests], runtime[2][num_tests];
    int n_steps[2][num_tests], n_isteps[2][num_tests];
    double h0 = h_imex;
    time_ns start_SBDF, end_SBDF;
    double tolerance= start_tol;       
    // Approximation of the intermediate step initial values 
    // using RK4 and very small stepsize
    Init_Vectors(order,neqn,  t0, h0,  IVP_list[problem_id], Y_init, Yf_init);   
    
    // Init SBDF solvers for both type of splittings (0 and 1)
    SBDF_Solver * SBDF_IVP[2];
    for (int j = 0; j <= 1; j++){            
      SBDF_IVP[j]=new SBDF_Solver(order, IVP_list[problem_id], j);
    }
    
    //******************************************************************************
    // Perform num_tests x 2 experiments using splitting type=0 and splitting type=1
    //******************************************************************************
    for (int i = 0; i < num_tests; i++){   
    //******************************************************************************
        for (int j = 0; j <= 1; j++){   
            start_SBDF = high_resolution_clock::now();
            SBDF_IVP[j]->Adaptive_dt_SBDF_Integrate0(t0, tf, h0, Y_init, Yf_init, Y1_variable,
                                          tolerance, &(n_steps[j][i]), &(n_isteps[j][i]));
            end_SBDF = high_resolution_clock::now();
            //SBDF_IVP.Store_result(tf, Y1_variable,error_tolerance);
            runtime[j][i] = duration_cast<nanoseconds>(end_SBDF - start_SBDF).count() * 1e-9;
            cblas_daxpy(neqn, -1.0, Y1_RK, 1, Y1_variable, 1);
            error[j][i] = cblas_dnrm2(neqn, Y1_variable, 1)/pow(neqn, 0.5);
            cout << "TOL=  "<<tolerance<<".....Splitting="<<j
                             <<"....NSTEPS= "   <<n_steps[j][i]
                             <<"......RUNTIME= "<< runtime[j][i]
                             <<"......ERROR = " << error[j][i]<< endl;  
            if (j == 0) {
                file1 << runtime[j][i] << endl;
                file2 << error[j][i] << endl;
            }
            else {
                file3 << runtime[j][i] << endl;
                file4 << error[j][i] << endl;
            }
        }
        tolerance=tolerance/10;
      
    }
    file1.close();
    file2.close();
    file3.close();
    file4.close();
    cout << endl << endl;
    cout << IVP_list[problem_id]->get_name()<<" with Neqn = "
             <<neqn<<"......VSSSBDF"<<order<<" .....  ";
    for (int j = 0; j < num_args; j++){
        cout<<argv[j]<<" ";
    }     
    cout<<endl<<endl;
    for (int i = 0; i < num_tests; i++){ 
        cout << endl<<"***** TOL=" <<setw(4) <<start_tol*pow(1.0e-1,i)<<"******"<<endl; 
        for (int splitting = 0; splitting <=1; splitting++){
            cout<<"ERR["<<splitting<<"]("<<n_steps[splitting][i]<<"," 
                 <<n_isteps[splitting][i]<<")= " 
                 << setw(4)<< error[splitting][i];
            cout <<"..TIME["<<splitting<<"]= "<<  setw(4)<<runtime[splitting][i];
            if (splitting==0) 
                cout<<"-----";
            else 
            cout<<" ---->  ERR-RAT= "<<(error[1][i])*100/error[0][i]<<"%  "
                                      <<"TIME-RAT= "<<(runtime[1][i])*100/runtime[0][i]<<"%"<<endl;         
        } 
    }

    // Free vectors in dynamic memory	
    Free_vectors(order, Y_init, Yf_init,Y1_RK,Y1_variable);
    for (int j = 0; j <= 1; j++){            
      delete SBDF_IVP[j];
    }

    //Finalize LIS environment 
    lis_finalize();

    cout << "Done!" << endl;
}
