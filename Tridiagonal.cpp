#include <iostream>
#include<vector>
#include <fstream>
#include <cstdlib>
#include<cmath>
#include<ctime>
#define MATRIX_SIZE 100;
#define ITERATION_MAX 2000000;
using namespace std;
class TridiagMatrix{
   
    public :
        int size ;
        // Constructor
        TridiagMatrix(int _size)
        {
            size = _size;
        }
   
    double getMatrixEntry(int i , int j );
    // Test Function to Print
    void PrintMatrix();
};
double TridiagMatrix::getMatrixEntry(int i , int j)
{
   double value  = 0;
    if ( i == j)
        value = 2.0;
    else if ( abs( i - j ) == 1)
        value = -1.0;
    else 
        value = 0;
   
    if(i >= size || j >= size ) cout << " I, J exceeds Matrix size" <<endl;
   return value ;
}


void TridiagMatrix::PrintMatrix()
{
    cout << " Sixe : " << size <<endl;
    for ( int i = 0 ;i < size ; i++)
    {
        for ( int j  = 0 ; j < size ; j++)
        {
            cout << getMatrixEntry(i,j) << "  " ;
        }
        cout << endl;
    }
}

void GaussSeidalWeightedTridiagonal(double* xo, double* b , TridiagMatrix* matrix, double*& sol ,int eig_mode,double omega )
{
    // NOte : the Second For loop must permit only two members inside the Computation loop
    // NOte : fr the 1st and Last element , the NUmber of elements enetering should be 1
    // Sol(1) =  1/a_11 * ( rhs(1)  - (x2*a_12 + x3*a_13))
    int mat_size = matrix->size;
   
    std::string filename = "w_gs";
    std::string num = std::to_string(eig_mode);
    std::string c;
    if(omega == 1.2 )  c = "01";
    if(omega == 1.4 ) c = "02";
    if(omega == 1.6 ) c = "03";
    if(omega == 1.8 )  c = "04";
    if(omega == 1.9 )  c = "05";
    if(omega == 1.999) c = "06";
    filename = filename + "_" +  to_string(mat_size) + "_" + num + "_" + c;

    //cout<< " FileName : " <<filename<<endl;
    int last_index = mat_size -1;
    int last_index_prev = mat_size - 2;
    double first_val = 0;
    double second_val = 0;
    double diag_val = 0.;
    unsigned int iterationCount = 0;
    ofstream file;
    file.open(filename);
    file.precision(6);
    file << "#  "<<"Iteration" << "," <<"Error" <<endl;
    const clock_t begin_time = clock();
    unsigned int ITER_MAX = ITERATION_MAX;
    double residual = 1000;
    double temp = 0;
    
        double residual_inf = 0;
    do
    {  
        residual_inf = 0;
        // Compute the Solution for 1st and last row
        temp = 0.5 * ( b[0]  - (-1.0 * xo[1]));
        sol[0] =  ( temp * omega )  + (1.0 - omega)*xo[0];
        // For the other rows
        for ( int i = 1 ; i < mat_size-1 ; i++)
        {
            temp  =  (0.5) *  ( b[i] + ( (sol[i-1]*1.0) + (xo[i+1]*1.0) ) );  
            sol[i] =  ( temp * omega ) + (1.0 - omega)*xo[i] ; 
        }

        temp = 0.5 * ( b[last_index] - (-1.0*sol[last_index_prev]) );
        sol[last_index] = ( temp * omega )  + (1.0 - omega)*xo[last_index];
        // Compute Residual
        //Compute AX

        // residual = 0;
        // residual += pow( b[0] - (2.0*sol[0] - 1.0*sol[1]) , 2 ) ;
        // residual += pow (b[last_index] - (2.0*sol[last_index] -1.0*sol[last_index_prev]) ,2) ;
        // for ( int i = 1 ; i < mat_size - 1 ; i++)
        //     residual +=  pow (b[i] - ( -1.0*sol[i-1] - 1.0*sol[i+1] + 2.0*sol[i] ) ,2);

        // residual = sqrt(residual);
        
        for ( int k = 0 ; k < mat_size ; k++)
            if(fabs(sol[k]) > residual_inf )
                residual_inf = fabs(sol[k]);
        
        //cout << " Res INf : " << residual_inf <<endl;
        // Transfer the Old Solution to new Solution
        for (int i = 0 ; i < mat_size ; i++)
            xo[i] = sol[i];

        file << iterationCount <<"," <<residual_inf<<endl;
        //cout << iterationCount << "\t" << residual <<endl;
        iterationCount ++;
    }
    while ( (iterationCount < ITER_MAX) && residual_inf>1e-6 );
    const clock_t end_time = clock();
   
    file.close();
    if(iterationCount < ITER_MAX  ){
        cout << "    Jacobi Iteration Completed **SUCCESSFULL***" <<endl;
        cout << "    TIME TAKEN : " <<  ( end_time - begin_time ) / CLOCKS_PER_SEC <<endl;
        cout << "    RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    else {
       // cout << "   JACOBI ITERATION HAS ***NOT*** CONVERGED " <<endl;
        cout << "   Weight GAUSS Seidal RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    //system("gnuplot -e \"plot 'jacobi.txt' using 1:2 with linespoints; pause -1 ; exit\" ");
}

void JacobiWeightedTridiagonal(double* xo, double* b , TridiagMatrix* matrix, double*& sol ,int eig_mode,double omega )
{
    // NOte : the Second For loop must permit only two members inside the Computation loop
    // NOte : fr the 1st and Last element , the NUmber of elements enetering should be 1
    // Sol(1) =  1/a_11 * ( rhs(1)  - (x2*a_12 + x3*a_13))
    int mat_size = matrix->size;
   
    std::string filename = "w_jacobi";
    std::string num = std::to_string(eig_mode);
    filename = filename + "_" +  to_string(mat_size) + "_" + num;
    //cout<< " FileName : " <<filename<<endl;
    int last_index = mat_size -1;
    int last_index_prev = mat_size - 2;
    double first_val = 0;
    double second_val = 0;
    double diag_val = 0.;
    unsigned int iterationCount = 0;
    ofstream file;
    file.open(filename);
    file.precision(6);
    file << "#  "<<"Iteration" << "," <<"Error" <<endl;
    const clock_t begin_time = clock();
    unsigned int ITER_MAX = ITERATION_MAX;
    double residual = 1000;
    double temp = 0;
    do
    {  
        // Compute the Solution for 1st and last row
        temp = 0.5 * ( b[0]  - (-1.0 * xo[1]));
        sol[0] =  ( temp * omega )  + (1.0 - omega)*xo[0];
        temp = 0.5 * ( b[last_index] - (-1.0*xo[last_index_prev]) );
        sol[last_index] = ( temp * omega )  + (1.0 - omega)*xo[last_index];
        // For the other rows
        for ( int i = 1 ; i < mat_size-1 ; i++)
        {
            temp  =  (0.5) *  ( b[i] + ( (xo[i-1]*1.0) + (xo[i+1]*1.0) ) );  
            sol[i] =  ( temp * omega ) + (1.0 - omega)*xo[i] ; 
        }
        // Compute Residual
        //Compute AX

        // residual = 0;
        // residual += pow( b[0] - (2.0*sol[0] - 1.0*sol[1]) , 2 ) ;
        // residual += pow (b[last_index] - (2.0*sol[last_index] -1.0*sol[last_index_prev]) ,2) ;
        // for ( int i = 1 ; i < mat_size - 1 ; i++)
        //     residual +=  pow (b[i] - ( -1.0*sol[i-1] - 1.0*sol[i+1] + 2.0*sol[i] ) ,2);

        // residual = sqrt(residual);
        
        double residual_inf = 0;
        for ( int k = 0 ; k < mat_size ; k++)
            if(fabs(sol[k]) > residual_inf )
                residual_inf = fabs(sol[k]);
        
        //cout << " Res INf : " << residual_inf <<endl;
        // Transfer the Old Solution to new Solution
        for (int i = 0 ; i < mat_size ; i++)
            xo[i] = sol[i];

        file << iterationCount <<"," <<residual_inf<<endl;
        //cout << iterationCount << "\t" << residual <<endl;
        iterationCount ++;
    }
    while ( (iterationCount < ITER_MAX)  );
    const clock_t end_time = clock();
   
    file.close();
    if(iterationCount < ITER_MAX  ){
        cout << "    Jacobi Iteration Completed **SUCCESSFULL***" <<endl;
        cout << "    TIME TAKEN : " <<  ( end_time - begin_time ) / CLOCKS_PER_SEC <<endl;
        cout << "    RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    else {
       // cout << "   JACOBI ITERATION HAS ***NOT*** CONVERGED " <<endl;
        cout << "   Weight jacobi RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    //system("gnuplot -e \"plot 'jacobi.txt' using 1:2 with linespoints; pause -1 ; exit\" ");
}


void JacobiTridiagonal(double* xo, double* b , TridiagMatrix* matrix, double*& sol ,int eig_mode )
{
    // NOte : the Second For loop must permit only two members inside the Computation loop
    // NOte : fr the 1st and Last element , the NUmber of elements enetering should be 1
    // Sol(1) =  1/a_11 * ( rhs(1)  - (x2*a_12 + x3*a_13))
    int mat_size = matrix->size;
   
    std::string filename = "jacobi";
    std::string num = std::to_string(eig_mode);
    filename = filename + "_" +  to_string(mat_size) + "_" + num;
    //cout<< " FileName : " <<filename<<endl;
    int last_index = mat_size - 1 ;
    int last_index_prev = mat_size - 2;
    double first_val = 0; double second_val = 0; double diag_val = 0.;
    unsigned int iterationCount = 0;
    ofstream file;
    file.open(filename);
    file.precision(6);
    file << "#  "<<"Iteration" << "," <<"Error" <<endl;
    const clock_t begin_time = clock();
    unsigned int ITER_MAX = ITERATION_MAX;
    double residual = 1000;
    do
    {  
        // Compute the Solution for 1st and last row
        sol[0] = 0.5 * ( b[0]  - (-1.0 * xo[1]));
        sol[last_index] = 0.5 * ( b[last_index] - (-1.0*xo[last_index_prev]) );
        // For the other rows
        for ( int i = 1 ; i < mat_size-1 ; i++)
        {
            sol[i] =  (0.5) *  ( b[i] + ( (xo[i-1]*1.0) + (xo[i+1]*1.0) ) );  
        }
        
        
        // Compute Residual
        //Compute AX
        
        // residual = 0;
        // residual += pow( b[0] - (2.0*sol[0] - 1.0*sol[1]) , 2 ) ;
        // residual += pow (b[last_index] - (2.0*sol[last_index] -1.0*sol[last_index_prev]) ,2) ;
        // for ( int i = 1 ; i < mat_size - 1 ; i++)
        //     residual +=  pow (b[i] - ( -1.0*sol[i-1] - 1.0*sol[i+1] + 2.0*sol[i] ) ,2);

        // residual = sqrt(residual);



        double residual_inf = 0;
        for ( int k = 0 ; k < mat_size ; k++)
            if(fabs(sol[k]) > residual_inf )
                residual_inf = fabs(sol[k]);
        

        
        // cout << " Res INf : " << residual_inf <<endl;
        // Transfer the Old Solution to new Solution
        for (int i = 0 ; i < mat_size ; i++)
            xo[i] = sol[i];

        file << iterationCount <<"," <<residual_inf<<endl;
        //cout << iterationCount << "\t" << residual <<endl;
        iterationCount ++;
    }
    while ( (iterationCount < ITER_MAX)  );
    const clock_t end_time = clock();
   
    file.close();
    if(iterationCount < ITER_MAX  ){
        cout << "    Jacobi Iteration Completed **SUCCESSFULL***" <<endl;
        cout << "    TIME TAKEN : " <<  ( end_time - begin_time ) / CLOCKS_PER_SEC <<endl;
        cout << "    RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    else {
       // cout << "   JACOBI ITERATION HAS ***NOT*** CONVERGED " <<endl;
        cout << "   Jacobi RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    //system("gnuplot -e \"plot 'jacobi.txt' using 1:2 with linespoints; pause -1 ; exit\" ");
}


void gaussSeidalTridiagonal(double* xo, double* b , TridiagMatrix* matrix, double*& sol ,int eig_mode )
{
    // NOte : the Second For loop must permit only two members inside the Computation loop
    // NOte : fr the 1st and Last element , the NUmber of elements enetering should be 1
    // Sol(1) =  1/a_11 * ( rhs(1)  - (x2*a_12 + x3*a_13))
    int mat_size = matrix->size;
   
    std::string filename = "gs";
    std::string num = std::to_string(eig_mode);
    filename = filename + "_" +  to_string(mat_size) + "_" + num;
    //cout<< " FileName : " <<filename<<endl;
    int last_index = mat_size - 1 ;
    int last_index_prev = mat_size - 2;
    double first_val = 0; double second_val = 0; double diag_val = 0.;
    unsigned int iterationCount = 0;
    ofstream file;
    file.open(filename);
    file.precision(6);
    file << "#  "<<"Iteration" << "," <<"Error" <<endl;
    const clock_t begin_time = clock();
    unsigned int ITER_MAX = ITERATION_MAX;
    double residual = 1000;
        double residual_inf = 0;
    do
    {  
        residual_inf = 0;
        // Compute the Solution for 1st and last row
        sol[0] = 0.5 * ( b[0]  - (-1.0 * xo[1]));

        // For the other rows
        for ( int i = 1 ; i < mat_size-1 ; i++ )
        {
            sol[i] =  (0.5) *  ( b[i] + ( (sol[i-1]*1.0) + (xo[i+1]*1.0) ) );  
        }
        sol[last_index] = 0.5 * ( b[last_index] - (-1.0*sol[last_index_prev]) );
        
        // Compute Residual
        //Compute AX
        
        // residual = 0;
        // residual += pow( b[0] - (2.0*sol[0] - 1.0*sol[1]) , 2 ) ;
        // residual += pow (b[last_index] - (2.0*sol[last_index] -1.0*sol[last_index_prev]) ,2) ;
        // for ( int i = 1 ; i < mat_size - 1 ; i++)
        //     residual +=  pow (b[i] - ( -1.0*sol[i-1] - 1.0*sol[i+1] + 2.0*sol[i] ) ,2);

        // residual = sqrt(residual);

        for ( int k = 0 ; k < mat_size ; k++)
            if(fabs(sol[k]) > residual_inf )
                residual_inf = fabs(sol[k]);
        

        
        // cout << " Res INf : " << residual_inf <<endl;
        // Transfer the Old Solution to new Solution
        for (int i = 0 ; i < mat_size ; i++)
            xo[i] = sol[i];

        file << iterationCount <<"," <<residual_inf<<endl;
        //cout << iterationCount << "\t" << residual <<endl;
        iterationCount ++;
    }
    while ( (iterationCount < ITER_MAX)  && fabs(residual_inf) > 1e-3);
    const clock_t end_time = clock();
   
    file.close();
    if(iterationCount < ITER_MAX  ){
        cout << "    gauss Seidal Iteration Completed **SUCCESSFULL***" <<endl;
        cout << "    TIME TAKEN : " <<  ( end_time - begin_time ) / CLOCKS_PER_SEC <<endl;
        cout << "    RESIDUAL : " <<fabs(residual_inf) << "  Iteration Taken : " << iterationCount <<endl;
    }
    else {
       // cout << "   JACOBI ITERATION HAS ***NOT*** CONVERGED " <<endl;
        cout << "   Gauss Seidal RESIDUAL : " <<residual << "  Iteration Taken : " << iterationCount <<endl;
    }
    //system("gnuplot -e \"plot 'jacobi.txt' using 1:2 with linespoints; pause -1 ; exit\" ");
}



int main (int argc, char** argv)
{
    // create a Object for the Matrix
    //int size = MATRIX_SIZE;
    int size = stoi(argv[1]);
    int eig_mode = stoi(argv[2]);
    double omega = stod(argv[3]);
    int solver = stoi(argv[4]);
    

    cout << " SIZE : "<< size << "   EIG MODE : "<< eig_mode << " OMEGA : "<< omega <<endl;


    if(argc < 4 ){
        cout << " INSUFFICIENT ARGUMENTS " <<endl;
        exit(0);
    }


    TridiagMatrix* matrix = new TridiagMatrix(size);
    // Create The Solution Array b
    double* b =  new double[size]();
    //Create The iniial Solution x0
    double* x0 =  new double[size]();
    // Fill The Vectors with 1st Eigen Value
    double norm = 0;
    for ( int i = 0 ; i < size ; i++)
    {
        double temp = cos(eig_mode*3.1415926535/(size+1));
        temp *=temp;
        temp = pow(temp,(i+1)/2);
        x0[i] = 1.0* ( sin(eig_mode * (i+1) * 3.1415926535 / (size+1 ) ) )   ;
        norm += x0[i]*x0[i];
    }

    // Not Normalising the Solution

    // norm = sqrt(norm);
    // // // NOrmalise the Vector
    // for ( int i = 0 ; i < size ; i++){
    //     x0[i] = x0[i] / norm;
    // }


    double* sol = new double[size]();
   
    // Call for JAcobi Solver
    if(solver == 1) JacobiTridiagonal(x0,b,matrix,sol,eig_mode);
    else if(solver==2) JacobiWeightedTridiagonal(x0,b,matrix,sol,eig_mode,omega);
    else if(solver==3) GaussSeidalWeightedTridiagonal(x0,b,matrix,sol,eig_mode,omega);
    else if(solver==4) gaussSeidalTridiagonal(x0,b,matrix,sol,eig_mode);
    else {cout<<" NO SOLVER" <<endl;}
    
 
    return 0;
}