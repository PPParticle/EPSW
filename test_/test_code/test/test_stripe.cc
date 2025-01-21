#include<iostream>
#include"ELLPACK_P.h"
#include<random>
#include<vector>
#define dim 10
using namespace std;
default_random_engine gen;
uniform_real_distribution<float> dis(-10,10);
int main(){
    Val *dense_matrix = new Val[dim*dim];
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            if(rand()&127!=0){// generate 0
                dense_matrix[i*dim+j]=dis(gen);
            }
            else{
                dense_matrix[i*dim+j]=0.0;
            }   
        }
    }

    ELLPACK_P M(dim, dim);
    M.ELLPACKinitWithDenseMatrix(dense_matrix, dim, dim);
    M.ELLPACK_P_vreorder();
    M.print_ELLPACK_P();
    M.stripe();
    M.print_strips();


    return 0;
}