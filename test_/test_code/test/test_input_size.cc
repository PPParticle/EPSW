#include<iostream>
#include<riscv-vector.h>
#include"ELLPACK_P_0.7.h"
#include<random>
#include"utilis.h"
#include <chrono>

default_random_engine gen;
uniform_real_distribution<float> dis(-10,10);
using namespace std;

// inner dataflow scalar kernel
void scalar(Val *a, const unsigned int ai
                        Val *b, const unsigned int bi 
                        Val *c, const unsigned int ci){
    
    for(int i=0; i<ai; i++)
        for(int j=0; j<bi; j++){
            c[j]=a[i]*b[j];
        }
}

void rvv(Val *a, const unsigned int ai, Val *b, const unsigned int bi, Val *c, const unsigned int ci) {
    size_t vl;
    float32xm8_t vtmp, vb, vc;
    vtmp = vfmvvf_float32xm8 (0.0, vsetvli_max(RVV_E32, RVV_M8));
    float *tmp_a = a;
    for (int j=bi; j > 0; j-=vl) {
        float *tmp_b = b;
        float *tmp_c = c;
        for(int i = 0; i < ai; i++){
            vl = vsetvli(j, RVV_E32, RVV_M8);
            vb = vlev_float32xm8(tmp_b, vl);
            vc = vlev_float32xm8(tmp_c, vl);
            vtmp = vfmulvf_float32xm8(vb, tmp_a, vl);
            vc = vfaddvv_float32xm8(vc, vtmp, vl);
            

            tmp_b += vl;
            tmp_c += vl;
        }
        vsev_float32xm8(tmp_c, vc, vl);
        
    }
    
}

int main(int argc, char **argv){
    int dim = atoi(argv[1]);

    // spmv rvv test
    struct ELLPACK_P M;
    M._rows = dim;
    M._rowPtr = new Val*[M._rows];
    M._coordPtr = new unsigned int*[M._rows];
    M._nnz_num = new unsigned int[M._rows];
    M._nnz_old_num = new unsigned int[M._rows];

    // construct original sparse matrix
    // cout << "original matrix" << endl;
    Val *dense_matrix = new Val[M._rows*dim];
    for(int i=0; i<M._rows; i++){
        for(int j=0; j<dim; j++){
            if(rand()&127!=0){// generate 0
                dense_matrix[i*dim+j]=dis(gen);
            }
            else{
                dense_matrix[i*dim+j]=0.0;
            }   
        }
    }
    // print(dense_matrix, dim, dim);
    M.ELLPACKinitWithDenseMatrix(dense_matrix, dim, dim);
    // M.print_rowidx();
    M.ELLPACK_P_vreorder();
    // M.print_rowidx();
    // M.print_ELLPACK_P();

    Val *b = new Val[dim*dim];
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            b[i*dim+j]=dis(gen);
        }
    }
    // print(b, dim, dim);
    
    Val *output_vec = new Val [dim*dim];
    Val *output_rvv_normal_red = new Val [dim*dim];
    Val *output_rvv_normal_row = new Val [dim*dim];
    Val *output_scalar_ELL = new Val [dim*dim];
    Val *output_scalar_normal = new Val [dim*dim];

    auto start_rvv = chrono::high_resolution_clock::now();
    spmm_rvv(M, b, dim, dim, output_vec, dim, dim);
    auto stop_rvv = chrono::high_resolution_clock::now();
    auto duration_rvv = chrono::duration_cast<chrono::microseconds>(stop_rvv - start_rvv);
    cout << "rvv " << duration_rvv.count() << endl;

    auto start_rvv_normal_red = chrono::high_resolution_clock::now();
    spmm_rvv_normal_red(dense_matrix, dim, dim, b, dim, dim, output_rvv_normal_red, dim, dim);    
    auto stop_rvv_normal_red = chrono::high_resolution_clock::now();
    auto duration_rvv_normal_red = chrono::duration_cast<chrono::microseconds>(stop_rvv_normal_red - start_rvv_normal_red);
    cout << "rvv_normal_red " << duration_rvv_normal_red.count() << endl;

    auto start_rvv_normal_row = chrono::high_resolution_clock::now();
    spmm_rvv_normal_row(dense_matrix, dim, dim, b, dim, dim, output_rvv_normal_row, dim, dim);    
    auto stop_rvv_normal_row = chrono::high_resolution_clock::now();
    auto duration_rvv_normal_row = chrono::duration_cast<chrono::microseconds>(stop_rvv_normal_row - start_rvv_normal_row);
    cout << "rvv_normal_row " << duration_rvv_normal_row.count() << endl;

    auto start_scalar_ELLP = chrono::high_resolution_clock::now();
    spmm_scalar_ELLP(M, b, dim, dim, output_scalar_ELL, dim, dim);
    auto stop_scalar_ELLP = chrono::high_resolution_clock::now();
    auto duration_scalar_ELLP = chrono::duration_cast<chrono::microseconds>(stop_scalar_ELLP - start_scalar_ELLP);
    cout << "scalar_ELLP " << duration_scalar_ELLP.count() << endl;

    auto start_scalar_normal = chrono::high_resolution_clock::now();
    spmm_scalar_normal(dense_matrix, dim, dim, b, dim, dim, output_scalar_normal, dim, dim);
    auto stop_scalar_normal = chrono::high_resolution_clock::now();
    auto duration_scalar_normal = chrono::duration_cast<chrono::microseconds>(stop_scalar_normal - start_scalar_normal);
    cout << "scalar_normal " << duration_scalar_normal.count() << endl;
    cout << endl;

    delete [] M._rowPtr;
    delete [] M._coordPtr;
    delete [] M._nnz_num;
    delete [] M._nnz_old_num;

}