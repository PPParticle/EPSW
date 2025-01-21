#include<iostream>
#include<riscv_vector.h>
#include"ELLPACK_P_INT8.h"
#include<random>
#include"utilis_INT8.h"

#define dim 5
default_random_engine gen;
uniform_int_distribution<int> dis(-10,10);
using namespace std;

void radix_sort(const Idx *src, Idx *dst, Idx *idx, const Idx n, Idx logic);

// inner dataflow scalar kernel
void spmm_scalar_normal(Val *a, const Idx ai, const Idx aj,
                        Val *b, const Idx bi, const Idx bj, 
                        Val *c, const Idx ci, const Idx cj){
    Idx m=ai;
    Idx p=aj;
    Idx n=bj;
    
    for(int i=0; i<ci; i++){
        for(int j=0; j<cj; j++){
            for(int k=0; k<p; k++){
                // if(i==0&&j==0){
                //     printf("%5.2d %5.2d %5.2d \n", (int32_t)c[i*n+j], (int32_t)a[i*p+k],
                //     (int32_t)b[k*n+j]);
                // }
                c[i*n+j]+=a[i*p+k]*b[k*n+j];
            }
        }
    } 
}

// row_based dataflow scalar kernel
// SpM m*p  M p*n 
void spmm_scalar_ELLP(ELLPACK_P M, Val *fvalues, const Idx rows, const Idx cols,
                        Val *c, const Idx ci, const Idx cj){
    Idx m=M._rows;
    Idx p=rows;
    Idx n=cols;

    for(int i=0; i<M._rows; i++){
        for(int j=0; j<M._nnz_num[i]; j++){
            Val brd = M._rowPtr[i][j];
            Idx brd_idx = M._coordPtr[i][j];
            Idx row_idx = M._rowidx[i];
            for(int k=0; k<n; k++){
                c[row_idx*n+k] += brd * fvalues[brd_idx*n+k];
            }
        }
    }
}

void spmm_rvv(ELLPACK_P M, Val *fvalues, const Idx rows, const Idx cols,
                Val *c, const Idx ci, const Idx cj){
    Idx m=M._rows;
    // for(int i=0 ; i<M._rows; i++){
    //     cout << M._nnz_num[i] << " ";
    // }
    // cout << endl;
    Idx p=rows;
    Idx n=cols;
    size_t vl;
    for(int i=0; i<M._rows; i++){
        for(int j=0; j<M._nnz_num[i]; j++){
            Val          brd     = M._rowPtr[i][j];
            Idx          brd_idx = M._coordPtr[i][j];
            Val          *tmp_b  = fvalues+brd_idx*n;
            Val          *tmp_c  = c+M._rowidx[i]*n;
            vint8m8_t vans = vmv_v_x_i8m8(0, vsetvlmax_e8m8());
            for(int k=n; k>0; k-=vl){
                vl=vsetvl_e8m8(k);
                vint8m8_t vb=vle8_v_i8m8(tmp_b, vl);
                vint8m8_t vtmp=vmul_vx_i8m8(vb, brd, vl);
                vb=vle8_v_i8m8(tmp_c, vl);
                vans = vadd_vv_i8m8(vtmp, vb, vl);
                vse8_v_i8m8(tmp_c, vans, vl);

                tmp_c += vl;
                tmp_b += vl;
            }
        }
    }
}

int main(){

    // spmv rvv test
    struct ELLPACK_P M;
    M._rows = dim;
    M._rowPtr = new Val*[M._rows];
    M._coordPtr = new Idx*[M._rows];
    M._nnz_num = new Idx[M._rows];
    M._nnz_old_num = new Idx[M._rows];

    // construct original sparse matrix
    cout << "original matrix" << endl;
    Val *dense_matrix = new Val[M._rows*dim];
    for(int i=0; i<M._rows; i++){
        for(int j=0; j<dim; j++){
            if(rand()&256!=0){// generate 0
                dense_matrix[i*dim+j]=dis(gen);
            }
            else{
                dense_matrix[i*dim+j]=0.0;
            }   
        }
    }
    print(dense_matrix, dim, dim);
    M.ELLPACKinitWithDenseMatrix(dense_matrix, dim, dim);
    // M.print_rowidx();

    M.ELLPACK_P_vreorder();
    // M.print_rowidx();

    M.print_ELLPACK_P();

    Val *b = new Val[dim*(dim-1)];
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            b[i*dim+j]=dis(gen);
        }
    }
    print(b, dim, dim-1);
    
    Val *output_vec = new Val[dim*(dim-1)]{0};
    Val *output_scalar_ELL = new Val[dim*(dim-1)]{0};
    Val *output_scalar_normal = new Val[dim*(dim-1)]{0};
    spmm_rvv(M, b, dim, (dim-1), output_vec, dim, (dim-1));
    spmm_scalar_ELLP(M, b, dim, (dim-1), output_scalar_ELL, dim, (dim-1));
    spmm_scalar_normal(dense_matrix, dim, dim, b, dim, (dim-1), output_scalar_normal, dim, (dim-1));
    
    for(int i=0; i<dim; i++){
        for(int j=0; j<(dim-1); j++){
            // if(output_vec[i*dim+j]==output_scalar_ELL[i*dim+j]&&output_vec[i*dim+j]==output_scalar_normal[i*dim+j]){
            //     ;
            // }
            // else{
            //     cout << "false" << endl;
            //     cout << i << " " << j << " " ;
                printf("%d %d %d     ",(int32_t)output_vec[i*dim+j], (int32_t)output_scalar_ELL[i*dim+j], (int32_t)output_scalar_normal[i*dim+j]); 
        //         return -1;
        //     }

        }
        cout << endl;
    }
    cout << "success" << endl;

    delete [] M._rowPtr;
    delete [] M._coordPtr;
    delete [] M._nnz_num;
    delete [] M._nnz_old_num;

}