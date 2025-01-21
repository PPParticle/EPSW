#include<iostream>
#include<riscv_vector.h>
#include<vector>
#include"ELLPACK_P.h"
#include<random>
#include"utilis.h"

#define dim 5
default_random_engine gen;
uniform_real_distribution<double> dis(-10,10);
using namespace std;

void radix_sort(const unsigned int *src, unsigned int *dst, unsigned int *idx, const unsigned int n, unsigned int logic);

// inner dataflow scalar kernel
void spmm_scalar_normal(Val *a, const unsigned int ai, const unsigned int aj,
                        Val *b, const unsigned int bi, const unsigned int bj, 
                        Val *c, const unsigned int ci, const unsigned int cj){
    unsigned int m=ai;
    unsigned int p=aj;
    unsigned int n=bj;
    
    for(int i=0; i<ci; i++){
        for(int j=0; j<cj; j++){
            for(int k=0; k<p; k++)
            c[i*n+j]+=a[i*p+k]*b[k*n+j];
        }
    } 
}

// row_based dataflow scalar kernel
// SpM m*p  M p*n 
void spmm_scalar_ELLP(ELLPACK_P M, Val *fvalues, const int rows, const int cols,
                        Val *c, const unsigned int ci, const unsigned int cj){
    unsigned int m=M._rows;
    unsigned int p=rows;
    unsigned int n=cols;

    for(int i=0; i<M._rows; i++){
        for(int j=0; j<M._nnz_num[i]; j++){
            Val brd = M._rowPtr[i][j];
            unsigned int brd_idx = M._coordPtr[i][j];
            unsigned int row_idx = M._rowidx[i];
            for(int k=0; k<n; k++){
                c[row_idx*n+k] += brd * fvalues[brd_idx*n+k];
            }
        }
    }
}

void spmm_rvv(ELLPACK_P M, Val *fvalues, const int rows, const int cols,
                Val *c, const unsigned int ci, const unsigned int cj){
    unsigned int m=M._rows;
    // for(int i=0 ; i<M._rows; i++){
    //     cout << M._nnz_num[i] << " ";
    // }
    // cout << endl;
    unsigned int p=rows;
    unsigned int n=cols;
    size_t vl;
    for(int i=0; i<M._rows; i++){
        for(int j=0; j<M._nnz_num[i]; j++){
            Val          brd     = M._rowPtr[i][j];
            unsigned int brd_idx = M._coordPtr[i][j];
            Val          *tmp_b  = fvalues+brd_idx*n;
            Val          *tmp_c  = c+M._rowidx[i]*n;
            
            for(int k=n; k>0; k-=vl){
                vl=vsetvl_e32m8(k);
                vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
                vfloat32m8_t vtmp=vfmul_vf_f32m8(vb, brd, vl);
                vb=vle32_v_f32m8(tmp_c, vl);
                vfloat32m8_t vans = vfadd_vv_f32m8(vtmp, vb, vl);
                // vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
                vse32_v_f32m8(tmp_c, vans, vl);

                tmp_c += vl;
                tmp_b += vl;
            }
        }
    }
}

void spmm_rvv_strips(ELLPACK_P M, Val *fvalues, const int rows, const int cols,
                Val *c, const unsigned int ci, const unsigned int cj){
    unsigned int m=M._rows;
    unsigned int p=rows;
    unsigned int n=cols;
    unsigned int strips = M._strips;
    size_t vl;
    for(int i=0; i<M._rows; i+=strips){
        // cout << i << endl;
        for(int j=0; j<M._nnz_num[i]; j++){
            for(int k=0; k<strips && i+k<m; k++){
                if(M._rowPtr[i+k][j]!=NULL){
                    // cout << "i+k" << i+k << endl;
                    Val          brd     = M._rowPtr[i+k][j];
                    unsigned int brd_idx = M._coordPtr[i+k][j];
                    Val          *tmp_b  = fvalues+brd_idx*n;
                    Val          *tmp_c  = c + M._rowidx[i+k]*n;
                    // cout << "brd: " << brd << " brd_idx: "<<brd_idx<<" *tmp_b: "<<*tmp_b <<" *tmp_c: " << *tmp_c << endl;;
                    for(int l=n; l>0; l-=vl){
                        vl=vsetvl_e32m8(l);
                        vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
                        vfloat32m8_t vtmp=vfmul_vf_f32m8(vb, brd, vl);
                        vb=vle32_v_f32m8(tmp_c, vl);
                        vfloat32m8_t vans = vfadd_vv_f32m8(vtmp, vb, vl);
                        // vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
                        vse32_v_f32m8(tmp_c, vans, vl);

                        tmp_c += vl;
                        tmp_b += vl;
                    }
                }
            }          
        }
    }
}


// void spmm_rvv_strips(ELLPACK_P M, Val *fvalues, const int rows, const int cols,
//                 Val *c, const unsigned int ci, const unsigned int cj){
//     unsigned int m=M._rows;
//     unsigned int p=rows;
//     unsigned int n=cols;
//     size_t vl;
//     for(int i=0; i<M._strips.size(); i++){
//         for(int j=0; j<M._strips[i].size(); j++){
//             for(int k=0; )
//             Val          brd     = M._rowPtr[i][j];
//             unsigned int brd_idx = M._coordPtr[i][j];
//             Val          *tmp_b  = fvalues+brd_idx*n;
//             Val          *tmp_c  = c+M._rowidx[i]*n;
            
//             for(int k=n; k>0; k-=vl){
//                 vl=vsetvl_e32m8(k);
//                 vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
//                 vfloat32m8_t vtmp=vfmul_vf_f32m8(vb, brd, vl);
//                 vb=vle32_v_f32m8(tmp_c, vl);
//                 vfloat32m8_t vans = vfadd_vv_f32m8(vtmp, vb, vl);
//                 // vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
//                 vse32_v_f32m8(tmp_c, vans, vl);

//                 tmp_c += vl;
//                 tmp_b += vl;
//             }
//         }
//     }
// }

int main(){

    // spmv rvv test
    struct ELLPACK_P M;
    M._rows = dim;
    M._rowPtr = new Val*[M._rows];
    M._coordPtr = new unsigned int*[M._rows];
    M._nnz_num = new unsigned int[M._rows];
    M._nnz_old_num = new unsigned int[M._rows];

    // construct original sparse matrix
    cout << "original matrix" << endl;
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
    print(dense_matrix, dim, dim);
    M.ELLPACKinitWithDenseMatrix(dense_matrix, dim, dim);
    M.ELLPACK_P_vreorder();
    M.print_ELLPACK_P();
    M.stripe_static(3);

    cout << "mattrix b" << endl;
    Val *b = new Val[dim*dim];
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            b[i*dim+j]=dis(gen);
        }
    }
    print(b, dim, dim);
    
    Val *output_vec = new Val [dim*dim];
    Val *output_vec_strips = new Val [dim*dim];
    Val *output_scalar_ELL = new Val [dim*dim];
    Val *output_scalar_normal = new Val [dim*dim];

    spmm_rvv(M, b, dim, dim, output_vec, dim, dim);
    spmm_rvv_strips(M, b, dim, dim, output_vec_strips, dim, dim);
    spmm_scalar_ELLP(M, b, dim, dim, output_scalar_ELL, dim, dim);
    spmm_scalar_normal(dense_matrix, dim, dim, b, dim, dim, output_scalar_normal, dim, dim);
    
    // cout << endl;
    // print(output_vec, dim, dim);
    for(int i=0; i<dim; i++){
        for(int j=0; j<dim; j++){
            // printf("%6.2lf %6.2lf %6.2lf      ",output_vec[i*dim+j], output_scalar_ELL[i*dim+j], output_scalar_normal[i*dim+j]); 
            if(output_vec[i*dim+j]==output_vec_strips[i*dim+j]&&output_vec[i*dim+j]==output_scalar_ELL[i*dim+j]&&output_vec[i*dim+j]==output_scalar_normal[i*dim+j]){
                ;
            }
            else{
                cout << "false" << endl;
                cout << i << " " << j << " " ;
                printf("%6.2lf %6.2lf %6.2lf %6.2lf \n",output_vec_strips[i*dim+j], output_vec[i*dim+j], output_scalar_ELL[i*dim+j], output_scalar_normal[i*dim+j]); 
                return -1;
            }

        }
        // cout << endl;
    }
    cout << "success" << endl;

    delete [] M._rowPtr;
    delete [] M._coordPtr;
    delete [] M._nnz_num;
    delete [] M._nnz_old_num;

}