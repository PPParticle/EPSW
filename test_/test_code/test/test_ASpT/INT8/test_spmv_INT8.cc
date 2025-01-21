#include<iostream>
#include<cmath>
#include<random>
#include<time.h>
#include<riscv_vector.h>

#define Val float
using namespace std;

struct ELLPACK_P{
public:
    /* values and coord has identical organization */  
    Val **_rowPtr;
    Idx**_coordPtr;

    /* the number of values of each row*/
    int *_nnz_num;
    int *_nnz_old_num;

    /* fault-only-first instructions are used to vectorize loops with data-dependent */
    int _nnz, _rows;

    /* strip accomodating cache size */
    int _strips;
};

Val *spmv_rvv(ELLPACK_P M, Val *fvalues, int rows){
    Val *output = new Val[rows];
    size_t vl;
    vfloat32m8_t va, vv, vtmp;
    vuint32m8_t vidx;
    Val *ans = output;

    for(int i=0; i<M._rows; i++){

        Val          *tmp_m  = M._rowPtr[i];
        int           j      = M._nnz_num[i];
        Idx*tmp_idx= M._coordPtr[i];
        Val          *tmp_v  = fvalues;

        vtmp = vfmv_v_f_f32m8(0.0, 32);
        for(; j>0; j-=vl){
            // cout << j << endl;
            vl = vsetvl_e32m8(j);
            va = vle32_v_f32m8(tmp_m, vl);
            vidx = vle32_v_u32m8(tmp_idx, vl);
            vidx = vsll(vidx, 2, vl);
            vv = vluxei32_v_f32m8(tmp_v, vidx, vl); 
            vtmp = vfmacc_vv_f32m8(vtmp, va, vv, vl);

            tmp_m+=vl;
            tmp_v+=vl;
        }

        vfloat32m1_t vans = vfmv_v_f_f32m1(0.0, 4);
        vans = vfredusum_vs_f32m8_f32m1(vans, vtmp, vans, 32);
        vse32_v_f32m1(ans, vans, 1);
        ans++;
    }
    return output;
}

Val *spmv_scalar(ELLPACK_P M, Val *fvalues, int rows){
    Val *output = new Val[rows];
    for(int i=0; i<M._rows; i++){
        for(int j=0; j<rows; j++){
            output[i]+=M._rowPtr[i][j]*fvalues[j];
        }
    }
    return output;
}

int main(){

    // spmv rvv test
    struct ELLPACK_P M;
    M._rows = 64;
    M._rowPtr = new Val*[M._rows];
    M._coordPtr = new unsigned int*[M._rows];
    M._nnz_num = new int[M._rows];
    for(int i=0; i<M._rows; i++)
    M._nnz_num[i]=64;

    // default_random_engine generator;
    // uniform_real_distribution<float> distrib(0, 256);
    for(int i=0; i<M._rows; i++){
        Val *tmp_row = new Val[64];
        Val val = 0.0;
        Idx*tmp_idx = new unsigned int[64];
        for(int j=0; j<64; j++){
            // tmp_row[j]=(double)rand() / (double)RAND_MAX + (double)(rand() % 10);
            tmp_row[j] = val + j;
            tmp_idx[j]=j;
        }
        M._rowPtr[i]=tmp_row;
        M._coordPtr[i]=tmp_idx;
    }

    Val *vec = new Val[64];
    Val val = 0.0;
    for(int i=0; i<64; i++){
        // vec[i]=(double)rand() / (double)RAND_MAX + (double)(rand() % 10);
        vec[i] = val + i;
    }

    cout << "matrix and vector initiation complete" << endl;

    Val *vec_scalar = spmv_scalar(M, vec, 64);
    Val *vec_rvv = spmv_rvv(M, vec, 64);
    
    for(int i=0; i<64; i++){
        if(vec_rvv[i]==vec_scalar[i]) ;
        else{
            cout << "false" << endl;
            cout << vec_rvv[i] << " " << vec_scalar[i] << endl;
            return -1;
        }
    }

    for(int i=0; i<64; i++){
        if(vec_rvv[i]==vec_scalar[i]) ;
        else{
            cout << "false" << endl;
            cout << vec_rvv[i] << " " << vec_scalar[i] << endl;
            return -1;
        }
    }
    cout << "success" << endl;

    for(int i=0; i<64; i++){
        cout << vec_rvv[i] << " ";
    }
    cout << endl;
    return 0;
}