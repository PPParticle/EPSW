#include <iostream>
#include <riscv_vector.h>
#include <assert.h>
#include <random>
#include <chrono>
#include "utilis.h"
#include "ELLPACK_P.h"
using namespace std;
default_random_engine generator;
normal_distribution<float> distribution(0.0, 4.0);

void print_ELLPACK_P(ELLPACK_P M)
{
    cout << "ELLPACK_P M: " << endl;
    for (Idxi = 0; i < M._rows; i++)
    {
        if (M._nnz_num[i] == 0)
            continue;
        for (Idxj = 0; j < M._nnz_num[i]; j++)
            cout << M._coordPtr[i][j] << " ";
        cout << "         ";
        if (M._rowPtr[i] != nullptr)
            for (Idxj = 0; j < M._nnz_num[i]; j++)
                cout << M._rowPtr[i][j] << " ";
        cout << endl;
    }
}

void sddmm_ELLP(Idx*c, const int ci, const int cj,
                float *a, const int ai, const int aj,
                float *b, const int bi, const int bj,
                float *o, float *o_s)
{
    ELLPACK_P M;
    M._rows = ci;
    M._rowPtr = new Val *[M._rows];
    M._coordPtr = new Idx*[M._rows];
    M._nnz_old_num = new unsigned int[M._rows];
    M._nnz_num = new unsigned int[M._rows];

    size_t vl;
    vfloat32m8_t va, v_a;
    vuint32m8_t vi, v_i;
    vuint32m8_t vidx;
    vbool4_t vmask;

    // for (int i = 0; i < ci; i++)
    // {
    //     int nnz = 0;
    //     for (int j = cj; j > 0; j--)
    //     {
    //         if (c[i * cj + j] != 0)
    //             nnz++;
    //     }
    //     if (nnz == 0)
    //         continue;

    //     M._nnz_old_num[i] = nnz;
    //     int num = 0;
    //     Idx*tmp_c = c + i * cj;
    //     unsigned *tmp_i = new unsigned int[nnz];
    //     M._coordPtr[i] = tmp_i;

    //     for (int j = cj; j > 0; j -= vl)
    //     {
    //         vl = vsetvl_e32m8(j);
    //         vidx = vle32_v_u32m8(tmp_c, vl);
    //         vmask = vmsne(vidx, 0, vl);
    //         v_i = vcompress_vm_u32m8(vmask, v_i, vidx, vl);
    //         num = vcpop_m_b4(vmask, vl);
    //         vse32_v_u32m8(tmp_i, v_i, num);

    //         tmp_i += num;
    //         tmp_c += vl;
    //     }
    // }
    // M._nnz_num = M._nnz_old_num;
    // print_ELLPACK_P(M);

    // two reasons for not transform sampled matrix into ELLPACK_P:
    // 1. the index span might out of maxvl because of extreme sparsity
    // 2. sampled matrix is knew on runtime, requiring timely transformation which has overhead
    // print(c, ci, cj);

    Idxm = ai;
    Idxp = bi;
    Idxn = bj;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if(c[i*n+j]==0) continue;
            for (int k = 0; k < p; k++)
            {
                o[i*n + j] += a[i*p + k] * b[k*n + j];
            }
        }
    }

    // print(o, n, n);

    for (int i = 0; i < m; i++)
    {   
        vfloat32m8_t vans=vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
        for (int j = 0; j < p; j++)
        {   
            Val *tmp_a = a + i * p + j;  // a[i][j]
            Val *tmp_o_s = o_s + i * n;  // o_s[i]
            Val *tmp_b = b + j * n;      // b[j]
            Idx*tmp_c = c + i * n;
            vuint32m8_t vc;
            vfloat32m8_t vb, vtmp;
            vbool4_t vmask;
            for (int k = n; k > 0; k -= vl)
            {   
                vl = vsetvl_e32m8(k); 
                vc = vle32_v_u32m8(tmp_c, vl);
                vb = vle32_v_f32m8(tmp_b, vl);
                vans = vle32_v_f32m8(tmp_o_s, vl);

                vmask = vmsne(vc, 0, vl);

                vtmp = vfmul_vf_f32m8_m(vmask, vtmp, vb, *tmp_a, vl);
                vans = vfadd_vv_f32m8_m(vmask, vans, vans, vtmp, vl);
                vse32_v_f32m8_m(vmask, tmp_o_s, vans, vl);
                
                tmp_c += vl;
                tmp_o_s += vl;
                tmp_b += vl;
            }

        }
    }
}

void initiation(float *a, Idxai, Idxaj,
                float *b, Idxbi, Idxbj,
                Idx*c, Idxci, Idxcj)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < aj; j++)
        {
            a[i * aj + j] = distribution(generator);
            b[i * bj + j] = distribution(generator);
            float tmp = distribution(generator);

            if (rand() & ci)
                c[i * cj + j] = 0;
            else
                c[i * cj + j] = j;
        }
    }
}

void matrixtoELLPACK(int *c, Idxci, Idxcj, int **c_new, int *size, int **c_d)
{

    for (int i = 0; i < ci; i++)
    {
        int nnz = 0;
        for (int j = 0; j < cj; j++)
        {
            if (c[i * cj + j])
            {
                nnz++;
            }
        }
        int *tmp = new int[nnz];   // tmp row space
        int *tmp_d = new int[nnz]; // tmp index space
        size[i] = nnz;
        c_new[i] = tmp;
        c_d[i] = tmp_d;
    }

    for (int i = 0; i < ci; i++)
    {
        for (int j = 0; j < cj; j++)
        {
            if (c[i * cj + j])
            {
                c_new[i][j] = 1;
                c_d[i][j] = j;
            }
        }
    }
}

void col_major(float **a, float *b, Idxbi, Idxbj)
{
    for (int i = 0; i < bi; i++)
    {
        for (int j = 0; j < bj; j++)
        {
            if (a[j][i] != 0.0)
                b[i * bj + j] = a[j][i];
            else
                b[i * bj + j] = 0.0;
        }
    }
}

int main(int argc, char *argv[])
{
    // if(argc<2) cout << "please type the d of square-matirx" << endl;
    int n = 0;
    cin >> n;
    // sddmm output = C * A * B
    float *a = new float[n * n];
    float *b = new float[n * n];
    Idx*c = new unsigned int[n * n];

    float *o = new float[n * n];
    float *o_s = new float[n * n];

    initiation(a, n, n, b, n, n, c, n, n);
    // matrixtoELLPACK(c,n,n,c_new,size,c_d);

    // cout << "matrix a" << endl;
    // print(a, n, n);
    // cout << "matrix b" << endl;
    // print(b, n, n);
    // cout << "matrix c" << endl;
    // print(c, n, n);

    // col_major(a_new, a_col_major, n, n);
    // col_major(a_d, n, size, a_col_d, n, n);

    // cout << "col_major_value" << endl;
    // print(a_col_major, n, n);
    // cout << "col_major index" << endl;
    // print(a_col_d, n ,n);

    auto start_in_row = chrono::high_resolution_clock::now();
    sddmm_ELLP(c, n, n, a, n, n, b, n, n, o, o_s);
    auto stop_in_row = chrono::high_resolution_clock::now();
    auto duration_in_row = chrono::duration_cast<chrono::microseconds>(stop_in_row - start_in_row);

    // cout << "original" << endl;
    // print(o, n, n);
    // cout << "sampled" << endl;
    // print(o_s, n, n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(o_s[i*n+j]==o[i*n+j]);
            else{ cout << i << j << " " << o_s[i*j+j] << " " << o[i*j+j] << endl;
            return -1;}
        }
    }

    return 0;
}