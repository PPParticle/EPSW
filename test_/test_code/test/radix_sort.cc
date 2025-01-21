/*     
    Author: fcqiao
    Data: 2023/6/14 19:45
    rvv radix sort algorithm
    test input: integer array dimension
    function input: *src, *dst, array dimension, logic (radix)
*/
#include <iostream>
#include <riscv_vector.h>
#include <random>
#define VLMAX 4 * 8
using namespace std;
default_random_engine generator;
uniform_int_distribution<int> distribution(0, 256);

// return merge value number
unsigned int merge(const unsigned int *src1, const unsigned int n1, const unsigned int *src2, const unsigned int n2, unsigned int *dst)
{
    size_t vl, k, num = 0;
    vl = 0;
    k = n1;
    // if(k>VLMAX){
    for (; k > 0; k -= vl)
    {
        vl = vsetvl_e32m8(k);
        vuint32m8_t v1 = vle32_v_u32m8(src1, vl);
        vse32_v_u32m8(dst, v1, vl);
        src1 += vl;
        dst += vl;
    }
    num += n1;
    // }
    vl = 0;
    k = n2;
    for (; k > 0; k -= vl)
    {
        vl = vsetvl_e32m8(k);
        vuint32m8_t v1 = vle32_v_u32m8(src2, vl);
        vse32_v_u32m8(dst, v1, vl);
        src2 += vl;
        dst += vl;
    }
    num += n2;
    return num;
}

void radix_sort(const unsigned int *src, unsigned int *dst, unsigned int *idx, const unsigned int n, unsigned int logic)
{   
    if (n > 1)
    {
        size_t vl;
        size_t n0 = 0, n1 = 0;
        unsigned int *dst_0 = new unsigned int[n];
        unsigned int *dst_1 = new unsigned int[n];
        unsigned int *idx_0 = new unsigned int[n];
        unsigned int *idx_1 = new unsigned int[n];

        unsigned int *tmp_0 = dst_0, *tmp_1 = dst_1;
        unsigned int *tmp_i_0 = idx_0, *tmp_i_1 = idx_1;

        const unsigned int *tmp_s = src;
        const unsigned int *tmp_i = idx;
        
        // cout << dst_1[1] << endl;
        for (unsigned int k = n; k > 0; k -= vl)
        {
            vl = vsetvl_e32m8(k);
            vuint32m8_t vtmp = vle32_v_u32m8(tmp_s, vl);
            vuint32m8_t vidx = vle32_v_u32m8(tmp_i, vl);

            vuint32m8_t v_and = vand_vx_u32m8(vtmp, logic, vl);
            vbool4_t v_neq = vmseq_vx_u32m8_b4(v_and, 0, vl);
            vbool4_t v_eq = vmnot_m_b4(v_neq, vl);
            size_t n_1 = vcpop_m_b4(v_eq, vl);
            size_t n_0 = vcpop_m_b4(v_neq, vl);
            // cout << n_1 << endl;

            vuint32m8_t v_1 = vcompress_vm_u32m8(v_eq, v_1, vtmp, vl);
            vuint32m8_t v_0 = vcompress_vm_u32m8(v_neq, v_0, vtmp, vl);
            vse32_v_u32m8(tmp_1, v_1, n_1);
            vse32_v_u32m8(tmp_0, v_0, n_0);

            tmp_0 += n_0;
            tmp_1 += n_1;

            v_1 = vcompress_vm_u32m8(v_eq, v_1, vidx, vl);
            v_0 = vcompress_vm_u32m8(v_neq, v_0, vidx, vl);
            vse32_v_u32m8(tmp_i_1, v_1, n_1);
            vse32_v_u32m8(tmp_i_0, v_0, n_0);

            tmp_i_0 += n_0;
            tmp_i_1 += n_1;


            n0 += n_0;
            n1 += n_1;
            
            tmp_s += vl;
            tmp_i += vl;

        }
        merge(idx_1, n1, idx_0, n0, idx);
        merge(dst_1, n1, dst_0, n0, dst);

        delete[] dst_0;
        delete[] dst_1;
        delete[] idx_0;
        delete[] idx_1;
        dst_0 = dst_1 = nullptr;

    }
}

// int main(int argc, char **argv)
// {   
//     size_t n = atoi(argv[1]);
//     unsigned int *src0 = new unsigned int[n];    
//     unsigned int *dst = new unsigned int[n];
//     for (unsigned int i = 0; i < n; i++)
//     {
//         src0[i] = distribution(generator);
//         cout << src0[i] << " ";
//     }
//     cout << endl;
//     unsigned int *src = src0;

//     unsigned int shift = 1;
//     unsigned int logic = 1;
//     for (; shift < 10; shift++, logic = logic << 1){
//         radix_sort(src, dst, n, logic);
//         src = dst;
//     }
//     // cout << "dst" << endl;
//     // for (int i = 0; i < 10; i++)
//     //     cout << dst[i] << " ";
//     for (int i = 0; i < n; i++)
//         cout << dst[i] << " ";
//     cout << endl;
//     return 0;
// }
