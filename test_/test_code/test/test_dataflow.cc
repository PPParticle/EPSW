#include<iostream>
#include<riscv_vector.h>
#include<time.h>
#include<sys/time.h>
#include<random>
using namespace std;
#define LEN_1D 512
#define LEN_2D 224
#define iterations 10

float32_t a[224][224]={};
float32_t b[224][224]={};
float32_t c[224][224]={};

struct args_t {
    struct timeval t1;
    struct timeval t2;
    void * __restrict__ arg_info;
};

typedef float(*test_function_t)(struct args_t *);

void time_function(test_function_t vector_func, void * arg_info)
{
    struct args_t func_args = {.arg_info=arg_info};

    double result = vector_func(&func_args);

    double tic=func_args.t1.tv_sec+(func_args.t1.tv_usec/1000000.0);
    double toc=func_args.t2.tv_sec+(func_args.t2.tv_usec/1000000.0);

    double taken = toc-tic;

    printf("%10.3f\t%f\n", taken, result);
}

void init_arr(){
    default_random_engine generator;
    uniform_real_distribution<float32_t> distrib(.0, 256.0 );
    for(int i=0; i<224; i++){
        for(int j=0; j<224; j++){
            a[i][j] = distrib(generator);
            b[i][j] = distrib(generator);
            c[i][j] = 0;
        }
    }
}

float sum(float arr[][224], int r){
    float ret = 0;
    for (int i = 0; i < r; i++)
        for(int j=0; j < LEN_1D; j++)
            ret += arr[i][j]; 
    return ret;
}

float matrix_inner_m1(struct args_t * func_args)
{
    init_arr();
    gettimeofday(&func_args->t1, NULL);
    for (int nl = 0; nl < iterations; nl++) {
        size_t vl;
        float* a_ptr = &a[0][0];
        float* b_ptr = &b[0][0];
        float *c_ptr = &c[0][0];
        for(int j=0; j<224; j++){
            for (int i = 0; i < 224; i=i+vl) {
                vl = vsetvl_e32m8(224-i);
                vfloat32m8_t va =  vle32_v_f32m8(a_ptr, vl);
                vfloat32m8_t vb =  vlse32_v_f32m8(b_ptr, 224*4, vl);
            
            
                va = vfmul_vv_f32m8(va, vb, vl);
                vfloat32m1_t vans = vfmv_s_f_f32m1(vans, 0.0, 1);
                vans = vfredusum_vs_f32m8_f32m1(vans, va, vans, vl);
                float32_t ans_c = vfmv_f_s_f32m1_f32(vans);
                *c_ptr += ans_c;
            }
            c_ptr++;
        }
    }

    gettimeofday(&func_args->t2, NULL);
    return sum(a, LEN_1D);
}

// float matrix_outer_m1(struct args_t * func_args)
// {
//     init_arr();
//     gettimeofday(&func_args->t1, NULL);
//     for (int nl = 0; nl < iterations; nl++) {
//         size_t vl;
//         float* a_ptr = &a[0][0];
//         float* b_ptr = &b[0][0];
//         float* c_ptr = &c[0][0];
//         for (int i = 0; i < LEN_1D; i=i+vl+2) {
//             vl = vsetvl_e32m8(LEN_1D-i);
//             vfloat32m8_t vb =  vle32_v_f32m8(b_ptr, vl);
//             vfloat32m8_t va =  vle32_v_f32m8(a_ptr, vl);
//             va = vfmul_vv_f32m8(va, vb, vl);
//             // if(LEN_1D-vl > 2){
//             a[i+vl] = b[i+vl] * a[i+vl];
//             a[i+vl+1] = b[i+vl+1] * a[i+vl+1];
//             a_ptr = a_ptr + vl + 2;
//             b_ptr = b_ptr + vl + 2;
//             // }
//             vse32_v_f32m8(a, va, vl);
//         }
//         a_ptr = &a[0];
//     }

//     gettimeofday(&func_args->t2, NULL);
//     return sum(a, LEN_1D);
// }

// float matrix_row_based_m1(struct args_t * func_args)
// {
//     init_arr();
//     gettimeofday(&func_args->t1, NULL);
//     for (int nl = 0; nl < iterations; nl++) {
//         size_t vl;
//         float* a_ptr = &a[0][0];
//         float* b_ptr = &b[0][0];
//         for (int i = 0; i < LEN_1D; i=i+vl+2) {
//             vl = vsetvl_e32m8(LEN_1D-i);
//             vfloat32m8_t vb =  vle32_v_f32m8(b_ptr, vl);
//             vfloat32m8_t va =  vle32_v_f32m8(a_ptr, vl);
//             va = vfmul_vv_f32m8(va, vb, vl);
//             // if(LEN_1D-vl > 2){
//             a[i+vl] = b[i+vl] * a[i+vl];
//             a[i+vl+1] = b[i+vl+1] * a[i+vl+1];
//             a_ptr = a_ptr + vl + 2;
//             b_ptr = b_ptr + vl + 2;
//             // }
//             vse32_v_f32m8(a, va, vl);
//         }
//         a_ptr = &a[0];
//         // dummy((float*)a);
//     }

//     gettimeofday(&func_args->t2, NULL);
//     return sum(a, LEN_1D);
// }


int main(){
    args_t my_args;
    // float row = matrix_row_based_m1(&my_args);
    // float out = matrix_outer_m1(&my_args);
    float in  = matrix_inner_m1(&my_args);
    // cout << row << " " << out << " " << in << endl;
    cout << in <<  endl;


}