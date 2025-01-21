#include<iostream>
#include<riscv_vector.h>
#include<assert.h>
#include<random>
#include<chrono>
#include"utilis.h"
using namespace std;
default_random_engine generator;
normal_distribution<float> distribution(0.0, 1.0);

// 预编译所执行的操作就是简单的 “文本” 替换
#define VSETVL vsetvl_e32m8
#define VSETVLMAX vsetvlmax_e32m8()
#define VFLOAT vfloat32m8_t
#define VLE_FLOAT vle32_v_f32m8
#define VLSE_FLOAT vlse32_v_f32m8
#define VSE_FLOAT vse32_v_f32m8
#define VFMACC_VV vfmacc_vv_f32m8
#define VFREDSUM vfredusum_vs_f32m8_f32m1

float *in_row_parallelism(float **a_new, unsigned int ai, unsigned int **a_d, int *size,
                            float *b, unsigned int bi,
                            float *c, unsigned int ci){
    size_t vl;
    VFLOAT va, vb, vtmp;
    vuint32m8_t vd;
    float* ans=c;
    for(int i=0; i<ai; i++){
        float* tmp_a=a_new[i];
        unsigned int* tmp_d=a_d[i];
        int s=size[i];
        float* tmp_b=b;
        vtmp=vfmv_v_f_f32m8(0.0, VSETVLMAX);
        for(int k=s; k>0; k-=vl){
            vl=VSETVL(k);
            va=VLE_FLOAT(tmp_a, vl);
            vd=vle32_v_u32m8(tmp_d, vl);
            vd=vsll(vd, 2, vl);
            vb=vluxei32_v_f32m8(tmp_b, vd, vl); 
            vtmp=VFMACC_VV(vtmp, va, vb, vl);
            tmp_a+=vl;
            tmp_b+=vl;
        } 
        
    // vfloat32m1_t vans=vfmv_v_f_f32m1(0.0, 1);
    vfloat32m1_t vans;
    vans=vfredusum_vs_f32m8_f32m1(vans, vtmp, vans, 32);
    vse32_v_f32m1(ans, vans, 1);
    ans++;
    }   
    return c;
}

float *cross_line_parallelism(float *a, unsigned int ai, unsigned int aj,
                            unsigned int *a_d, unsigned int adi, unsigned int adj,
                            float *b, unsigned int bi,
                            float *c, unsigned int ci){
    size_t vl;
    VFLOAT va, vb, vtmp, vc;
    vuint32m8_t vidx;
    vbool4_t mask;
    float* tmp_b=b;
    for(int i=0; i<ai; i++){
        float *tmp_a=a+i*aj;
        unsigned int *tmp_i=a_d+i*aj;
        float *tmp_b=b;
        float *tmp_c=c;
        vtmp=vfmv_v_f_f32m8(0.0, VSETVLMAX);
        for(int i=0; i<ai; i+=vl){
            vl=VSETVL(ai-i);
            va=VLE_FLOAT(tmp_a, vl);
            vidx=vle32_v_u32m8(tmp_i, vl);
            mask=vmsne_vx_u32m8_b4(vidx, 0, vl);
            vidx=vsub_vx_u32m8(vidx, 1, vl);
            vidx=vsll_vx_u32m8(vidx, 2, vl);
            vb=vluxei32_v_f32m8_m(mask, vb, tmp_b, vidx, vl);
            vtmp=vfmul_vv_f32m8_m(mask, vtmp, va, vb, vl);

            vc=vluxei32_v_f32m8_m(mask, vc, tmp_c, vidx, vl);
            vtmp=vfadd_vv_f32m8_m(mask, vtmp, vtmp, vc, vl);
            vsuxei32_v_f32m8_m(mask, tmp_c, vidx, vtmp, vl);
            tmp_a+=vl;
            tmp_i+=vl;
            tmp_b+=vl;
            tmp_c+=vl;
        }
    }
    return c;
}

void initiation(float *a, unsigned int ai, unsigned int aj, float *v, unsigned int vi){
    for(int i=0; i<ai; i++){
        for(int j=0; j<aj; j++){
            a[i*aj+j]=distribution(generator);
        }
    }

    float val=0.0;
    for(int i=0; i<vi; i++){
        v[i]=distribution(generator);
    }

}
void zero(float *a, unsigned int ai){
    for(int i=0; i<ai; i++){
        if(a[i]>-1.0&&a[i]<1.0) a[i]=0.0;
    }
}
void zero(float *a, unsigned int ai, unsigned int aj){
    for(int i=0; i<ai; i++){
        for(int j=0; j<aj; j++){
            if(a[i*aj+j]>-1.0&&a[i*aj+j]<1.0) a[i*aj+j]=0.0;
        }
    }
}
void permute(float *a, unsigned int ai, unsigned int aj, float **a_new, unsigned int **a_d, int *size){

    for(int i=0; i<ai; i++){
        int nnz=0;
        for(int j=0; j<aj; j++){
            if(a[i*aj+j]!=0.0){
                nnz++;
            }
        }
        float *tmp=new float[nnz];
        unsigned int *tmp_d=new unsigned int[nnz];
        size[i]=nnz;
        a_new[i]=tmp;
        a_d[i]=tmp_d;
    }

    for(int i=0; i<ai; i++){
        int tmp=0;
        for(int j=0; j<aj; j++){
            if(a[i*aj+j]!=0.0){
                a_new[i][tmp]=a[i*aj+j];
                a_d[i][tmp]=j;
                tmp++;
            }
        }
    }
}
void col_major(float **a, float *b, unsigned int bi, unsigned int bj){
    for(int i=0; i<bi; i++){
        for(int j=0; j<bj; j++){
            if(a[j][i]!=0.0)
            b[i*bj+j]=a[j][i];
            else b[i*bj+j]=0.0;
        }
    }
}
void col_major(unsigned int **a_d, unsigned int a_di, int *size, unsigned int *b, unsigned int bi, unsigned int bj){
    for(int j=0; j<bj; j++){
        for(int i=0; i<bi; i++){
            if(i<size[j])
            b[i*bj+j]=a_d[j][i]+1;
            else b[i*bj+j]=0;
        }
    }

}

int main(int argc, char*argv[]){
    // if(argc<2) cout << "please type the d of square-matirx" << endl;
    int n=0;
    cin >> n;
    float *a=new float[n*n];
    float *a_col_major=new float[n*n];
    float **a_new=new float*[n];
    unsigned int **a_d=new unsigned int*[n];
    unsigned int *a_col_d=new unsigned int[n*n];
    int *size=new int[n];
    float *b=new float[n];
    float *c_in_row=new float[n];
    float *c_cross_line=new float[n];

    initiation(a,n,n,b,n);
    zero(a,n,n);
    zero(b,n);

    // cout << "matrix" << endl;
    // print(a,n,n);
    // cout << "vector" << endl;
    // print(b,n);

    permute(a,n,n,a_new,a_d,size);

    // cout << "permute matrix" << endl;
    // print(a_new,n,size);
    // cout << "permute vector" << endl;
    // print(a_d,n,size);
    // cout << "permute nnz" << endl;
    // for(int i=0; i<n; i++){
    //     cout <<size[i] << " ";
    // }
    // cout << endl;


    // col_major(a_new, a_col_major, n, n);
    // col_major(a_d, n, size, a_col_d, n, n);

    // cout << "col_major_value" << endl;
    // print(a_col_major, n, n);
    // cout << "col_major index" << endl;
    // print(a_col_d, n ,n);
    auto start_in_row=chrono::high_resolution_clock::now();
    c_in_row=in_row_parallelism(a_new,n,a_d,size,b,n,c_in_row,n);
    auto stop_in_row=chrono::high_resolution_clock::now();
    auto duration_in_row=chrono::duration_cast<chrono::microseconds>(stop_in_row-start_in_row);

    auto start_cross_line=chrono::high_resolution_clock::now();
    c_cross_line=cross_line_parallelism(a_col_major,n,n,a_col_d,n,n,b,n,c_in_row,n);
    auto stop_cross_line=chrono::high_resolution_clock::now();
    auto duration_cross_line=chrono::duration_cast<chrono::microseconds>(stop_cross_line-start_cross_line);

    // print(c_in_row, n);
    // print(c_cross_line, n);

    for(int i=0; i<n; i++){
        if(c_in_row[i]!=c_cross_line[i]){
            cout << "error" << endl;
            cout << i << " " << c_in_row[i] << " " << c_cross_line[i] << endl;
            return -1;
        }
    }

    float *distribution=new float[n];
    for(int i=0; i<n; i++){
        float tmp=0.0;
        for(int j=0; j<size[i]; j++){
            if(a[i*n+j]!=0.0) tmp+=1.0;
        }
        distribution[i]=tmp/float(size[i]);
    }
    print(distribution, n);
    cout << "in_row" << " " << duration_in_row.count() << endl;
    cout << "cross_line" << " " << duration_cross_line.count() << endl;

    return 0;
}