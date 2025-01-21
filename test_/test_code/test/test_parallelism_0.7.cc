#include<iostream>
#include<riscv-vector.h>
#include<assert.h>
#include<random>
#include<chrono>
using namespace std;

default_random_engine generator;

// 预编译所执行的操作就是简单的 “文本” 替换
#define VSETVL vsetvli
#define VSETVLMAX vsetvli_max(RVV_E32, RVV_M8)
#define VFLOAT float32xm8_t
#define VLE_FLOAT vlev_float32xm8
#define VLSE_FLOAT vlsev_float32xm8
#define VSE_FLOAT vsev_float32xm8
#define VFMACC_VV vfmaccvv_float32xm8
#define VFREDSUM vfredsumvs_float32xm8

float *in_row_parallelism(float **a_new,  int ai,  int **a_d, int *size,
                            float *b,  int bi,
                            float *c,  int ci){
    size_t vl;
    VFLOAT va, vb, vtmp;
    int32xm8_t vd;
    float* ans=c;
    for(int i=0; i<ai; i++){
        float* tmp_a=a_new[i];
         int* tmp_d=a_d[i];
        int s=size[i];
        float* tmp_b=b;
        float* test=new float[ai];
        vtmp=vfmvvf_float32xm8(0.0, VSETVLMAX);
        for(int k=s; k>0; k-=vl){
            vl=VSETVL(k, RVV_E32, RVV_M8);
            va=VLE_FLOAT(tmp_a, vl);
            vd=vlev_int32xm8(tmp_d, vl);
            vd=vsllvi_int32xm8(vd, 2, vl);
            vb=vlxev_float32xm8(tmp_b, vd, vl); 
            vtmp=VFMACC_VV(vtmp, va, vb, vl);
            tmp_a+=vl;
            tmp_b+=vl;
        } 
        
    // vfloat32m1_t vans=vfmvvf_float32xm1(0.0, 1);
    float32xm8_t vans;
    vans=VFREDSUM(vtmp, vans, VSETVLMAX);
    vsev_float32xm8(ans, vans, 1);
    ans++;
    }   
    return c;
}

float *cross_line_parallelism(float *a,  int ai,  int aj,
                             int *a_d,  int adi,  int adj,
                            float *b,  int bi,
                            float *c,  int ci){
    size_t vl;
    VFLOAT va, vb, vtmp, vc;
    int32xm8_t vidx;
    e32xm8_t mask;
    float* tmp_b=b;
    for(int i=0; i<ai; i++){
        float *tmp_a=a+i*aj;
         int *tmp_i=a_d+i*aj;
        float *tmp_b=b;
        float *tmp_c=c;
        vtmp=vfmvvf_float32xm8(0.0, VSETVLMAX);
        for(int i=0; i<ai; i+=vl){
            vl=VSETVL(ai-i, RVV_E32, RVV_M8);
            va=VLE_FLOAT(tmp_a, vl);
            vidx=vlev_int32xm8(tmp_i, vl);
            mask=vmsnevi_e32xm8_int32xm8(vidx, 0, vl);
            vidx=vsubvx_int32xm8(vidx, 1, vl);
            vidx=vsllvi_int32xm8(vidx, 2, vl);
            vb=vlxev_mask_float32xm8(vb, tmp_b, vidx, mask, vl);
            vtmp=vfmulvv_mask_float32xm8(vtmp, va, vb, mask, vl);

            vc=vlxev_mask_float32xm8(vc, tmp_c, vidx, mask, vl);
            vtmp=vfaddvv_mask_float32xm8(vtmp, vtmp, vc, mask, vl);
            vsxev_mask_float32xm8(tmp_c, vidx, vtmp, mask, vl);
            tmp_a+=vl;
            tmp_i+=vl;
            tmp_b+=vl;
            tmp_c+=vl;
        }
    }
    return c;
}

void print(float *a,  int ai){
    for(int i=0; i<ai; i++){

            cout << a[i] << " ";
        
    }
    cout << endl;
}

void print(int *a,  int ai){
    for(int i=0; i<ai; i++){

            cout << a[i] << " ";
        
    }
    cout << endl;
}

void print(float *a,  int ai,  int aj){
    for(int i=0; i<ai; i++){
        for(int j=0; j<aj; j++){
            cout << a[i*aj+j] << " ";
        }
        cout << endl;
    }
}
void print(int *a,  int ai,  int aj){
    for(int i=0; i<ai; i++){
        for(int j=0; j<aj; j++){
            cout << a[i*aj+j] << " ";
        }
        cout << endl;
    }
}

void print(float **a,  int ai, int *size){
    for(int i=0; i<ai; i++){
        for(int j=0; j<size[i]; j++){
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}
void print(int **a,  int ai, int *size){
    for(int i=0; i<ai; i++){
        for(int j=0; j<size[i]; j++){
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}
// void initiation(float *a,  int ai,  int aj, float *v,  int vi){
//     for(int i=0; i<ai; i++){
//         for(int j=0; j<aj; j++){
//             a[i*aj+j]=distribution(generator);
//         }
//     }

//     for(int i=0; i<vi; i++){
//         v[i]=distribution(generator);
//     }

// }
void zero(float *a,  int ai){
    for(int i=0; i<ai; i++){
        if(a[i]>-1.0&&a[i]<1.0) a[i]=0.0;
    }
}
void zero(float *a,  int ai,  int aj){
    for(int i=0; i<ai; i++){
        for(int j=0; j<aj; j++){
            if(a[i*aj+j]>-1.0&&a[i*aj+j]<1.0) a[i*aj+j]=0.0;
        }
    }
}
void permute(float *a,  int ai,  int aj, float **a_new,  int **a_d, int *size){

    for(int i=0; i<ai; i++){
        int nnz=0;
        for(int j=0; j<aj; j++){
            if(a[i*aj+j]!=0.0){
                nnz++;
            }
        }
        float *tmp=new float[nnz];
         int *tmp_d=new  int[nnz];
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
void col_major(float **a, float *b,  int bi,  int bj){
    for(int i=0; i<bi; i++){
        for(int j=0; j<bj; j++){
            if(a[j][i]!=0.0)
            b[i*bj+j]=a[j][i];
            else b[i*bj+j]=0.0;
        }
    }
}
void col_major( int **a_d,  int a_di, int *size,  int *b,  int bi,  int bj){
    for(int j=0; j<bj; j++){
        for(int i=0; i<bi; i++){
            if(i<size[j])
            b[i*bj+j]=a_d[j][i]+1;
            else b[i*bj+j]=0;
        }
    }

}

int main(int argc, char*argv[]){
    if(argc<2) cout << "please type ./parallelism mean deviation dimension" << endl;


    normal_distribution<float> distribution(atof(argv[1]), atof(argv[2]));

    int n=atoi(argv[3]);
    float *a=new float[n*n];
    float *a_col_major=new float[n*n];
    float **a_new=new float*[n];
     int **a_d=new  int*[n];
     int *a_col_d=new  int[n*n];
    int *size=new int[n];
    float *b=new float[n];
    float *c_in_row=new float[n];
    float *c_cross_line=new float[n];

    // initiation(a,n,n,b,n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            a[i*n+j]=distribution(generator);
        }
    }

    for(int i=0; i<n; i++){
        b[i]=distribution(generator);
    }

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

    float *dis=new float[n];
    for(int i=0; i<n; i++){
        float tmp=0.0;
        for(int j=0; j<size[i]; j++){
            if(a[i*n+j]!=0.0) tmp+=1.0;
        }
        dis[i]=tmp/float(size[i]);
    }
    print(dis, n);
    cout << "in_row" << " " << duration_in_row.count() << endl;
    cout << "cross_line" << " " << duration_cross_line.count() << endl;

    return 0;
}