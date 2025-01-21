#include<iostream>
#include<riscv_vector.h>
#include"ELLPACK_P.h"
#include<random>
#include"utilis.h"
#include<string.h>

default_random_engine gen;
uniform_real_distribution<double> dis(-10,10);
using namespace std;

void SpMSpM(ELLPACK_P M1, ELLPACK_P M2, ELLPACK_P M3){
    
    if(M3._rows != M1._rows) {
        cerr << "M3 rows dismatches M1 rows." << endl; 
    }

    Idxcol = 0;
    Idx*tmp_nnz = new Idx[M2._rows];
    Idx*tmp_row = new Idx[M2._rows];
    Idx**tmp_col = new unsigned int* [M2._rows];
    Val **tmp_v = new Val* [M2._rows];

    for(Idxi=0; i<M2._rows; i++){
        col = max(col, M2._coordPtr[i][M2._nnz_num[i]-1]);
        tmp_row[i] = i;
        tmp_col[M2._rowidx[i]] = M2._coordPtr[i];
        tmp_v[M2._rowidx[i]] = M2._rowPtr[i];
        tmp_nnz[M2._rowidx[i]] = M2._nnz_num[i];
    }
    delete [] M2._rowPtr;
    delete [] M2._coordPtr;
    delete [] M2._rowidx;
    delete [] M2._nnz_old_num;
    
    M2._rowidx = tmp_row;
    M2._coordPtr = tmp_col;
    M2._rowPtr = tmp_v;
    M2._nnz_old_num = M2._nnz_num;
    M2._nnz_num = tmp_nnz;
    
    // register M3 space
    col++;
    for(Idxi=0; i<M3._rows; i++){
        M3._nnz_old_num[i] = col;
        Val *tmp_c_v = new Val [col]{0.0f};
        Idx*tmp_c_col = new Idx[col]{0};
        M3._rowidx[i] = i;
        M3._rowPtr[i] = tmp_c_v;
        M3._coordPtr[i] = tmp_c_col;
    }

    size_t vl;
    vfloat32m8_t vb, vc, vab, vtmp;
    vuint32m8_t vbidx, vcidx;         

    for(Idxi=0; i<M1._rows; i++){

        Idxb_row = M2._rowidx[i];
        
        for(Idxj=0; j<M1._nnz_num[i]; j++){
            
            Val tmp_a_v = M1._rowPtr[i][j];
            Idxtmp_a_row = M1._rowidx[i];
            Idxtmp_a_col = M1._coordPtr[i][j];
            
            Val *tmp_b_v = M2._rowPtr[tmp_a_col];
            Idx*tmp_b_col = M2._coordPtr[tmp_a_col];

            Val *tmp_c_v = M3._rowPtr[tmp_a_row];
            Idx*tmp_c_col = M3._coordPtr[tmp_a_row];

            for(Idxl=M2._nnz_num[tmp_a_col]; l>0; l-=vl){
                vl = vsetvl_e32m8(l);
                vb = vle32_v_f32m8(tmp_b_v, vl);
                vbidx = vle32_v_u32m8(tmp_b_col, vl);
                vcidx = vmv_v_v_u32m8(vbidx, vl);
                vbidx = vsll_vx_u32m8(vbidx, 2, vl);
                vc = vluxei32_v_f32m8(tmp_c_v, vbidx, vl);
                vab = vfmul_vf_f32m8(vb, tmp_a_v, vl);               
                vc = vfadd_vv_f32m8(vc, vab, vl);
                vsuxei32_v_u32m8(tmp_c_col, vbidx, vcidx, vl);
                vsuxei32_v_f32m8(tmp_c_v, vbidx, vc, vl);
                
                tmp_b_v += vl;
                tmp_b_col += vl;
            }
        }
    }
    for(Idxi=0; i<M3._rows; i++){
        Idxnnz=0;
        for(Idxj=0; j<col; j++){
            if(M3._rowPtr[i][j]!=0){
                nnz++;
            }
        }
        M3._nnz_old_num[i] = nnz;
        M3._nnz_num[i] = M3._nnz_old_num[i];

        int tmp=0;
        Val *tmp_v = new Val [nnz]{0.0f};
        Idx*tmp_col = new Idx[nnz]{0};
        for(Idxj=0; j<col; j++){
            if(M3._rowPtr[i][j]!=0){
                tmp_v[tmp]=M3._rowPtr[i][j]; 
                // cout << M3._rowPtr[i][j] << " ";
                tmp_col[tmp]=M3._coordPtr[i][j]; 
                // cout << M3._coordPtr[i][j] << " ";
                tmp++;
            }
        } 

        delete [] M3._rowPtr[i];
        delete [] M3._coordPtr[i];
        M3._rowPtr[i] = tmp_v;
        M3._coordPtr[i] = tmp_col;
    }
    M3.print_nnz();    
}

void scalar_normal(Val *a, const Idxai, const Idxaj,
                     Val *b, const Idxbi, const Idxbj,
                     Val *c, const Idxci, const Idxcj){
    Idxm=ai;
    Idxp=aj;
    Idxn=bj;

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            for(int k=0; k<p; k++)
            c[i*n+j]+=a[i*p+k]*b[k*n+j];
        }
    }   
}

int main(int argc, char **argv){

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <integer>" << endl;
        return 1;
    }

    int dim = atoi(argv[1]);
    Idxd = static_cast<unsigned int>(dim);
    
    ELLPACK_P M1(d), M2(d), M3(d);

    Val *a = new Val[d*d];
    Val *b = new Val[d*d];
    Val *c = new Val[d*d];
    for(int i=0; i<d; i++){
        for(int j=0; j<d; j++){
            if(rand()&127!=0){// generate 0
                a[i*d+j]=dis(gen);
                b[i*d+j]=dis(gen);

            }
            else{
                a[i*d+j]=0.0;
                b[i*d+j]=0.0;
            }   
        }
    }

    M1.ELLPACKinitWithDenseMatrix(a, d, d);
    M2.ELLPACKinitWithDenseMatrix(b, d, d);

    M1.ELLPACK_P_vreorder();
    M2.ELLPACK_P_vreorder();  

    scalar_normal(a, d, d, b, d, d, c, d, d);
    SpMSpM(M1, M2, M3);

    ELLPACK_P M4(d);
    M4.ELLPACKinitWithDenseMatrix(c, d, d);
    M4.print_ELLPACK_P();


    M3.ELLPACK_P_vreorder();
    M4.ELLPACK_P_vreorder();
    
    M3.print_ELLPACK_P();
    M4.print_ELLPACK_P();
    
    return 0;
}
