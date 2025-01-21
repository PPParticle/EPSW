/*
 * ELLPACK.h
 *
 *  Created on: May 6, 2023
 *      Author: fcqiao
 */
#ifndef ELLPACK_P_H_
#define ELLPACK_P_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <riscv-vector.h>

#define Val float
using namespace std;

struct ELLPACK_P
{
public:
    /* values and coord has identical organization */
    Val **_rowPtr;
    unsigned int **_coordPtr;

    /* the number of values of each row, int8/int16 */
    unsigned int *_nnz_num;
    unsigned int *_nnz_old_num;

    unsigned int _nnz, _rows, _cols;

    unsigned int *_rowidx;

    /* strip accomodating cache size */
    unsigned int _strips;

    ELLPACK_P()
    {
        this->_rowPtr = NULL;
        this->_coordPtr = NULL;
        this->_nnz_num = NULL;
        this->_nnz_old_num = NULL;
        this->_rows = 0;
        this->_cols = 0;
        this->_rowidx = NULL;
        this->_nnz = 0;
        this->_strips = NULL;
    }

    ELLPACK_P(unsigned int rows, unsigned int cols)
    {
        this->_rows = rows;
        this->_cols = cols;
        this->_rowPtr = new Val*[this->_rows];
        this->_coordPtr = new unsigned int*[this->_rows];
        this->_nnz_num = new unsigned int[this->_rows];
        this->_nnz_old_num = new unsigned int[this->_rows];
        this->_rowidx = new unsigned int[this->_rows];

        for(int i=0; i<this->_rows; i++){
            this->_nnz_old_num[i] = 0;
            this->_nnz_num[i] = this->_nnz_old_num[i];
            this->_rowidx[i] = i;
            this->_rowPtr[i] = nullptr;
            this->_coordPtr[i] = nullptr;
        }

        this->_nnz = 0;
        this->_strips = 0;

    }

    void init(unsigned int rows, unsigned int cols, Val **rowPtr, unsigned int **coordPtr, unsigned int *nnz_num, unsigned int *nnz_old_num,
              unsigned int *rowidx, unsigned int nnz, unsigned int strips){
        this->_rows = rows;
        this->_cols = cols;
        this->_rowPtr = rowPtr;
        this->_coordPtr = coordPtr;
        this->_nnz_num = nnz_num;
        this->_nnz_old_num = nnz_old_num;
        this->_rowidx = rowidx;
        this->_nnz = nnz;
        this->_strips = strips;
    }

    ELLPACK_P(unsigned int rows, unsigned int cols, Val **rowPtr, unsigned int **coordPtr, unsigned int *nnz_num, unsigned int *nnz_old_num,
              unsigned int *rowidx, unsigned int nnz, unsigned int strips)
    {
        init(rows, cols, rowPtr, coordPtr, nnz_num, nnz_old_num, rowidx, nnz, strips);
    }

    void CSRinitWithDenseMatrix(const Val *fvalues, const int rows, const int cols, 
    int* csr_v, int* csr_e, Val* csr_ev){
        csr_v = new int[rows+1];
        csr_v[0] = 0;
        int cnt=0;
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                if(fvalues[i*cols+j]!=0.0){
                    cnt++;
                }
            }
            csr_v[i+1] = cnt;
        }
        csr_e = new int [cnt];
        csr_ev = new Val [cnt];
        cnt=0;
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                if(fvalues[i*cols+j]!=0.0){
                    csr_e[cnt] = j;
                    csr_ev[cnt++] = fvalues[i*cols+j];
                }
            }
        }   
    }

    void ELLPACKinitWithDenseMatrix(const Val *fvalues, const int rows, const int cols)
    {
        this->_rows = rows;
        this->_cols = cols;
        this->_rowPtr = new Val *[this->_rows];
        this->_coordPtr = new unsigned int*[this->_rows];
        this->_rowidx = new unsigned int [this->_rows];
        for (int i = 0; i < rows; i++)
        {
            int n = 0;
            for (int j = 0; j < cols; j++)
            {
                Val val = fvalues[i * cols + j];
                if (val == 0)
                {
                    continue;
                }
                n++;
            }
            this->_nnz+=n;
            this->_rowidx[i] = i;
            this->_nnz_old_num[i] = n;
            this->_rowPtr[i] = new Val[n];
            this->_coordPtr[i] = new unsigned int[n];
        } 

        for (int i = 0; i < this->_rows; i++)
        {   
            this->_nnz_num[i] = this->_nnz_old_num[i];
            int tmp=0;
            for (int j = 0; j < cols; j++)
            {
                Val val = fvalues[i * cols + j];
                if (val == 0)
                {
                    continue;
                }
                _rowPtr[i][tmp] = val;
                _coordPtr[i][tmp] = j;
                tmp++;
            }
        } 
    }

    void ELLPACKinitWithCSR(const int* csr_v, const int* csr_e, const Val* csr_ev, 
        const unsigned int nr, const unsigned int nc, const unsigned int ne){
            this->_rows = nr;
	        this->_cols = nc;
            this->_nnz = 0;
            this->_rowPtr = new Val*[this->_rows];
            this->_coordPtr = new unsigned int*[this->_rows];
            this->_nnz_num = new unsigned int[this->_rows];
            this->_nnz_old_num = new unsigned int[this->_rows];
            this->_rowidx = new unsigned int[this->_rows];
            
            for(int i=1; i<nr+1; i++){
                this->_rowidx[i-1] = i-1;
                this->_nnz_num[i-1] = this->_nnz_num[i-1] = csr_v[i] - csr_v[i-1];
                Val *tmp_v = new Val[this->_nnz_num[i-1]];
                unsigned int *tmp_i = new unsigned int [this->_nnz_num[i-1]];
                this->_nnz += this->_nnz_num[i-1];
                this->_rowPtr[i-1] = tmp_v;
                this->_coordPtr[i-1] = tmp_i;
            }
            for(int i=1; i<nr+1; i++){
                for(int k=csr_v[i-1]; k<csr_v[i]; k++){
                    this->_coordPtr[i-1][k-csr_v[i-1]] = csr_e[k];
                    this->_rowPtr[i-1][k-csr_v[i-1]] = csr_ev[k];
                }
            }
    }

    void ELLPACKinitWithCSR(const int* csr_v, const int* csr_e, const int* csr_ev,
        const unsigned int nr, const unsigned int ne){
            this->_rows = nr;
            this->_nnz = 0;
            this->_rowPtr = new Val*[this->_rows];
            this->_coordPtr = new unsigned int*[this->_rows];
            this->_nnz_num = new unsigned int[this->_rows];
            this->_nnz_old_num = new unsigned int[this->_rows];
            this->_rowidx = new unsigned int[this->_rows];

            for(int i=1; i<nr+1; i++){
                this->_rowidx[i-1] = i-1;
                this->_nnz_num[i-1] = this->_nnz_num[i-1] = csr_v[i] - csr_v[i-1];
                Val *tmp_v = new Val[this->_nnz_num[i-1]];
                unsigned int *tmp_i = new unsigned int [this->_nnz_num[i-1]];
                this->_nnz += this->_nnz_num[i-1];
                this->_rowPtr[i-1] = tmp_v;
                this->_coordPtr[i-1] = tmp_i;
            }

            // for(int i=0; i<this->_rows; i++){
            //     cout << this->_nnz_num[i] << " ";
            // }
            // cout << endl;

            for(int i=1; i<nr+1; i++){
                for(int k=csr_v[i-1]; k<csr_v[i]; k++){
                    // this->_rowPtr[i-1][k-csr_v[i-1]] = csr_ev[csr_v[k]];
                    this->_coordPtr[i-1][k-csr_v[i-1]] = csr_e[k];
                    this->_rowPtr[i-1][k-csr_v[i-1]] = csr_ev[k];
                }
            }
    }
    
     void EPSWToCSR(int* csr_v, int* csr_e, Val* csr_ev){
        csr_v[0]=0;
        for(int i=0; i<this->_rows; i++){
            csr_v[i+1] = this->_nnz_num[i];
            csr_v[i+1] += csr_v[i];
        }
        int tmp=0;
        for(int i=0; i<this->_rows; i++){
            for(int j=0; j<this->_nnz_num[i]; j++){
                    csr_ev[tmp]=this->_rowPtr[i][j];
                    csr_e[tmp++] = this->_coordPtr[i][j];
            }
        }
    }
    
    void reorder()
    {   
        unsigned int *tmp_rowidx = new unsigned int [this->_rows];
        unsigned int *tmp_nnz_old_num = new unsigned int[this->_rows];
        unsigned int **tmp_coordPtr = new unsigned int*[this->_rows];
        Val **tmp_rowPtr = new Val*[this->_rows]; 
        
        for (int i = 0; i < this->_rows; i++)
        {  
            tmp_nnz_old_num[i] = this->_nnz_old_num[i];
            tmp_rowidx[i] = this->_rowidx[i];
            tmp_coordPtr[i] = this->_coordPtr[i];
            tmp_rowPtr[i] = this->_rowPtr[i];
        }
        // _rows times scan tmp and pick the biggest one
        for (int i = 0; i < _rows; i++)
        {   
            int max_idx = 0;           // idx of the longest row
            for (int j = 0; j < _rows; j++)
            {
                if (tmp_nnz_old_num[j] > tmp_nnz_old_num[max_idx])
                {
                    max_idx = j;
                }
            } 
            tmp_nnz_old_num[max_idx] = 0;
            this->_rowidx[i] = tmp_rowidx[max_idx];
            this->_nnz_num[i] = _nnz_old_num[max_idx];
            this->_coordPtr[i] = tmp_coordPtr[max_idx];
            this->_rowPtr[i] = tmp_rowPtr[max_idx];
        }
        delete [] tmp_rowidx;
        delete [] tmp_nnz_old_num;
        delete [] tmp_rowPtr;
        delete [] tmp_coordPtr;
    }

unsigned int merge(const unsigned int *src1, const unsigned int n1, const unsigned int *src2, const unsigned int n2, unsigned int *dst)
    {
        size_t vl, k, num = 0;
        vl = 0;
        k = n1;
        // if(k>VLMAX){
        for (; k > 0; k -= vl)
        {
            vl = vsetvli(k, RVV_E32, RVV_M8);
            uint32xm8_t v1 = vlev_uint32xm8(src1, vl);
            vsev_uint32xm8(dst, v1, vl);
            src1 += vl;
            dst += vl;
        }
        num += n1;
        // }
        vl = 0;
        k = n2;
        for (; k > 0; k -= vl)
        {
            vl = vsetvli(k, RVV_E32, RVV_M8);
            uint32xm8_t v1 = vlev_uint32xm8(src2, vl);
            vsev_uint32xm8(dst, v1, vl);
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
                vl = vsetvli(k, RVV_E32, RVV_M8);
                uint32xm8_t vtmp = vlev_uint32xm8(tmp_s, vl);
                uint32xm8_t vidx = vlev_uint32xm8(tmp_i, vl);

                uint32xm8_t v_and = vandvx_uint32xm8(vtmp, logic, vl);
                e32xm8_t v_neq = vmseqvx_e32xm8_uint32xm8(v_and, 0, vl);
                e32xm8_t v_eq = vmnotm_e32xm8(v_neq, vl);
                size_t n_1 = vmpopcm_e32xm8(v_eq, vl);
                size_t n_0 = vmpopcm_e32xm8(v_neq, vl);
                // cout << n_1 << endl;

                uint32xm8_t v_1 = vcompressvm_uint32xm8_e32xm8(vtmp, v_eq, vl);
                uint32xm8_t v_0 = vcompressvm_uint32xm8_e32xm8(vtmp, v_neq, vl);
                vsev_uint32xm8(tmp_1, v_1, n_1);
                vsev_uint32xm8(tmp_0, v_0, n_0);

                tmp_0 += n_0;
                tmp_1 += n_1;

                v_1 = vcompressvm_uint32xm8_e32xm8(vidx, v_eq, vl);
                v_0 = vcompressvm_uint32xm8_e32xm8(vidx, v_neq, vl);
                vsev_uint32xm8(tmp_i_1, v_1, n_1);
                vsev_uint32xm8(tmp_i_0, v_0, n_0);

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

    void ELLPACK_P_vreorder()
    {   
        unsigned int *tmp_rowidx = new unsigned int[this->_rows];
        unsigned int* tmp_old_num = new unsigned int[this->_rows];
        Val **tmp_rowptr = new Val *[this->_rows];
        unsigned int **tmp_coordptr = new unsigned int*[this->_rows];

        for (unsigned int i = 0; i < this->_rows; i++)
        {
            tmp_rowidx[i] = this->_rowidx[i];
            this->_nnz_old_num[i] = this->_nnz_num[i];
            tmp_old_num[i] = this->_nnz_num[i];
            tmp_coordptr[i] = this->_coordPtr[i];
            tmp_rowptr[i] = this->_rowPtr[i];
        }

        unsigned int shift = 1;
        unsigned int logic = 1;
        for (; shift < 32; shift++, logic = logic << 1)
        {   
            radix_sort(tmp_old_num, this->_nnz_num, tmp_rowidx, this->_rows, logic);
            tmp_old_num = this->_nnz_num;
        }
        
        for (unsigned int i = 0; i < this->_rows; i++)
        {
            tmp_rowptr[i] = this->_rowPtr[tmp_rowidx[i]];
            tmp_coordptr[i] = this->_coordPtr[tmp_rowidx[i]];

        }

        delete[] this->_rowPtr;
        delete[] this->_coordPtr;
        delete[] this->_rowidx;

        this->_rowidx = tmp_rowidx;
        this->_rowPtr = tmp_rowptr;
        this->_coordPtr = tmp_coordptr;
    }

    void stripe_static(int stripe)
    {   
        // int strips = (this->_rows-this->_rows%stripe)/stripe;
        this->_strips = stripe;
    }
   
    vector<vector<unsigned int>> stripe_adaptive(){
        vector<vector<unsigned int>> stripes;
        for(int i=0; i<this->_rows;){
            vector<unsigned int> s;
            int cnt=1;
            int tmp = this->_nnz_num[i];
            s.push_back(this->_rowidx[i]);
            if(tmp>(this->_cols+2)/2){
                while(tmp>cnt*(this->_cols+2)/2){
                    if(cnt==10){
                        break;
                    }
                    i++;
                    tmp+=this->_nnz_num[i];
                    s.push_back(this->_rowidx[i]);
                    cnt++;
                    // cout << tmp << " " << cnt*(this->_cols+2)/2 << endl;
                }
                if(cnt>1&&cnt<10){
                    cnt--;
                    i--;
                    s.pop_back();
                }
            }else if(tmp==0){
                break;
            }
            stripes.push_back(s);
            i++;
        }
        // cout << "stripe success" << endl;
        // for(auto x:stripes){
        //     cout << x.size() << " " << endl;
        // }
        return stripes;
    };

    Val *spmv_rvv(const ELLPACK_P M, const Val *fvalues, const int rows);
    Val *spmv_scalar(const ELLPACK_P M, const Val *fvalues, const int rows);
    Val *spmm_rvv(ELLPACK_P M, Val *fvalues, const int rows, const int cols);
    Val *spmm_scalar_ELLP(ELLPACK_P M, Val *fvalues, const int rows, const int cols);
    Val *spmm_scalar_normal(ELLPACK_P M, Val *fvalues, const int rows, const int cols);


    void print_nnz(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_nnz_num[i] << " ";
        }
        cout << endl;
    }
    void print_nnz_old(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_nnz_old_num[i] << " ";
        }
        cout << endl;
    }
    void print_ELLPACK_P()
    {
        cout << "ELLPACK_P M: " << endl;
        for (unsigned int i = 0; i < this->_rows; i++)
        {
            for (unsigned int j = 0; j < this->_nnz_num[i]; j++)
                // cout << this->_rowPtr[i][j] << " ";
                printf("%4.2lf ",this->_rowPtr[i][j]);
            cout << endl;

            for (unsigned int j = 0; j < this->_nnz_num[i]; j++)
                // cout << this->_coordPtr[i][j] << " ";
                printf("%d ",this->_coordPtr[i][j]);
            cout << endl;
        }
    }
    
    void print_rowidx(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_rowidx[i] << " ";
        }
        cout << endl;
    }
    void print_rowPtr(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_rowPtr[i] << " ";
        }
        cout << endl;
    }
    void print_value(){
        for(int i=0; i<this->_rows; i++){
            for(int j=0; j<this->_nnz_num[i]; j++)
            cout << this->_rowPtr[i][j] << " ";
            cout << endl;
        }
    }
    void print_coordPtr(){
        for(int i=0; i<this->_rows; i++){
            for(int j=0; j<this->_nnz_num[i]; j++)
            cout << this->_coordPtr[i][j] << " ";
            cout << endl;
        }
    }
};
#endif