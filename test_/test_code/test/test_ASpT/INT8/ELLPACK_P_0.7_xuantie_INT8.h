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
#include <riscv_vector.h>

#define Val int8_t
#define Idx Idx
using namespace std;

struct ELLPACK_P
{
public:
    /* values and coord has identical organization */
    Val **_rowPtr;
    Idx **_coordPtr;

    /* the number of values of each row*/
    Idx *_nnz_num;
    Idx *_nnz_old_num;

    /* fault-only-first instructions are used to vectorize loops with data-dependent */
    Idx _nnz, _rows;

    Idx *_rowidx;

    /* strip accomodating cache size */
    Idx _strips;

    ELLPACK_P()
    {
        this->_rowPtr = NULL;
        this->_coordPtr = NULL;
        this->_nnz_num = NULL;
        this->_nnz_old_num = NULL;
        this->_rows = 0;
        this->_rowidx = NULL;
        this->_nnz = 0;
        this->_strips = 0;
    }

    ELLPACK_P(Idx n)
    {
        this->_rows = n;
        this->_rowPtr = new Val*[this->_rows];
        this->_coordPtr = new Idx*[this->_rows];
        this->_nnz_num = new Idx[this->_rows];
        this->_nnz_old_num = new Idx[this->_rows];
        this->_rowidx = new Idx[this->_rows];
        this->_nnz = 0;
        this->_strips = 0;
    }

    void init(Idxrows, Val **rowPtr, Idx **coordPtr, Idx *nnz_num, Idx *nnz_old_num,
              Idx *rowidx, Idx nnz, Idx strips){
        this->_rows = rows;
        this->_rowPtr = rowPtr;
        this->_coordPtr = coordPtr;
        this->_nnz_num = nnz_num;
        this->_nnz_old_num = nnz_old_num;
        this->_rowidx = rowidx;
        this->_nnz = nnz;
        this->_strips = strips;
    }

    ELLPACK_P(Idx rows, Val **rowPtr, Idx **coordPtr, Idx *nnz_num, Idx *nnz_old_num,
              Idx *rowidx, Idx nnz, Idx strips)
    {
        init(rows, rowPtr, coordPtr, nnz_num, nnz_old_num, rowidx, nnz, strips);
    }

    void ELLPACKinitWithDenseMatrix(const Val *fvalues, const Idx rows, const Idx cols)
    {
        this->_rows = rows;
        this->_rowPtr = new Val*[this->_rows];
        this->_coordPtr = new Idx*[this->_rows];
        this->_rowidx = new Idx[this->_rows];
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
            this->_rowidx[i] = i;
            this->_nnz_old_num[i] = n;
            this->_rowPtr[i] = new Val[n];
            this->_coordPtr[i] = new Idx[n];
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

    void print_nnz(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_nnz_num[i] << " ";
        }
        cout << endl;
    }
    void print_ELLPACK_P()
    {
        cout << "ELLPACK_P M: " << endl;
        for (Idxi = 0; i < this->_rows; i++)
        {
            for (Idxj = 0; j < this->_nnz_num[i]; j++)
                // cout << this->_rowPtr[i][j] << " ";
                printf("%8.2lf ", (int32_t)this->_rowPtr[i][j]);
            cout << "         ";
            for (Idxj = 0; j < this->_nnz_num[i]; j++)
                // cout << this->_coordPtr[i][j] << " ";
                printf("%d ", (int32_t)this->_coordPtr[i][j]);
            cout << endl;
        }
    }
    void print_rowidx(){
        for(int i=0; i<this->_rows; i++){
            cout << this->_rowidx[i] << " ";
        }
        cout << endl;
    }

    void reorder()
    {   
        Idx *tmp_rowidx = new Idx[this->_rows];
        Idx *tmp_nnz_old_num = new Idx[this->_rows];
        Idx **tmp_coordPtr = new Idx*[this->_rows];
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

    
Idx merge(const Idx *src1, const Idx n1, const Idx *src2, const Idx n2, Idx *dst)
    {
        size_t vl, k, num = 0;
        vl = 0;
        k = n1;
        for (; k > 0; k -= vl)
        {
            vl = vsetvl_e32m8(k);
            vuint32m8_t v1 = vle32_v_u32m8(src1, vl);
            vse32_v_u32m8(dst, v1, vl);
            src1 += vl;
            dst += vl;
        }
        num += n1;
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
void radix_sort(const Idx *src, Idx *dst, Idx *idx, const Idx n, Idx logic)
    {
        if (n > 1)
        {
            size_t vl;
            size_t n0 = 0, n1 = 0;
            Idx *dst_0 = new Idx[n];
            Idx *dst_1 = new Idx[n];
            Idx *idx_0 = new Idx[n];
            Idx *idx_1 = new Idx[n];

            Idx *tmp_0 = dst_0, *tmp_1 = dst_1;
            Idx *tmp_i_0 = idx_0, *tmp_i_1 = idx_1;

            const Idx *tmp_s = src;
            const Idx *tmp_i = idx;

            // cout << dst_1[1] << endl;
            for (Idx k = n; k > 0; k -= vl)
            {
                vl = vsetvl_e32m8(k);
                vuint32m8_t vtmp = vle32_v_u32m8(tmp_s, vl);
                vuint32m8_t vidx = vle32_v_u32m8(tmp_i, vl);

                vuint32m8_t v_and = vand_vx_u32m8(vtmp, logic, vl);
                vbool4_t v_neq = vmseq_vx_u32m8_b4(v_and, 0, vl);
                vbool4_t v_eq = vmnot_m_b4(v_neq, vl);
                size_t n_1 = vmpopc_m_b4(v_eq, vl);
                size_t n_0 = vmpopc_m_b4(v_neq, vl);
                // cout << n_1 << endl;

                vuint32m8_t v_1 = vcompress_vm_u32m8(v_eq, v_1, vtmp, vl);
                vuint32m8_t v_0 = vcompress_vm_u32m8(v_neq, v_0, vtmp, vl);
                vse32_v_u32m8(tmp_1, v_1, n_1);
                vse32_v_u32m8(tmp_0, v_0, n_0);

                tmp_0 += n_0;
                tmp_1 += n_1;

                v_1 = vcompress_vm_u32m8(v_eq, v_1, vidx,  vl);
                v_0 = vcompress_vm_u32m8(v_neq, v_0, vidx,  vl);
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
    void ELLPACK_P_vreorder()
    {
        Idx *tmp_rowidx = new Idx[this->_rows];
        Idx *tmp_old_num = new Idx[this->_rows];
        Val **tmp_rowptr = new Val*[this->_rows];
        Idx **tmp_coordptr = new Idx*[this->_rows];

        for (Idxi = 0; i < this->_rows; i++)
        {
            tmp_rowidx[i] = this->_rowidx[i];
            tmp_old_num[i] = this->_nnz_old_num[i];
            tmp_coordptr[i] = this->_coordPtr[i];
            tmp_rowptr[i] = this->_rowPtr[i];
        }

        Idxshift = 1;
        Idxlogic = 1;
        for (; shift < 10; shift++, logic = logic << 1)
        {
            radix_sort(tmp_old_num, this->_nnz_num, tmp_rowidx, this->_rows, logic);
            tmp_old_num = this->_nnz_num;
        }
        
        for (Idxi = 0; i < this->_rows; i++)
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

    // 待做
    void stripe(int strip)
    {
        _strips = strip;

        // the number of striped ELLPACK rows
        // int row_s = ceil(float(_rows)/float(_strips));
    }
   
    

    Val *spmv_rvv(const ELLPACK_P M, const Val *fvalues, const int rows);
    Val *spmv_scalar(const ELLPACK_P M, const Val *fvalues, const int rows);
    Val *spmm_rvv(ELLPACK_P M, Val *fvalues, const int rows, const int cols);
    Val *spmm_scalar_ELLP(ELLPACK_P M, Val *fvalues, const int rows, const int cols);
    Val *spmm_scalar_normal(ELLPACK_P M, Val *fvalues, const int rows, const int cols);

    
};
#endif