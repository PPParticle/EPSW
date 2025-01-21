#include<iostream>
#include<riscv_vector.h>
using namespace std;

// http://pllab.cs.nthu.edu.tw/~jklee/papers/ICPPEMS2022.pdf
// icpp'20 plus_scan arr[0, 1, 2, 3, 4, 5, 6]
// output arr[0, 0-1, 0-2, 0-3, 0-4, 0-5, 0-6]
void plus_scan_ui (int n , unsigned int * src ) {
    size_t vl ;
    size_t vlmax = vsetvlmax_e32m1 () ;
    unsigned int carry = 0;
    vuint32m1_t x, y, vec_zero;
    vec_zero = vmv_v_x_u32m1 (0 , vlmax ) ;
    for (; n > 0; n -= vl ) {
        vl = vsetvl_e32m1 ( n ) ;
        x = vle32_v_u32m1 ( src , vl ) ;
        for ( size_t offset = 1; offset < vl ; offset = offset << 1) {
            // cout << " offset is "  << offset << endl;
            y = vslideup_vx_u32m1 ( vec_zero, x, offset, vl) ;
            x = vadd_vv_u32m1 (x , y , vl ) ;
        }
        x = vadd_vx_u32m1 (x , carry , vl ) ;
        vse32 ( src , x , vl ) ;
        carry = src [ vl - 1];
        src += vl ;
    }
}

// icpp'20 split radix sort
// code is not complete 
// lack p_add(), p_select(), permute(), get_flags()
// enumerate is an operation that takes a vector of flags as input and outputs a vector of integer that the value
// of an element is set to i if its flag is the i^{th} true flag.
unsigned int enumerate(int n, unsigned int *flags, unsigned int *dst, bool setBit){
    size_t vl;
    unsigned int count=0;
    unsigned int carry;
    for(; n>0; n-=vl){
        vl=vsetvl_e32m1(n);
        vuint32m1_t v=vle32_v_u32m1(flags, vl);
        vbool32_t mask = vmseq(v, setBit, vl);
        v=viota_m_u32m1(mask, vl);
        v=vadd(v, count, vl);
        vse32(dst, v, vl);
        count+=vcpop(mask, vl);
        flags+=vl;
        dst+=vl;
    }
    return count;
}

void p_add(int n, unsigned int* src, unsigned int x){
    size_t vl;
    for(; n>0; n-=vl){
        vl=vsetvl_e32m1(n);
        vuint32m1_t va=vle32_v_u32m1(src, vl);
        va=vadd(va, x, vl);
        vse32(src, va, vl);
        src+=vl;
    }
}

void permute(int n, unsigned int* src, unsigned int* dst, unsigned int* index){
    size_t vl;
    for(; n>0; n-=vl){
        vl=vsetvl_e32m1(n);
        vuint32m1_t x=vle32_v_u32m1(src, vl);
        vuint32m1_t vidx=vle32_v_u32m1(index, vl);
        vidx = vsll(vidx, 2, vl);
        vsuxei32(dst, vidx, x, vl);
        src+=vl;
        index+=vl;
    }
}

// the split is a form of permutation that uses a vector of flags to split the source vector into two halves.
void split(int n, unsigned int *src, unsigned int *dst, unsigned int *flags){
    unsigned int *i_up=new unsigned int[n];
    unsigned int *i_down=new unsigned int[n];
    unsigned int count=enumerate(n, flags, i_up, 0);
    enumerate(n, flags, i_down, 1);
    p_add(n, i_down, count);
    p_select(n, flags, i_down, i_up);
    permute(n, src, dst, i_up);
    delete [] i_up;
    delete [] i_down;
}

// split_radix_sort 基数排序函数
void split_radix_sort ( int n , unsigned int * src ) {
    unsigned int * buffer=new unsigned int[n];
    unsigned int * flags=new unsigned int[n];
    for ( int i = 0; i < 32; i ++) {
        get_flags (n , src , flags , i ) ;
        split (n , src , buffer , flags ) ;
        // swap src and buffer
        unsigned int * tmp = src ;
        src = buffer ;
        buffer = tmp ;
        }
    delete [] buffer;
    delete [] flags;
    return;
}

// icpp'20 segmented scan 
// it can be primitives for quick sort
// 3 segment-descriptor : head-flags, lengths, and head-pointers
// head-flags : an array of flags to indicate the starting index of each segment
// lengths : the length of each segments
// head-pointers : an array of pointers pointing to the starting element of each segment
// segmented plus-scan
void seg_plus_scan_ui(int n , unsigned int *src , unsigned int *head_flags){
    size_t vl;
    size_t vlmax=vsetvlmax_e32m1();
    unsigned int carry=0;
    vuint32m1_t x, y, vec_zero, vec_one;
    vuint32m1_t flags, flags_slideup;
    vbool32_t mask, carry_mask;
    vec_zero=vmv_v_x_u32m1(0, vlmax);
    vec_one=vmv_v_x_u32m1(1, vlmax) ;
    for (; n > 0; n-=vl){
        vl=vsetvl_e32m1(n);
        x=vle32_v_u32m1(src, vl);
        flags=vle32_v_u32m1(head_flags, vl);
        mask=vmsne_vx_u32m1_b32(flags, 0, vl);
        carry_mask=vmsbf(mask, vl); // set 1 before the first 1 in vsr
        flags=vmv_s_x_u32m1(flags, 1, vl); // only modify the first element of the destination vector
        for(size_t offset=1; offset < vl; offset=offset << 1){
            mask=vmsne_vx_u32m1_b32(flags, 1, vl);
            y=vslideup_vx_u32m1(vec_zero, x, offset, vl);
            x=vadd_vv_u32m1_m(mask, x, x, y, vl);
            flags_slideup=vslideup_vx_u32m1(vec_one, flags, offset, vl);
            flags=vor_vv_u32m1(flags, flags_slideup, vl);
        }
        x=vadd_vx_u32m1_m(carry_mask, x, x, carry, vl);
        vse32(src, x, vl);
        carry=src[vl-1];
        src += vl ;
        head_flags += vl ;
    }
}

int main(){
    unsigned int src[9]={1,2,3,4,5,6,7,8,9};
    for(auto x:src)
    cout << x << " ";
    cout << endl;
    plus_scan_ui(9,src);
    for(auto x:src)
    cout << x << " ";
    cout << endl;


}