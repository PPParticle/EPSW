#include <iostream>
#include <riscv-vector.h>
#include <cstring>
#include <vector>
#include <algorithm>
#include <optional>

typedef unsigned char uchar;

template <typename rvv_T>
void print_rvv(rvv_T rvv, char* name, size_t vl){
    uint32_t tmp[vl];
    vsev_uint32xm8(tmp, rvv, vl);
    std::cout << name << " ";
    for(int i=0; i<vl; i++){
        std::cout << (int32_t)tmp[i] << " ";
    }   
    std::cout << std::endl;
}

// c906 doesn't support vector mask load/store instructions
// template <typename rvv_T>
// void print_rvv_m(rvv_T rvv, char* name, size_t vl){
//     uchar tmp[vl];
//     __builtin_riscv_vsmbool4(tmp, rvv);
//     std::cout << name << " ";
//     for(int i=0; i<vl; i++){
//         std::cout << (int32_t)tmp[i] << " ";
//     }   
//     std::cout << std::endl;
// }

void Partition1_rvv(uint32_t A[], size_t size) {
    size_t vlmax = vsetvli_max(RVV_E32, RVV_M8);
    size_t vl = 0;
    const int pivot = A[0];
    uint32xm8_t vlow, vhigh, vnum, vans, vtmp, v1;
    int slide = 0;
    e32xm8_t mask, mask_n;
    size_t pop_cnt[(size+vlmax)/vlmax];

    for(size_t i=0; i<size; i+=vl) {
            vl = vsetvli(size-i, RVV_E32, RVV_M8);
            vlow = vmvsx_uint32xm8(0, vl);
            vhigh = vmvsx_uint32xm8(0, vl);
            vnum = vlev_uint32xm8(A + i, vl);

            mask =  vmsleuvx_e32xm8_uint32xm8(vnum, pivot, vl);
            vlow = vcompressvm_uint32xm8_e32xm8(vnum, mask, vl); // less than pivot
            print_rvv(vlow, "vlow", vl);
            slide = vmpopcm_e32xm8(mask, vl); 
            pop_cnt[i/vlmax] = slide;
            
            mask_n = vmnotm_e32xm8(mask, vl);
            vhigh = vcompressvm_uint32xm8_e32xm8(vnum, mask_n, vl); // greater than pivot
            vhigh = vslideupvx_uint32xm8(vhigh, slide, vl);
            print_rvv(vhigh, "vhigh", vl);

            // vnum = vmergevvm_mask_int32xm8(vlow, vhigh, mask_n, vl);
            mask = vmsetm_e32xm8(6);
            v1 = vmvsx_uint32xm8(0, vl);
            vtmp = vmvsx_uint32xm8(0, vl);
            vtmp = vlev_mask_uint32xm8(v1, A+i, mask, vl);
            print_rvv(vtmp, "vmset", vl);
            vans = vrgathervv_mask_uint32xm8(vans, vhigh, vlow, mask, vl);
            vsev_uint32xm8(A + i, vnum, vl);      
            print_rvv(vans, "vans", vl);
    }

    // for(size_t i=0; i<(size+vlmax)/vlmax; i++) {std::cout << pop_cnt[i] << " ";}
    // std::cout << std::endl;

    int32_t B[size];
    size_t idx_b = 0;
    size_t pop_a = 0;
    int32xm8_t vmv;
    for(size_t i=0; i<(size+vlmax)/vlmax; i++){
        unsigned int idx = i*vlmax + pop_cnt[i];
        pop_a += pop_cnt[i];
        memcpy(B + idx_b, A + idx, vlmax - pop_cnt[i]);
        idx_b += vlmax - pop_cnt[i];
        vmv = vlev_uint32xm8(A + (i+1)*vlmax, pop_cnt[i+1]);
        vsev_uint32xm8(A + pop_a, vmv, pop_cnt[i+1]);
    }
    memcpy(A+pop_a, B, size-pop_a);
}


struct PrintElement {
    void operator()(int x) const {
        std::cout << x << " ";
    }
};

int main(){
    std::vector<uint32_t> A = {6,7,3,1,9,5,8,4,2};
    std::for_each(A.begin(), A.end(), [](int x){std::cout << x << " ";});
    std::cout << std::endl;
    Partition1_rvv(A.data(), A.size());
    std::for_each(A.begin(), A.end(), [](int x){std::cout << x << " ";});
    std::cout << std::endl;
    return 0;
}