// c910 长指令和短指令组合加载数据和存储数据的测试情况
// 总是 长指令速度>短指令组合 
#include <iostream>
#include <riscv_vector.h>
#include <assert.h>
#include <random>
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <malloc.h>

// Bytes
#define cache_v 32 * 1024 
#define num 128

using namespace std;
int main(int argc, char **argv)
{   

    size_t n = cache_v * atoi(argv[1]);
    cout << n << " bytes data" << endl;
    n = n/4;

    ifstream f("test_vl.txt", std::ios::binary);
    f.seekg(0, std::ios::end);
    std::streamsize size = f.tellg(); // signed integer
    _Float32* buffer = new _Float32[n]{0};
    int64_t num_t = 0;
    // 重置读取位置到文件开头
    f.seekg(0, std::ios::beg);
    _Float32 number;
    while (f >> number && num_t < n) {
        buffer[num_t++]=static_cast<_Float32>(number); // 将浮点数添加到数组中
    }
    f.close(); // 关闭文件
    // struct mallinfo mi = mallinfo();
    // cout << "Maximum space allocated from system: " << mi.hblkhd << " bytes" << endl;

    size_t vl, vl_max;
    vfloat32m8_t v8, vtmp;
    vfloat32m4_t v40, v41, vtmp40, vtmp41;
    vfloat32m2_t v20, v21, v22, v23, vtmp20, vtmp21, vtmp22, vtmp23;
    vfloat32m1_t v10, v11, v12, v13, v14, v15, v16, v17, 
                vtmp10, vtmp11, vtmp12, vtmp13, vtmp14, vtmp15, vtmp16, vtmp17;
    _Float32 *e32m8, *e32m4, *e32m2, *e32m1;



    e32m8 = new _Float32[n];
    _Float32 *tmp_8 = e32m8;
    _Float32 *tmp = &tmp[0];
    auto start_e32m8 = chrono::high_resolution_clock::now();
    vtmp = vfmv_v_f_f32m8 (0, vsetvlmax_e32m8());
    for (int j = 0; j < num; j++)
    {   
        tmp_8 = e32m8;
        tmp =  &buffer[0];
        for (int i = n; i > 0; i -= vl)
        {
            vl = vsetvl_e32m8(i);
            v8 = vle32_v_f32m8( tmp, vl);
            vtmp = vfmul_vf_f32m8(v8, tmp[0], vl);
            vse32_v_f32m8( tmp_8, vtmp, vl);
            tmp_8 += vl;
            tmp += vl;
        }
    }
    auto stop_e32m8 = chrono::high_resolution_clock::now();
    auto duration_e32m8 = chrono::duration_cast<chrono::microseconds>(stop_e32m8 - start_e32m8);
    cout << "e32m8 " << duration_e32m8.count()/num << endl;
    delete []e32m8;

    sleep(1);
    e32m4 = new _Float32[n];
    _Float32 *tmp_4 = e32m4;
    vl_max = vsetvlmax_e32m4();
    auto start_e32m4 = chrono::high_resolution_clock::now();
    vtmp40 = vfmv_v_f_f32m4 (0, vsetvlmax_e32m4());
    vtmp41 = vfmv_v_f_f32m4 (0, vsetvlmax_e32m4());
    for (int j = 0; j < num; j++)
    {   
        tmp_4 = e32m4;
        tmp = &buffer[0];
        for (int i = n; i > 0;) 
        {   
            if (i - 2 * vl_max > 0)
            {   
                v40 = vle32_v_f32m4(tmp, vl);
                v41 = vle32_v_f32m4((tmp + vl), vl);
                vtmp40 = vfmul_vf_f32m4(v40,  tmp[0], vl);
                vtmp41 = vfmul_vf_f32m4(v41,  tmp[0], vl);
                vse32_v_f32m4( tmp_4, vtmp40, vl);
                vse32_v_f32m4( (tmp_4 + vl), vtmp41, vl);
                tmp_4 += 2 * vl;
                tmp += 2 * vl;
                i -= 2 * vl;
            }
            else
            {   
                vl = vsetvl_e32m4(i);
                v40 = vle32_v_f32m4( tmp, vl);
                vtmp40 = vfmul_vf_f32m4(v40,  tmp[0], vl);
                vse32_v_f32m4( tmp_4, vtmp40, vl);
                tmp_4 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m4 = chrono::high_resolution_clock::now();
    auto duration_e32m4 = chrono::duration_cast<chrono::microseconds>(stop_e32m4 - start_e32m4);
    cout << "e32m4 " << duration_e32m4.count()/num << endl;
    delete []e32m4;

    sleep(1);
    e32m2 = new _Float32[n];
    _Float32 *tmp_2 = e32m2;
    vl_max = vsetvlmax_e32m2();
    auto start_e32m2 = chrono::high_resolution_clock::now();
    vtmp20 = vfmv_v_f_f32m2 (0, vsetvlmax_e32m2());
    vtmp21 = vfmv_v_f_f32m2 (0, vsetvlmax_e32m2());
    vtmp22 = vfmv_v_f_f32m2 (0, vsetvlmax_e32m2());
    vtmp23 = vfmv_v_f_f32m2 (0, vsetvlmax_e32m2());
    vl_max = vsetvlmax_e32m2();
    for (int j = 0; j < num; j++)
    {   
        tmp_2 = e32m2;
        tmp =  &buffer[0];
        for (int i = n; i > 0;)
        {   
            if (i - 4 * vl_max > 0)
            {   
                v20 = vle32_v_f32m2( tmp, vl);
                v21 = vle32_v_f32m2( (tmp + vl), vl);
                v22 = vle32_v_f32m2( (tmp + 2 * vl), vl);
                v23 = vle32_v_f32m2( (tmp + 3 * vl), vl);
                vtmp20 = vfmul_vf_f32m2(v20,  tmp[0], vl);
                vtmp21 = vfmul_vf_f32m2(v21,  tmp[0], vl);
                vtmp22 = vfmul_vf_f32m2(v22,  tmp[0], vl);
                vtmp23 = vfmul_vf_f32m2(v23,  tmp[0], vl);
                vse32_v_f32m2( tmp_2, vtmp20, vl);
                vse32_v_f32m2( (tmp_2 + vl), vtmp21, vl);
                vse32_v_f32m2( (tmp_2 + 2 * vl), vtmp22, vl);
                vse32_v_f32m2( (tmp_2 + 3 * vl), vtmp23, vl);
                tmp_2 += 4 * vl;
                tmp += 4 * vl;
                i -= 4 * vl;
            }
            else
            {   
                vl = vsetvl_e32m2(i);
                v20 = vle32_v_f32m2( tmp, vl);
                vtmp20 = vfmul_vf_f32m2(v20,  tmp[0], vl);
                vse32_v_f32m2( tmp_2, vtmp20, vl);
                tmp_2 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m2 = chrono::high_resolution_clock::now();
    auto duration_e32m2 = chrono::duration_cast<chrono::microseconds>(stop_e32m2 - start_e32m2);
    cout << "e32m2 " << duration_e32m2.count()/num << endl;
    delete []e32m2;


    sleep(1);
    e32m1 = new _Float32[n];
    _Float32 *tmp_1 = e32m1;
    auto start_e32m1 = chrono::high_resolution_clock::now();
    vtmp10 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp11 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp12 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp13 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp14 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp15 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp16 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vtmp17 = vfmv_v_f_f32m1 (0, vsetvlmax_e32m1());
    vl_max = vsetvlmax_e32m1();
    for (int j = 0; j < num; j++)
    {   
        tmp_1 = e32m1;
        tmp =  &buffer[0];
        for (int i = n; i > 0;)
        {   
            if (i - 8 * vl_max > 0)
            {
                v10 = vle32_v_f32m1( tmp, vl);
                v11 = vle32_v_f32m1( (tmp + vl), vl);
                v12 = vle32_v_f32m1( (tmp + 2 * vl), vl);
                v13 = vle32_v_f32m1( (tmp + 3 * vl), vl);
                v14 = vle32_v_f32m1( (tmp + 4 * vl), vl);
                v15 = vle32_v_f32m1( (tmp + 5 * vl), vl);
                v16 = vle32_v_f32m1( (tmp + 6 * vl), vl);
                v17 = vle32_v_f32m1( (tmp + 7 * vl), vl);
                vtmp10 = vfmul_vf_f32m1(v10,  tmp[0], vl);
                vtmp11 = vfmul_vf_f32m1(v11,  tmp[0], vl);
                vtmp12 = vfmul_vf_f32m1(v12,  tmp[0], vl);
                vtmp13 = vfmul_vf_f32m1(v13,  tmp[0], vl);
                vtmp14 = vfmul_vf_f32m1(v14,  tmp[0], vl);
                vtmp15 = vfmul_vf_f32m1(v15,  tmp[0], vl);
                vtmp16 = vfmul_vf_f32m1(v16,  tmp[0], vl);
                vtmp17 = vfmul_vf_f32m1(v17,  tmp[0], vl);
                vse32_v_f32m1( tmp_1, vtmp10, vl);
                vse32_v_f32m1( (tmp_1 + vl), vtmp11, vl);
                vse32_v_f32m1( (tmp_1 + 2 * vl), vtmp12, vl);
                vse32_v_f32m1( (tmp_1 + 3 * vl), vtmp13, vl);
                vse32_v_f32m1( (tmp_1 + 4 * vl), vtmp14, vl);
                vse32_v_f32m1( (tmp_1 + 5 * vl), vtmp15, vl);
                vse32_v_f32m1( (tmp_1 + 6 * vl), vtmp16, vl);
                vse32_v_f32m1( (tmp_1 + 7 * vl), vtmp17, vl);
                tmp_1 += 8 * vl;
                tmp += 8 * vl;
                i -= 8 * vl;
            }
            else
            {   
                vl = vsetvl_e32m1(i);
                v10 = vle32_v_f32m1( tmp, vl);
                vtmp10 = vfmul_vf_f32m1(v10,  tmp[0], vl);
                vse32_v_f32m1(tmp_1, v10, vl);
                tmp_1 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m1 = chrono::high_resolution_clock::now();
    auto duration_e32m1 = chrono::duration_cast<chrono::microseconds>(stop_e32m1 - start_e32m1);
    cout << "e32m1 " << duration_e32m1.count()/num << endl;
    delete []e32m1;

    

    return 0;
}