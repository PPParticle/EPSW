// c906 长指令和短指令组合加载数据和存储数据的测试情况
// 数据量 <= 3*L1 cache 长指令速度>短指令组合
// 3*L1 cache < 数据量 e32m2速度>长指令 
#include <iostream>
#include <riscv-vector.h>
#include <assert.h>
#include <random>
#include <chrono>
#include <fstream>
#include <unistd.h>

// Bytes
#define cache_v 32 * 1024 
#define num 128

using namespace std;
int main(int argc, char **argv)
{

    ifstream fp("test_vl.txt", std::ios::binary);
    fp.seekg(0, std::ios::end);
    std::streamsize size = fp.tellg(); // signed integer

    // 重置读取位置到文件开头
    fp.seekg(0, std::ios::beg);

    // 读取文件数据并存储到 std::vector 中
    char* buffer = new char[size];
    if (fp.read(buffer, size))
    {
        std::cout << "Read " << size << " bytes from file." << std::endl;
    }
    else
    {
        std::cerr << "Error reading file." << std::endl;
    }
    fp.close(); // 关闭文件
    
    size_t vl;
    int32xm8_t v8, vtmp;
    int32xm4_t v40, v41, vtmp40, vtmp41;
    int32xm2_t v20, v21, v22, v23, vtmp20, vtmp21, vtmp22, vtmp23;
    int32xm1_t v10, v11, v12, v13, v14, v15, v16, v17, 
                vtmp10, vtmp11, vtmp12, vtmp13, vtmp14, vtmp15, vtmp16, vtmp17;
    int8_t *e32m8, *e32m4, *e32m2, *e32m1;

    size_t n = cache_v * atoi(argv[1]);
    cout << n << " bytes data" << endl;
    // struct mallinfo mi = mallinfo();
    // cout << "Maximum space allocated from system: " << mi.hblkhd << " bytes" << endl;

    e32m8 = new int8_t[n];
    int8_t *tmp_8 = e32m8;
    int8_t *tmp = (int8_t*)&buffer[0];
    auto start_e32m8 = chrono::high_resolution_clock::now();
    vtmp = vmvvx_int32xm8 (0.0, vsetvli_max(RVV_E32, RVV_M8));
    for (int i = 0; i < num; i++)
    {   
        tmp_8 = e32m8;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0; i -= vl)
        {
            vl = vsetvli(i, RVV_E32, RVV_M8);
            v8 = vlev_int32xm8((int32_t*)tmp, vl);
            vtmp = vmulvx_int32xm8(v8, (int32_t)buffer[0], vl);
            vsev_int32xm8((int32_t*)tmp_8, vtmp, vl);
            tmp_8 += vl;
            tmp += vl;
        }
    }
    
    auto stop_e32m8 = chrono::high_resolution_clock::now();
    auto duration_e32m8 = chrono::duration_cast<chrono::microseconds>(stop_e32m8 - start_e32m8);
    cout << "e32m8 " << duration_e32m8.count()/num << endl;

    sleep(1);
    e32m4 = new int8_t[n];
    int8_t *tmp_4 = e32m4;
    auto start_e32m4 = chrono::high_resolution_clock::now();
    vtmp40 = vmvvx_int32xm4 (0.0, vsetvli_max(RVV_E32, RVV_M4));
    vtmp41 = vmvvx_int32xm4 (0.0, vsetvli_max(RVV_E32, RVV_M4));
    for (int i = 0; i < num; i++)
    {   
        tmp_4 = e32m4;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0;)
        {   
            vl = vsetvli(i, RVV_E32, RVV_M4);
            if (i - 2 * vl > 0)
            {
                v40 = vlev_int32xm4((int32_t*)tmp, vl);
                v41 = vlev_int32xm4((int32_t*)(tmp + vl), vl);
                vtmp40 = vmulvx_int32xm4(v40, (int32_t)buffer[0], vl);
                vtmp41 = vmulvx_int32xm4(v41, (int32_t)buffer[0], vl);
                vsev_int32xm4((int32_t*)tmp_4, vtmp40, vl);
                vsev_int32xm4((int32_t*)(tmp_4 + vl), vtmp41, vl);
                tmp_4 += 2 * vl;
                tmp += 2 * vl;
                i -= 2 * vl;
            }
            else
            {
                v40 = vlev_int32xm4((int32_t*)tmp, vl);
                vtmp40 = vmulvx_int32xm4(v40, (int32_t)buffer[0], vl);
                vsev_int32xm4((int32_t*)tmp_4, vtmp40, vl);
                tmp_4 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m4 = chrono::high_resolution_clock::now();
    auto duration_e32m4 = chrono::duration_cast<chrono::microseconds>(stop_e32m4 - start_e32m4);
    cout << "e32m4 " << duration_e32m4.count()/num << endl;

    sleep(1);
    e32m2 = new int8_t[n];
    int8_t *tmp_2 = e32m2;
    auto start_e32m2 = chrono::high_resolution_clock::now();
    vtmp20 = vmvvx_int32xm2 (0.0, vsetvli_max(RVV_E32, RVV_M2));
    vtmp21 = vmvvx_int32xm2 (0.0, vsetvli_max(RVV_E32, RVV_M2));
    vtmp22 = vmvvx_int32xm2 (0.0, vsetvli_max(RVV_E32, RVV_M2));
    vtmp23 = vmvvx_int32xm2 (0.0, vsetvli_max(RVV_E32, RVV_M2));
    for (int i = 0; i < num; i++)
    {
        tmp_2 = e32m2;
        tmp = (int8_t*)&buffer[0];

        for (int i = n/4; i > 0;)
        {   
            vl = vsetvli(i, RVV_E32, RVV_M2);
            if (i - 4 * vl > 0)
            {
                v20 = vlev_int32xm2((int32_t*)tmp, vl);
                v21 = vlev_int32xm2((int32_t*)(tmp + vl), vl);
                v22 = vlev_int32xm2((int32_t*)(tmp + 2 * vl), vl);
                v23 = vlev_int32xm2((int32_t*)(tmp + 3 * vl), vl);
                vtmp20 = vmulvx_int32xm2(v20, (int32_t)buffer[0], vl);
                vtmp21 = vmulvx_int32xm2(v21, (int32_t)buffer[0], vl);
                vtmp22 = vmulvx_int32xm2(v22, (int32_t)buffer[0], vl);
                vtmp23 = vmulvx_int32xm2(v23, (int32_t)buffer[0], vl);
                vsev_int32xm2((int32_t*)tmp_2, vtmp20, vl);
                vsev_int32xm2((int32_t*)(tmp_2 + vl), vtmp21, vl);
                vsev_int32xm2((int32_t*)(tmp_2 + 2 * vl), vtmp22, vl);
                vsev_int32xm2((int32_t*)(tmp_2 + 3 * vl), vtmp23, vl);
                tmp_2 += 4 * vl;
                tmp += 4 * vl;
                i -= 4 * vl;
            }
            else
            {
                v20 = vlev_int32xm2((int32_t*)tmp, vl);
                vtmp20 = vmulvx_int32xm2(v20, (int32_t)buffer[0], vl);
                vsev_int32xm2((int32_t*)tmp_2, vtmp20, vl);
                tmp_2 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m2 = chrono::high_resolution_clock::now();
    auto duration_e32m2 = chrono::duration_cast<chrono::microseconds>(stop_e32m2 - start_e32m2);
    cout << "e32m2 " << duration_e32m2.count()/num << endl;

    sleep(1);
    e32m1 = new int8_t[n];
    int8_t *tmp_1 = e32m1;
    auto start_e32m1 = chrono::high_resolution_clock::now();
    vtmp10 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp11 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp12 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp13 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp14 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp15 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp16 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    vtmp17 = vmvvx_int32xm1 (0.0, vsetvli_max(RVV_E32, RVV_M1));
    for (int i = 0; i < num; i++)
    {   
        tmp_1 = e32m1;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0;)
        {
            vl = vsetvli(i, RVV_E32, RVV_M1);
            if (i - 8 * vl > 0)
            {
                v10 = vlev_int32xm1((int32_t*)tmp, vl);
                v11 = vlev_int32xm1((int32_t*)(tmp + vl), vl);
                v12 = vlev_int32xm1((int32_t*)(tmp + 2 * vl), vl);
                v13 = vlev_int32xm1((int32_t*)(tmp + 3 * vl), vl);
                v14 = vlev_int32xm1((int32_t*)(tmp + 4 * vl), vl);
                v15 = vlev_int32xm1((int32_t*)(tmp + 5 * vl), vl);
                v16 = vlev_int32xm1((int32_t*)(tmp + 6 * vl), vl);
                v17 = vlev_int32xm1((int32_t*)(tmp + 7 * vl), vl);
                vtmp10 = vmulvx_int32xm1(v10, (int32_t)buffer[0], vl);
                vtmp11 = vmulvx_int32xm1(v11, (int32_t)buffer[0], vl);
                vtmp12 = vmulvx_int32xm1(v12, (int32_t)buffer[0], vl);
                vtmp13 = vmulvx_int32xm1(v13, (int32_t)buffer[0], vl);
                vtmp14 = vmulvx_int32xm1(v14, (int32_t)buffer[0], vl);
                vtmp15 = vmulvx_int32xm1(v15, (int32_t)buffer[0], vl);
                vtmp16 = vmulvx_int32xm1(v16, (int32_t)buffer[0], vl);
                vtmp17 = vmulvx_int32xm1(v17, (int32_t)buffer[0], vl);
                vsev_int32xm1((int32_t*)tmp_1, vtmp10, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + vl), vtmp11, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 2 * vl), vtmp12, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 3 * vl), vtmp13, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 4 * vl), vtmp14, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 5 * vl), vtmp15, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 6 * vl), vtmp16, vl);
                vsev_int32xm1((int32_t*)(tmp_1 + 7 * vl), vtmp17, vl);
                tmp_1 += 8 * vl;
                tmp += 8 * vl;
                i -= 8 * vl;
            }
            else
            {
                v10 = vlev_int32xm1((int32_t*)tmp, vl);
                vtmp10 = vmulvx_int32xm1(v16, (int32_t)buffer[0], vl);
                vsev_int32xm1((int32_t*)tmp_1, vtmp10, vl);
                tmp_1 += vl;
                tmp += vl;
                i -= vl;
            }
        }
    }
    auto stop_e32m1 = chrono::high_resolution_clock::now();
    auto duration_e32m1 = chrono::duration_cast<chrono::microseconds>(stop_e32m1 - start_e32m1);
    cout << "e32m1 " << duration_e32m1.count()/num << endl;
    

    return 0;
}