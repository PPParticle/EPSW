// c910 长指令和短指令组合加载数据和存储数据的测试情况
// 总是 长指令速度>短指令组合 
#include <iostream>
#include <riscv_vector.h>
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
    char *buffer = new char [size];
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
    vint32m8_t v8, vtmp;
    vint32m4_t v40, v41, vtmp40, vtmp41;
    vint32m2_t v20, v21, v22, v23, vtmp20, vtmp21, vtmp22, vtmp23;
    vint32m1_t v10, v11, v12, v13, v14, v15, v16, v17, 
                vtmp10, vtmp11, vtmp12, vtmp13, vtmp14, vtmp15, vtmp16, vtmp17;
    int8_t *e32m8, *e32m4, *e32m2, *e32m1;

    size_t n = cache_v * atoi(argv[1]);
    cout << n << " bytes data" << endl;

    e32m8 = new int8_t[n];
    int8_t *tmp_8 = e32m8;
    int8_t *tmp = (int8_t*)&buffer[0];
    auto start_e32m8 = chrono::high_resolution_clock::now();
    vtmp = vmv_v_x_i32m8 (0, vsetvlmax_e32m8());
    for (int i = 0; i < num; i++)
    {   
        tmp_8 = e32m8;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0; i -= vl)
        {
            vl = vsetvl_e32m8(i);
            v8 = vle32_v_i32m8((int32_t*)tmp, vl);
            vtmp = vmul_vx_i32m8(v8, (int32_t)buffer[0], vl);
            vse32_v_i32m8((int32_t*)tmp_8, vtmp, vl);
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
    vtmp40 = vmv_v_x_i32m4 (0, vsetvlmax_e32m4());
    vtmp41 = vmv_v_x_i32m4 (0, vsetvlmax_e32m4());
    for (int i = 0; i < num; i++)
    {   
        tmp_4 = e32m4;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0;)
        {   
            vl = vsetvl_e32m4(i);
            if (i - 2 * vl > 0)
            {   
                v40 = vle32_v_i32m4((int32_t*)tmp, vl);
                v41 = vle32_v_i32m4((int32_t*)(tmp + vl), vl);
                vtmp40 = vmul_vx_i32m4(v40, (int32_t)buffer[0], vl);
                vtmp41 = vmul_vx_i32m4(v41, (int32_t)buffer[0], vl);
                vse32_v_i32m4((int32_t*)tmp_4, vtmp40, vl);
                vse32_v_i32m4((int32_t*)(tmp_4 + vl), vtmp41, vl);
                tmp_4 += 2 * vl;
                tmp += 2 * vl;
                i -= 2 * vl;
            }
            else
            {   
                v40 = vle32_v_i32m4((int32_t*)tmp, vl);
                vtmp40 = vmul_vx_i32m4(v40, (int32_t)buffer[0], vl);
                vse32_v_i32m4((int32_t*)tmp_4, vtmp40, vl);
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
    vtmp20 = vmv_v_x_i32m2 (0, vsetvlmax_e32m2());
    vtmp21 = vmv_v_x_i32m2 (0, vsetvlmax_e32m2());
    vtmp22 = vmv_v_x_i32m2 (0, vsetvlmax_e32m2());
    vtmp23 = vmv_v_x_i32m2 (0, vsetvlmax_e32m2());
    for (int i = 0; i < num; i++)
    {   
        tmp_2 = e32m2;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0;)
        {   
            vl = vsetvl_e32m2(i);
            if (i - 4 * vl > 0)
            {   
                v20 = vle32_v_i32m2((int32_t*)tmp, vl);
                v21 = vle32_v_i32m2((int32_t*)(tmp + vl), vl);
                v22 = vle32_v_i32m2((int32_t*)(tmp + 2 * vl), vl);
                v23 = vle32_v_i32m2((int32_t*)(tmp + 3 * vl), vl);
                vtmp20 = vmul_vx_i32m2(v20, (int32_t)buffer[0], vl);
                vtmp21 = vmul_vx_i32m2(v21, (int32_t)buffer[0], vl);
                vtmp22 = vmul_vx_i32m2(v22, (int32_t)buffer[0], vl);
                vtmp23 = vmul_vx_i32m2(v23, (int32_t)buffer[0], vl);
                vse32_v_i32m2((int32_t*)tmp_2, vtmp20, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + vl), vtmp21, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + 2 * vl), vtmp22, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + 3 * vl), vtmp23, vl);
                tmp_2 += 4 * vl;
                tmp += 4 * vl;
                i -= 4 * vl;
            }
            else
            {
                v20 = vle32_v_i32m2((int32_t*)tmp, vl);
                vtmp20 = vmul_vx_i32m2(v20, (int32_t)buffer[0], vl);
                vse32_v_i32m2((int32_t*)tmp_2, vtmp20, vl);
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
    vtmp10 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp11 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp12 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp13 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp14 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp15 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp16 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    vtmp17 = vmv_v_x_i32m1 (0, vsetvlmax_e32m1());
    for (int i = 0; i < num; i++)
    {   
        tmp_1 = e32m1;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0;)
        {   
            vl = vsetvl_e32m1(i);
            if (i - 8 * vl > 0)
            {
                v10 = vle32_v_i32m1((int32_t*)tmp, vl);
                v11 = vle32_v_i32m1((int32_t*)(tmp + vl), vl);
                v12 = vle32_v_i32m1((int32_t*)(tmp + 2 * vl), vl);
                v13 = vle32_v_i32m1((int32_t*)(tmp + 3 * vl), vl);
                v14 = vle32_v_i32m1((int32_t*)(tmp + 4 * vl), vl);
                v15 = vle32_v_i32m1((int32_t*)(tmp + 5 * vl), vl);
                v16 = vle32_v_i32m1((int32_t*)(tmp + 6 * vl), vl);
                v17 = vle32_v_i32m1((int32_t*)(tmp + 7 * vl), vl);
                vtmp10 = vmul_vx_i32m1(v10, (int32_t)buffer[0], vl);
                vtmp11 = vmul_vx_i32m1(v11, (int32_t)buffer[0], vl);
                vtmp12 = vmul_vx_i32m1(v12, (int32_t)buffer[0], vl);
                vtmp13 = vmul_vx_i32m1(v13, (int32_t)buffer[0], vl);
                vtmp14 = vmul_vx_i32m1(v14, (int32_t)buffer[0], vl);
                vtmp15 = vmul_vx_i32m1(v15, (int32_t)buffer[0], vl);
                vtmp16 = vmul_vx_i32m1(v16, (int32_t)buffer[0], vl);
                vtmp17 = vmul_vx_i32m1(v17, (int32_t)buffer[0], vl);
                vse32_v_i32m1((int32_t*)tmp_1, vtmp10, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + vl), vtmp11, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 2 * vl), vtmp12, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 3 * vl), vtmp13, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 4 * vl), vtmp14, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 5 * vl), vtmp15, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 6 * vl), vtmp16, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 7 * vl), vtmp17, vl);
                tmp_1 += 8 * vl;
                tmp += 8 * vl;
                i -= 8 * vl;
            }
            else
            {
                v10 = vle32_v_i32m1((int32_t*)tmp, vl);
                vtmp10 = vmul_vx_i32m1(v10, (int32_t)buffer[0], vl);
                vse32_v_i32m1((int32_t*)tmp_1, v10, vl);
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