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
    vint32m8_t v8;
    vint32m4_t v40, v41;
    vint32m2_t v20, v21, v22, v23;
    vint32m1_t v10, v11, v12, v13, v14, v15, v16, v17;
    int8_t *e32m8, *e32m4, *e32m2, *e32m1;

    size_t n = cache_v * atoi(argv[1]);
    cout << n << " bytes data" << endl;

    e32m8 = new int8_t[n];
    int8_t *tmp_8 = e32m8;
    int8_t *tmp = (int8_t*)&buffer[0];
    auto start_e32m8 = chrono::high_resolution_clock::now();
    for (int i = 0; i < num; i++)
    {   
        tmp_8 = e32m8;
        tmp = (int8_t*)&buffer[0];
        for (int i = n/4; i > 0; i -= vl)
        {
            vl = vsetvl_e32m8(i);
            v8 = vle32_v_i32m8((int32_t*)tmp, vl);
            vse32_v_i32m8((int32_t*)tmp_8, v8, vl);
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
                vse32_v_i32m4((int32_t*)tmp_4, v40, vl);
                vse32_v_i32m4((int32_t*)(tmp_4 + vl), v41, vl);
                tmp_4 += 2 * vl;
                tmp += 2 * vl;
                i -= 2 * vl;
            }
            else
            {   
                v40 = vle32_v_i32m4((int32_t*)tmp, vl);
                vse32_v_i32m4((int32_t*)tmp_4, v40, vl);
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
                vse32_v_i32m2((int32_t*)tmp_2, v20, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + vl), v21, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + 2 * vl), v22, vl);
                vse32_v_i32m2((int32_t*)(tmp_2 + 3 * vl), v23, vl);
                tmp_2 += 4 * vl;
                tmp += 4 * vl;
                i -= 4 * vl;
            }
            else
            {
                v20 = vle32_v_i32m2((int32_t*)tmp, vl);
                vse32_v_i32m2((int32_t*)tmp_2, v20, vl);
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
                vse32_v_i32m1((int32_t*)tmp_1, v10, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + vl), v11, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 2 * vl), v12, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 3 * vl), v13, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 4 * vl), v14, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 5 * vl), v15, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 6 * vl), v16, vl);
                vse32_v_i32m1((int32_t*)(tmp_1 + 7 * vl), v17, vl);
                tmp_1 += 8 * vl;
                tmp += 8 * vl;
                i -= 8 * vl;
            }
            else
            {
                v10 = vle32_v_i32m1((int32_t*)tmp, vl);
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