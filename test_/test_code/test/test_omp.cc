#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <chrono>

#define N 1000

int main() {
    int sum = 0;
    auto start = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < N; i++) {
        sum += i;
    }
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    int sum1 = 0;
    auto start1 = std::chrono::high_resolution_clock::now();
    for (int i=0; i < N; i++)
        sum1 +=i;
    auto end1 = std::chrono::high_resolution_clock::now();
    double duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start1).count();
    printf("Sum: %d %d\n", sum, sum1);
    std::cout << "sum duration " << duration << std::endl;
    std::cout << "sum1 duration " << duration1 << std::endl;

    return 0;
}
