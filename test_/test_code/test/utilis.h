#ifndef _UTILS_H_
#define _UTILS_H_

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

void print(float *a, unsigned int ai, unsigned int aj)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < aj; j++)
        {
            // cout << a[i*aj+j] << " ";
            printf("%8.2lf ", a[i * aj + j]);
        }
        cout << endl;
    }
}
void print(int *a, unsigned int ai, unsigned int aj)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < aj; j++)
        {
            cout << a[i * aj + j] << " ";
        }
        cout << endl;
    }
}

void print(unsigned int *a, unsigned int ai, unsigned int aj)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < aj; j++)
        {
            cout << a[i * aj + j] << " ";
        }
        cout << endl;
    }
}

void print(float **a, unsigned int ai, int *size)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < size[i]; j++)
        {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}
void print(int **a, unsigned int ai, int *size)
{
    for (int i = 0; i < ai; i++)
    {
        for (int j = 0; j < size[i]; j++)
        {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }
}

void print(float *a, unsigned int ai)
{
    for (int i = 0; i < ai; i++)
    {

        cout << a[i] << " ";
    }
    cout << endl;
}

void print(int *a, unsigned int ai)
{
    for (int i = 0; i < ai; i++)
    {

        cout << a[i] << " ";
    }
    cout << endl;
}

#endif
