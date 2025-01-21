#include<iostream>
#include<riscv_vector.h>
#include <sys/time.h>
using namespace std;

float *m1, *m2, *m_i, *m_r1, *m_r2, *m_o;
int main(int argc, char* argv[]){
    int n = atoi(argv[1]);
    m1 = new float[n*n];
    m2 = new float[n*n];
    m_r1 = new float[n*n];
    m_r2 = new float[n*n];
    m_i = new float[n*n];
    m_o = new float[n*n];

    for(int i=0; i<n*n; i++){
        m1[i] = float(rand() % 1048576) / 1048576;
        m2[i] = float(rand() % 1048576) / 1048576;
    }

    struct timeval starti, endi, starto, endo, startr1, endr1, startr2, endr2;

#define ITER 10
//inner-product
	gettimeofday(&starti, NULL);
	for(int loop=0;loop<ITER;loop++) {
        size_t vl;
        vfloat32m8_t vtmp, vb, va, vc;
        vfloat32m1_t vans;
        vans = vfmv_v_f_f32m1(0.0, vsetvlmax_e32m8());
		vtmp = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                vans = vfmv_v_f_f32m1(0.0, vsetvlmax_e32m8());
		        vtmp = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
                for(int k=0; k<n; k+=vl){
                    vl = vsetvl_e32m8(n-k);
                    va = vle32_v_f32m8(m1+i*n+k, vl);
                    vb = vlse32_v_f32m8(m2+j+k*n, n*4, vl);
                    vc = vfmul_vv_f32m8(va, vb, vl);
                    vtmp = vfadd_vv_f32m8(vtmp, vc, vl);
                }
                vans = vfredosum_vs_f32m8_f32m1(vans, vtmp, vans, vsetvlmax_e32m8());
                vse32_v_f32m1(m_i+i*n+j, vans, 1);
            }
        }
    }
    gettimeofday(&endi, NULL);

// outer-product
    gettimeofday(&starto, NULL);
	for(int loop=0;loop<ITER;loop++) {
        size_t vl;
        vfloat32m8_t vans, vtmp, vb, va, vc;
        vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
        for(int j=0; j<n; j++){
            for(int i=0; i<n; i++){
                float brd = m1[i*n+j];
                vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
                for(int k=0; k<n; k+=vl){
                    vl = vsetvl_e32m8(n-k);
                    vb = vle32_v_f32m8(m2+j*n+k, vl);
                    vc = vfmul_vf_f32m8(vb, brd, vl);
                    vans = vle32_v_f32m8(m_o+i*n+k, vl);
                    vans = vfadd_vv_f32m8(vans, vc, vl);
                    vse32_v_f32m8(m_o+i*n+k, vans, vl);
                }
            }
        }
    }
    gettimeofday(&endo, NULL);

    gettimeofday(&startr1, NULL);
	for(int loop=0;loop<ITER;loop++) {
        size_t vl;
        vfloat32m8_t vans, vtmp, vb, va, vc;
        vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
		vtmp = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                float brd = m1[i*n+j];
                vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
                for(int k=0; k<n; k+=vl){
                    vl = vsetvl_e32m8(n-k);
                    vb = vle32_v_f32m8(m2+j*n+k, vl);
                    vc = vfmul_vf_f32m8(vb, brd, vl);
                    vans = vle32_v_f32m8(m_r1+i*n+k, vl);
                    vans = vfadd_vv_f32m8(vans, vc, vl);
                    vse32_v_f32m8(m_r1+i*n+k, vans, vl);
                }
            }
        }

    }
    gettimeofday(&endr1, NULL);

    gettimeofday(&startr2, NULL);
	for(int loop=0;loop<ITER;loop++) {
        size_t vl;
        vfloat32m8_t vans, vtmp, vb, va, vc;
        vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
        for(int i=0; i<n; i++){
            for(int k=0; k<n; k+=vl){
                vans = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
                vl = vsetvl_e32m8(n-k);
                for(int j=0; j<n; j++){
                    float brd = m1[i*n+j];
                    vb = vle32_v_f32m8(m2+j*n+k, vl);
                    vc = vfmul_vf_f32m8(vb, brd, vl);
                    vans = vfadd_vv_f32m8(vans, vc, vl);
                }
                vse32_v_f32m8(m_r2+i*n+k, vans, vl);
            }
        }
    }
    gettimeofday(&endr2, NULL);

    double elapsed[4];
    elapsed[0] = ((endi.tv_sec-starti.tv_sec)*1000000 + endi.tv_usec-starti.tv_usec);
    elapsed[1] = ((endo.tv_sec-starto.tv_sec)*1000000 + endo.tv_usec-starto.tv_usec);
    elapsed[2] = ((endr1.tv_sec-startr1.tv_sec)*1000000 + endr1.tv_usec-startr1.tv_usec);
    elapsed[3] = ((endr2.tv_sec-startr2.tv_sec)*1000000 + endr2.tv_usec-startr2.tv_usec);

    #define VALIDATE
#if defined VALIDATE
	for (int i = 0; i < n*n; i++)
	{
		float p0 = m_i[i]/10;
		float p1 = m_o[i]/10;
		float p2 = m_r1[i]/10;
        float p3 = m_r2[i]/10;


		if (p0 < 0)
			p0 *= -1;
		if (p1 < 0)
			p1 *= -1;
		if (p2 < 0)
			p2 *= -1;
		if (p0!=p1 || p1!=p2 || p2!=p0 || p0!=p3 || p1!=p3 || p2!=p3)
		{
		cout << "("<<i<<", "<<p0<<", "<<p1<<", "<<p2 <<", "<<p3 <<")";
		cout << endl;
		return -1;
		}
	}
#endif
	fprintf(stdout, "%f, %f, %f, %f\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, elapsed[2]/(double)ITER, elapsed[3]/(double)ITER);
}