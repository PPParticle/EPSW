#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <riscv-vector.h>
#include <iostream>
using namespace std;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float

int sc, nr, nc, ne;

int *csr_v, *csr_e;
FTYPE *csr_ev;

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;
	sc = atoi(argv[2]);
	//fprintf(stdout, "%s, \n", argv[1]);
	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);
	printf("%d,%d,%d,%d ", nr, nc, ne, sc);
	

	csr_v = new int[nr + 1]{0};	  // csr 行数组
	csr_e = new int[ne]{0};	  // csr 列索引
	csr_ev = new FTYPE[ne]{0.0}; // csr 元素数组
	// csr 赋值
	for(i=0; i<nr+1; i++){
		fscanf(fp, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e[i]);
		csr_ev[i] = (FTYPE)(rand() % 1048576) / 1048576;
	}
	fclose(fp);
}

void mprocess()
{
	double elapsed[2];
	FTYPE *vin, *vout, *vout_rvv;
	vin = new FTYPE[nc*sc]{0.0};
	vout = new FTYPE[nr*sc]{0.0};
	vout_rvv = new FTYPE[nr*sc]{0.0};

	struct timeval starttime_n, endtime_n,
					starttime_rvv, endtime_rvv;
	for (int i = 0; i < nc * sc; i++)
	{
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// vin[i] = i;
	}

#define ITER 10
#ifdef ITER
	gettimeofday(&starttime_n, NULL);
	for (int loop = 0; loop < ITER; loop++)
	{	
		for(int i=0; i<nr; i++){
			for(int j=0; j<sc; j++){
				for(int k=csr_v[i]; k<csr_v[i+1]; k++){
					vout[i*sc+j] += csr_ev[k] * vin[csr_e[k]*sc + j];
				}
			}
		}
	}
	gettimeofday(&endtime_n, NULL);

	gettimeofday(&starttime_rvv, NULL);
	// rvv
	for (int loop = 0; loop < ITER; loop++)
	{
		float32xm8_t va, vb, vtmp, vans;
		int32xm8_t vidx;
		size_t vl = vsetvli_max(RVV_E32,RVV_M2);
		for(int i=0; i<nr+1; i++){
			for(int j=0; j<sc; j++){
					int num = csr_v[i+1]-csr_v[i];
					if(num > 32){
						int itr = num/vl;
						int remain = csr_v[i] + vl * itr;
						for(int l=0; l<itr; l++){
							vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32,RVV_M8));
							vans = vfmvvf_float32xm8(0.0, vl);
							va = vlev_float32xm8(csr_ev+csr_v[i]+l*vl, vl);
							vidx = vlev_int32xm8(csr_e+csr_v[i]+l*vl, vl);
							vidx = vmulvx_int32xm8(vidx, sc, vl);
							vidx = vaddvx_int32xm8(vidx, j, vl);
							vidx = vsllvx_int32xm8(vidx, 2, vl);
							vb = vlxev_float32xm8(vin, vidx, vl);
							vtmp = vfmulvv_float32xm8(va, vb, vl);
							vans = vfredsumvs_float32xm8(vtmp, vans, vl);
							vout_rvv[i*sc+j] += vfmvfs_float32xm8(vans, 1);
						}
						// cout << remain << endl;
						for(int k=remain; k<csr_v[i+1]; k++){
							vout_rvv[i*sc+j] += csr_ev[k] * vin[csr_e[k]*sc + j];
						}
					}
					else{
						for(int k=csr_v[i]; k<csr_v[i+1]; k++){
							vout_rvv[i*sc+j] += csr_ev[k] * vin[csr_e[k]*sc + j];
						}
					}
					
				// }
			}
		}
	}
		gettimeofday(&endtime_rvv, NULL);
#endif

elapsed[0]  = ((endtime_n.tv_sec-starttime_n.tv_sec)*1000000 + endtime_n.tv_usec-starttime_n.tv_usec)/1000000.0;
elapsed[1]  = ((endtime_rvv.tv_sec-starttime_rvv.tv_sec)*1000000 + endtime_rvv.tv_usec-starttime_rvv.tv_usec)/1000000.0;
#define VALIDATE
#if defined VALIDATE
	for (int i = 0; i < nr * sc; i++)
	{
		FTYPE p1 = vout[i];
		FTYPE p2 = vout_rvv[i];

		if (p1 < 0)
			p1 *= -1;
		if (p2 < 0)
			p2 *= -1;
		FTYPE diff;
		diff = p1 - p2;
		if (diff < 0)
			diff *= -1;
		if (diff / MAX(p1, p2) > 0.01)
		{	
			cout << "("<<i<<","<<"."<<p1<<","<<p2<<") ";
			break;
		}

	}
#endif
	fprintf(stdout, "%f,%f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	mprocess();
}
