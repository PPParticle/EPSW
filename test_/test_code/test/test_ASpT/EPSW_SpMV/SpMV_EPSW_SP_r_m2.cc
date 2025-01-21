#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "ELLPACK_P_0.7.h"
#include <riscv-vector.h>
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float

int sc, nr, nc, ne;
int *csr_v;
int *csr_e;
FTYPE *csr_ev;
ELLPACK_P M;
char *dir;
char *app;

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;

	app = argv[0];
	dir = argv[1];
	fp = fopen(dir, "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);

	csr_v = new int[nr + 1]{0};	  
	csr_e = new int[ne]{0};	   
	csr_ev = new FTYPE[ne]{0.0}; 

	// csr 赋值
	for (i = 0; i < nr+1; i++){
		fscanf(fp, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne; i++){
		fscanf(fp, "%d", &csr_e[i]);
		csr_ev[i] = (FTYPE)(rand() & 1048575) / 1048576;
		// csr_ev[i] = 1;
	}
	fclose(fp);
	fprintf(stdout, "%d %d %d ", ne, nc, (ne+nc)*4/1024);
}

void gen()
{
	// spmv rvv test    
	M.ELLPACKinitWithCSR(csr_v, csr_e, csr_ev, nr, nc, ne);
    M.ELLPACK_P_vreorder();
}

void mprocess()
{
	double elapsed[3];
	FTYPE *vin = new FTYPE[nc]{0};
	FTYPE *vout_scalar = new FTYPE[nr]{0};	
	FTYPE *vout_EPSW = new FTYPE[nr]{0.0};

	for (int i = 0; i < nc; i++)
		vin[i] = (FTYPE)(rand() & 1048575) / 1048576;
		// vin[i] = i;

	struct timeval starttime_scalar, endtime_scalar, starttime_EPSW, endtime_EPSW;
	#define ITER 100
	////begin
	gettimeofday(&starttime_EPSW, NULL);
	// M.print_nnz();
	for(int loop=0;loop<ITER;loop++) {
		for (int i = 0; i < M._rows; i++){   
			int num = M._nnz_num[i];
			int ridx = M._rowidx[i];

			size_t vl;
			float32xm2_t va, vb, vtmp, vans;
			uint32xm2_t vcidx;
			vans = vfmvvf_float32xm2(0.0, 1);
			vtmp = vfmvvf_float32xm2(0.0, vsetvli_max(RVV_E32, RVV_M2));
			for(int j=0; j<num; j+=vl){
				vl = vsetvli(num-j, RVV_E32, RVV_M2);
				va = vlev_float32xm2(M._rowPtr[i]+j, vl);
				vcidx = vlev_uint32xm2(M._coordPtr[i]+j, vl);
				vcidx = vsllvx_uint32xm2(vcidx, 2, vl);
				vb = vlxev_float32xm2(vin, vcidx, vl);
				vtmp = vfmaccvv_float32xm2(vtmp, va, vb, vsetvli_max(RVV_E32, RVV_M2));
				// if(ridx==0){
				// float *test = new float[vsetvli_max(RVV_E32, RVV_M2)];
				// vsev_float32xm2(test, vtmp, vsetvli_max(RVV_E32, RVV_M2));
				// for(int i=0; i<vsetvli_max(RVV_E32, RVV_M2); i++){
				// 	cout << test[i] << " ";
				// }
				// cout << endl;
				// }
			}
			vans = vfredsumvs_float32xm2(vtmp, vans, vsetvli_max(RVV_E32, RVV_M2));
			// if(ridx==0){cout << vfmvfs_float32xm2(vans, 1) << endl;}
			vout_EPSW[ridx] += vfmvfs_float32xm2(vans, 1);		
		}
	}
	gettimeofday(&endtime_EPSW, NULL);

	gettimeofday(&starttime_scalar, NULL);
		for (int loop = 0; loop < ITER; loop++){
		for (int i = 0; i < nr; i++){
			for (int k = csr_v[i]; k < csr_v[i + 1]; k++){
				vout_scalar[i] += csr_ev[k] * vin[csr_e[k]];	
			}
		}
	}
	gettimeofday(&endtime_scalar, NULL);
	
	elapsed[0] = ((endtime_EPSW.tv_sec-starttime_EPSW.tv_sec)*1000000 + endtime_EPSW.tv_usec-starttime_EPSW.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_scalar.tv_sec-starttime_scalar.tv_sec)*1000000 + endtime_scalar.tv_usec-starttime_scalar.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	
	for (int i = 0; i < nr; i++){	

			FTYPE p1 = vout_EPSW[i]/ITER;
			FTYPE p2 = vout_scalar[i]/ITER;

			if (p1 < 0)
				p1 *= -1;
			if (p2 < 0)
				p2 *= -1;
			FTYPE diff;
			diff = p1 - p2;
			if (diff < 0)
				diff *= -1;
			if (MAX(p1,p2) !=0 && diff / MAX(p1, p2) > 0.01){
				fprintf(stdout, "%d %f %f %s %s\n", i, p1, p2, app, dir);
				return;
			}
	}
	
#endif

	fprintf(stdout, "%f %f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
