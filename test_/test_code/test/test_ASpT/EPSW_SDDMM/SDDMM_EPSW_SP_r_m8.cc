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
	sc = atoi(argv[2]);
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
	fprintf(stdout, "%d %d %d %d ", ne, nr, nc, sc);
}

void gen()
{
	M.ELLPACKinitWithCSR(csr_v, csr_e, csr_ev, nr, nc, ne);
    M.ELLPACK_P_vreorder();
}

void mprocess()
{
	double elapsed[3];
	FTYPE *vin, *vout, *vout_scalar;
	vin = new FTYPE[nc * sc]{0};
	vout = new FTYPE[nr * sc]{0};
	vout_scalar = new FTYPE[ne]{0};
	
	FTYPE *vout_EPSW = new FTYPE[nr * nc]{0.0};
	// vout_EPSW = new FTYPE[ne]{0};

	for (int i = 0; i < nc * sc; i++)
		vin[i] = (FTYPE)(rand() & 1048575) / 1048576;
		// vin[i] = i;

	for (int i = 0; i < nr * sc; i++)
		vout[i] = (FTYPE)(rand() & 1048575) / 1048576;
		// vout[i] = 1;

	struct timeval starttime_scalar, endtime_scalar, starttime_EPSW, endtime_EPSW;
	#define ITER 10
	////begin
	gettimeofday(&starttime_EPSW, NULL);
	for(int loop=0;loop<ITER;loop++) {
			size_t vl;
			float32xm8_t va, vc, vo;
			int32xm8_t vidx, vidx_t, vidx_s;

		 for (int i = 0; i < M._rows; i++){   
			int num = M._nnz_num[i];
			int m = M._rowidx[i];

			for(int j=0; j<num; j+=vl){
				vl = vsetvli(num-j, RVV_E32, RVV_M8);
				va = vlev_float32xm8(M._rowPtr[i]+j, vl);
				vidx = vlev_int32xm8((int *)M._coordPtr[i]+j, vl);
				vidx_s = vsllvx_int32xm8(vidx, 2, vl);
				vo = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
				vo = vlxev_float32xm8(vout_EPSW+m*nc, vidx_s, vl);
				for(int k=0; k<sc; k++){
					vidx_t = vmulvx_int32xm8(vidx, sc, vl);
					vidx_t = vaddvx_int32xm8(vidx_t, k, vl);
					vidx_t = vsllvx_int32xm8(vidx_t, 2, vl);
					vc = vlxev_float32xm8(vin, vidx_t, vl);
					vo = vfmaccvf_float32xm8(vo, vout[m*sc+k], vc, vl);
				}
				va = vfmulvv_float32xm8(va, vo, vl);
				vsxev_float32xm8(vout_EPSW+m*nc, vidx_s, va, vl);

			}
		}
	}
	gettimeofday(&endtime_EPSW, NULL);

	gettimeofday(&starttime_scalar, NULL);
		for(int loop=0;loop<ITER;loop++) {
        	for(int i=0;i<nr;i++){
                for(int j=csr_v[i]; j<csr_v[i+1]; j++){
                    for(int k=0; k<sc; k++){
                        vout_scalar[j] += vin[csr_e[j]*sc + k] * vout[i*sc + k];
					}
					vout_scalar[j] *= csr_ev[j];
            	}
        	}
		}
	gettimeofday(&endtime_scalar, NULL);
	
	elapsed[0] = ((endtime_EPSW.tv_sec-starttime_EPSW.tv_sec)*1000000 + endtime_EPSW.tv_usec-starttime_EPSW.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_scalar.tv_sec-starttime_scalar.tv_sec)*1000000 + endtime_scalar.tv_usec-starttime_scalar.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	
	FTYPE *csr_eo = new FTYPE[ne]{0.0};
	int init = 0;
	for(int i=0; i<nr; i++){
		for(int j=csr_v[i]; j<csr_v[i+1]; j++)
			csr_eo[init++] = vout_EPSW[i*nc+csr_e[j]];		
	}

	

	for (int i = 0; i < ne; i++){	

			FTYPE p1 = csr_eo[i];
			FTYPE p2 = vout_scalar[i];

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

	fprintf(stdout, "%f|%f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
