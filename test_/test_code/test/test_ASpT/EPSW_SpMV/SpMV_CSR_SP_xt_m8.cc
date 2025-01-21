#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <riscv_vector.h>
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float

int sc, nr, nc, ne;
int *csr_v;
int *csr_e;
FTYPE *csr_ev;
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
	for(int loop=0;loop<ITER;loop++) {

		for (int i = 0; i < nr; i++){
			size_t vl;
			vfloat32m8_t va, vb, vtmp;
			vfloat32m1_t vans;
			vuint32m8_t vcidx;
			vans = vfmv_v_f_f32m1(0.0, 1);
			vtmp = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
			for (int j = csr_v[i]; j < csr_v[i+1]; j+=vl){
				vl = vsetvl_e32m8(csr_v[i+1] - j);
				va = vle32_v_f32m8(&csr_ev[j], vl);
				vcidx = vle32_v_u32m8((unsigned int*)&csr_e[j], vl);
				vcidx = vsll_vx_u32m8(vcidx, 2, vl);
				vb = vloxei32_v_f32m8(vin, vcidx, vl);
				vtmp = vfmacc_vv_f32m8(vtmp, va, vb, vsetvlmax_e32m8());
			}
			vans = vfredusum_vs_f32m8_f32m1(vans, vtmp, vans, vsetvlmax_e32m8());
			vout_EPSW[i] += vfmv_f_s_f32m1_f32(vans);		
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
	mprocess();
}
