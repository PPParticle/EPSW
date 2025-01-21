#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <riscv_vector.h>
#include <iostream>
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float

int sc, nr, nc, ne;
int *csr_v;
int *csr_e;
FTYPE *csr_ev;
char* dir;
char *app;

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;

	app = argv[0];
	sc = atoi(argv[2]);
	dir = argv[1];
	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);
	// printf("%d, %d, %d", nr, nc, ne);

	csr_v = new int[nr + 1]{0};	  
	csr_e = new int[ne]{0};	   
	csr_ev = new FTYPE[ne]{0.0}; 

	// csr 赋值
	for (i = 0; i < nr+1; i++)
	{
		fscanf(fp, "%d ", &csr_v[i]);
	}
	// for(;i<nr+1;i++){
	// 	csr_v[i] = csr_v[i-1];
	// }
	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e[i]);
		csr_ev[i] = (FTYPE)(rand() % 1048576) / 1048576;
	}
	fclose(fp);
	fprintf(stdout, "%d %d %d %d ", ne, nr, nc, sc);
}

void mprocess()
{
	double elapsed[3];
	FTYPE *vin, *vout, *vout_rvv, *vout_scalar;
	vin = new FTYPE[nc*sc]{0.0};
	vout = new FTYPE[nr*sc]{0.0};
	vout_scalar = new FTYPE[ne]{0.0};
	vout_rvv = new FTYPE[ne]{0.0};

	struct timeval starttime_scalar, endtime_scalar,
					starttime_rvv, endtime_rvv;

	for (int i = 0; i < nc * sc; i++)
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
	
	for (int i = 0; i < nr * sc; i++)
		vout[i] = (FTYPE)(rand() % 1048576) / 1048576;
	
#define ITER 10
#ifdef ITER
	// printf("begin normal\n");
	gettimeofday(&starttime_scalar, NULL);
	for(int loop=0;loop<ITER;loop++) {
        for(int i=0;i<nr;i++){
			for(int j=csr_v[i]; j<csr_v[i+1]; j++){
        		for(int k=0; k<sc; k++)
                    vout_scalar[j] += vin[csr_e[j]*sc + k] * vout[i*sc + k];
				vout_scalar[j] *= csr_ev[j];
            }
        }
	}
	gettimeofday(&endtime_scalar, NULL);

	// 	// o*a*b = c
	gettimeofday(&starttime_rvv, NULL);
	for (int loop = 0; loop < ITER; loop++){

		for(int i=0; i<nr; i++){

			int loc1 = csr_v[i], loc2 = csr_v[i+1];
			int interm = loc1 + (((loc2-loc1)>>3)<<3);
			long j;

			size_t vl=vsetvlmax_e32m2();
			vfloat32m2_t va, vb;
			vuint32m2_t vidx, vidx_t;
			vb = vfmv_v_f_f32m2(0.0, vl);
			for(j=loc1; j<interm; j+=8){
				vidx = vle32_v_u32m2((unsigned int *)(csr_e+j), vl);
				vidx = vmul_vx_u32m2(vidx, sc, vl);
				for(int k=0; k<sc; k++){
					vb = vle32_v_f32m2(vout_rvv+j, vl);
					vidx_t = vadd_vx_u32m2(vidx, k, vl);
					vidx_t = vsll_vx_u32m2(vidx_t, 2, vl);
					va = vloxei32_v_f32m2(vin, vidx_t, vl);
					vb = vfmacc_vf_f32m2(vb, vout[i*sc+k], va, vl);
					vse32_v_f32m2(vout_rvv+j, vb, vl);
				}
				for(int k=0; k<8; k++) 
					vout_rvv[j+k] *= csr_ev[j+k];
			}
			for(; j<loc2; j++){
				for(int k=0; k<sc; k++) 
               		vout_rvv[j] += vin[csr_e[j]*sc + k] * vout[i*sc + k];
				vout_rvv[j] *= csr_ev[j];
			}
			
		}

		
	}
#endif
	gettimeofday(&endtime_rvv, NULL);

elapsed[0]  = ((endtime_rvv.tv_sec-starttime_rvv.tv_sec)*1000000 + endtime_rvv.tv_usec-starttime_rvv.tv_usec)/1000000.0;
elapsed[1]  = ((endtime_scalar.tv_sec-starttime_scalar.tv_sec)*1000000 + endtime_scalar.tv_usec-starttime_scalar.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	int num_diff = 0;
	for (int i = 0; i < ne; i++)
	{
		FTYPE p1 = vout_scalar[i];
		FTYPE p2 = vout_rvv[i];

		if (p1 < 0)
			p1 *= -1;
		if (p2 < 0)
			p2 *= -1;
		FTYPE diff;
		diff = p1 - p2;
		if (diff < 0)
			diff *= -1;
		if (MAX(p1, p2) !=0 && diff / MAX(p1, p2) > 0.01){	
			fprintf(stdout, "%f %f %s %s\n", p1, p2, app, dir);
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
