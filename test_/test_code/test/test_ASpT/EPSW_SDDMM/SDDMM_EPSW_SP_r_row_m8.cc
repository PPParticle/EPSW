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
	// fprintf(stdout, "%d %d %d %d ", ne, nr, nc, sc);
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
	#define EXCUTION 2
	////begin
	gettimeofday(&starttime_EPSW, NULL);

	for(int loop=0;loop<ITER;loop++) {
	#if EXCUTION == 2
		vector<vector<unsigned int>> stripes;
		stripes = M.stripe_adaptive();
		// o*a*b
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int strips = M._strips;
    	size_t vl;
		int row=0;
    	for(int i=0; i<stripes.size(); i++){
        	for(int j=0; j<M._nnz_num[row]; j++){
            	for(int k=0; k<stripes[i].size() && row+k<m; k++){
                	if(j<M._nnz_num[row+k]){
                    	Val          brd   = M._rowPtr[row+k][j];
                    	unsigned int brd_j = M._coordPtr[row+k][j];
						unsigned int brd_i = M._rowidx[row+k];
						Val 		 *tmp_a  = vout + brd_i * sc;
                    	Val          *tmp_b  = vin + brd_j*sc;

						float32xm8_t vans = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
						float32xm8_t vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<sc; l+=vl){
	                        vl=vsetvli(sc-l, RVV_E32, RVV_M8);;
							float32xm8_t va=vlev_float32xm8(tmp_a, vl);
                        	float32xm8_t vb=vlev_float32xm8(tmp_b, vl); 
							vtmp=vfmulvv_float32xm8(va, vb, vl);	
                        	vans = vfaddvv_float32xm8(vans, vtmp, vl);
					
                        	tmp_a += vl;
                        	tmp_b += vl;
                    	}
						
						float32xm8_t tmp = vfmvsf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
						tmp =  vfredsumvs_float32xm8(vans, tmp, vl);
						vout_EPSW[brd_i*nc + brd_j] +=  vfmvfs_float32xm8(tmp, 1);
						vout_EPSW[brd_i*nc + brd_j] *=  brd;
                	}
            	}          
        	}
			row+=stripes[i].size();
    	}
	#elif EXCUTION == 1
	// M.print_rowidx();
		M.stripe_static(3);
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	unsigned int strips = M._strips;
    	size_t vl;
   	 	for(int i=0; i<M._rows; i+=3){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	for(int k=0; k<strips && i+k<m; k++){
                	if(j<M._nnz_num[i+k]){
                    	Val          brd   = M._rowPtr[i+k][j];
                    	unsigned int brd_j = M._coordPtr[i+k][j];
						unsigned int brd_i = M._rowidx[i+k];
						Val 		 *tmp_a  = vout + brd_i * sc; // vout nr*sc
                    	Val          *tmp_b  = vin + brd_j*sc;    // vin nc*sc

						float32xm8_t vans = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
						float32xm8_t vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<sc; l+=vl){
	                        vl=vsetvli(sc-l, RVV_E32, RVV_M8);
							float32xm8_t va=vlev_float32xm8(tmp_a, vl);
                        	float32xm8_t vb=vlev_float32xm8(tmp_b, vl); 
							vtmp=vfmulvv_float32xm8(va, vb, vl);	
                        	vans = vfaddvv_float32xm8(vans, vtmp, vl);
					
                        	tmp_a += vl;
                        	tmp_b += vl;
                    	}
						
						float32xm8_t tmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
						tmp =  vfredsumvs_float32xm8(vans, tmp, vl);
						vout_EPSW[brd_i*nc + brd_j] +=  vfmvfs_float32xm8(tmp, 1);
						vout_EPSW[brd_i*nc + brd_j] *=  brd;
                	}
            	}          
        	}
			// row+=stripes[i].size();
    	}
	#elif EXCUTION == 0
		 for (int i = 0; i < M._rows; i++){   
			int num = M._nnz_num[i];
			int m = M._rowidx[i];

			size_t vl;
			float32xm8_t va, vb;
			uint32xm8_t vidx, vidx_t, vidx_s;

			for(int j=0; j<num; j+=vl){
				vl = vsetvli(num-j, RVV_E32, RVV_M8);
				vidx = vlev_uint32xm8(M._coordPtr[i]+j, vl);
				vidx_s = vsllvx_uint32xm8(vidx, 2, vl);

				for(int k=0; k<sc; k++){
					vb = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
					vidx_t = vmulvx_int32xm8(vidx, sc, vl);
					vidx_t = vaddvx_uint32xm8(vidx_t, k, vl);
					vidx_t = vsllvx_uint32xm8(vidx_t, 2, vl);
					va = vlxev_float32xm8(vin, vidx_t, vl);
				
					vb = vlxev_float32xm8(vout_EPSW+m*nc, vidx_s, vl);
					vb = vfmaccvf_float32xm8(vb, vout[m*sc+k], va, vl);
					
					vsxev_float32xm8(vout_EPSW+m*nc, vidx_s, vb, vl);
				}
			}
			for(int j=0; j<num; j++){
				vout_EPSW[m*nc + M._coordPtr[i][j]] *= M._rowPtr[i][j];
			}
		}
	#endif
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
