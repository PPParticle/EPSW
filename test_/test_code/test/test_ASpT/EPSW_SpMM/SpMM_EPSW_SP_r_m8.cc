#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "ELLPACK_P_0.7.h"
#include <riscv-vector.h>
using namespace std;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float
#define EXECUTION 2

int sc, nr, nc, ne;
int *csr_v, *csr_e;
FTYPE *csr_ev;
struct ELLPACK_P M;

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;

	//fprintf(stdout, "%s, \n", argv[1]);
	sc = atoi(argv[2]);

	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);
	printf("%d,%d,%d,%d ", nr, nc, ne, sc);


	csr_v = new int[nr + 1]{0};	  
	csr_e = new int[ne]{0};	 
	csr_ev = new FTYPE[ne]{0.0}; 

	for(i=0; i<nr+1; i++){
		fscanf(fp, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e[i]);
		csr_ev[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// csr_ev[i] = (FTYPE)i;
	}
	fclose(fp);
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
	FTYPE *vin, *vout, *vout_EPSW;
	vin = new FTYPE[nc * sc]{0};
	vout = new FTYPE[nr * sc]{0};
	vout_EPSW = new FTYPE[nr * sc]{0};

	for (int i = 0; i < nc * sc; i++)
	{
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// vin[i] = (FTYPE)i;
	}

	struct timeval starttime, endtime, starttime_EPSW, endtime_EPSW;
	#define ITER 2
	////begin
	gettimeofday(&starttime_EPSW, NULL);
	for(int loop=0;loop<ITER;loop++) {
	#if EXECUTION == 2
		vector<vector<unsigned int>> stripes;
		stripes = M.stripe_adaptive();
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	unsigned int strips = M._strips;
    	size_t vl;
		int row=0;
    	for(int i=0; i<stripes.size(); i++){
        	for(int j=0; j<M._nnz_num[row]; j++){
            	for(int k=0; k<stripes[i].size() && row+k<m; k++){
                	if(j<M._nnz_num[row+k]){
                    	Val          brd     = M._rowPtr[row+k][j];
                    	unsigned int brd_idx = M._coordPtr[row+k][j];
                    	Val          *tmp_c  = vout_EPSW + M._rowidx[row+k]*n;
						float32xm8_t vb;
						float32xm8_t vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<n; l+=vl){
	                        vl=vsetvli(n-l, RVV_E32, RVV_M8);
                        	vb=vlev_float32xm8(vin+brd_idx*n+l, vl);
                        	vtmp=vlev_float32xm8(tmp_c+l, vl);
                        	vtmp=vfmaccvf_float32xm8(vtmp, brd, vb, vsetvli_max(RVV_E32, RVV_M8));
                        	vsev_float32xm8(tmp_c+l, vtmp, vl);
                    	}
                	}
            	}          
        	}
			row+=stripes[i].size();

    	}
	#elif EXECUTION == 1
	// M.print_rowidx();
		M.stripe_static(3);
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	unsigned int strips = M._strips;
    	size_t vl;
   	 	for(int i=0; i<M._rows; i+=strips){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	for(int k=0; k<strips && i+k<m; k++){
                	if(j<M._nnz_num[i+k]){
                		Val          brd     = M._rowPtr[i+k][j];
                    	unsigned int brd_idx = M._coordPtr[i+k][j];
                    	Val          *tmp_c  = vout_EPSW + M._rowidx[i+k]*n;
						float32xm8_t vb;
						float32xm8_t vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<n; l+=vl){
	                        vl=vsetvli(n-l, RVV_E32, RVV_M8);
                        	vb=vlev_float32xm8(vin+brd_idx*n+l, vl);
                        	vtmp=vlev_float32xm8(tmp_c+l, vl);
                        	vtmp=vfmaccvf_float32xm8(vtmp, brd, vb, vsetvli_max(RVV_E32, RVV_M8));
                        	vsev_float32xm8(tmp_c+l, vtmp, vl);
                    	}
                	}
            	}          
        	}
    	}
	#elif EXECUTION == 0
		unsigned int m=M._rows;
    	// for(int i=0 ; i<M._rows; i++){
    	//     cout << M._nnz_num[i] << " ";
    	// }
    	// cout << endl;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	size_t vl;
    	for(int i=0; i<M._rows; i++){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	Val          brd     = M._rowPtr[i][j];
            	unsigned int brd_idx = M._coordPtr[i][j];
            	Val          *tmp_c  = vout_EPSW+M._rowidx[i]*n;
				float32xm8_t vb;
				float32xm8_t vtmp = vfmvsf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
            	for(int l=0; l<n; l+=vl){
	                vl=vsetvli(n-l, RVV_E32, RVV_M8);
                   	vb=vlev_float32xm8(vin+brd_idx*n+l, vl);
                  	vtmp=vlev_float32xm8(tmp_c+l, vl);
                    vtmp=vfmaccvf_float32xm8(vtmp, brd, vb, vsetvli_max(RVV_E32, RVV_M8));
                 	vsev_float32xm8(tmp_c+l, vtmp, vl);
                }
        	}
    	}
	#endif
	}
	gettimeofday(&endtime_EPSW, NULL);

	gettimeofday(&starttime, NULL);
	for(int loop=0;loop<ITER;loop++) {
		for(int i=0; i<nr; i++){
			for(int j=0; j<sc; j++){
				for(int k=csr_v[i]; k<csr_v[i+1]; k++){
					vout[i*sc+j] += csr_ev[k] * vin[csr_e[k]*sc + j];
				}
			}
		}
	}
	gettimeofday(&endtime, NULL);
	
	elapsed[0] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_EPSW.tv_sec-starttime_EPSW.tv_sec)*1000000 + endtime_EPSW.tv_usec-starttime_EPSW.tv_usec)/1000000.0;


#define VALIDATE
#if defined VALIDATE

	for (int i = 0; i < nr * sc; i++)
	{
		FTYPE p1 = vout[i];
		FTYPE p2 = vout_EPSW[i];

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
			cout << i << " " << p1 << " " << p2 << endl;
			break;
		}
	}
#endif
	fprintf(stdout, "%f,%f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
