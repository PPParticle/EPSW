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
}

void gen()
{
	M.ELLPACKinitWithCSR(csr_v, csr_e, csr_ev, nr, nc, ne);
}

void mprocess()
{
	double elapsed[5];
	FTYPE *vin, *vout, *vout_EPSW_s, *vout_EPSW_ada, *vout_ELLPACK, *vout_ELLPACK_reorder;
	vin = new FTYPE[nc * sc]{0};
	vout =  new FTYPE[nr * sc]{0};
	vout_ELLPACK = new FTYPE[nr * sc]{0};
	vout_ELLPACK_reorder =    new FTYPE[nr * sc]{0};
	vout_EPSW_s =    new FTYPE[nr * sc]{0};
	vout_EPSW_ada =     new FTYPE[nr * sc]{0};
	

	for (int i = 0; i < nc * sc; i++)
	{
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// vin[i] = (FTYPE)i;
	}

	struct timeval starttime, endtime, starttime_EPSW_s, endtime_EPSW_s,
	starttime_ELLPACK,endtime_ELLPACK,starttime_ELLPACK_reorder, endtime_ELLPACK_reorder,starttime_EPSW_ada, endtime_EPSW_ada;

	#define ITER 1
	////begin
	gettimeofday(&starttime_ELLPACK, NULL);
	for(int loop=0;loop<ITER;loop++) {
    	unsigned int p=M._cols;
    	unsigned int n=sc;
		float32xm8_t vtmp, vb;
    	size_t vl;
    	for(int i=0; i<M._rows; i++){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	Val          brd   = M._rowPtr[i][j];
            	unsigned int brd_j = M._coordPtr[i][j];
				unsigned int brd_i = M._rowidx[i];
            	Val          *tmp_b  = vin+brd_j*n;
            	Val          *tmp_c  = vout_ELLPACK+brd_i*n;
				vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
            	for(int k=0; k<n; k+=vl){
                	vl = vsetvli(n-k, RVV_E32, RVV_M8);
                	vb = vlev_float32xm8(tmp_b+k, vl);
					vtmp = vlev_float32xm8(tmp_c+k, vl);
                	vtmp = vfmaccvf_float32xm8(vtmp, brd, vb, vl);
                	vsev_float32xm8(tmp_c+k, vtmp, vl);
            	}
        	}
    	}
	}
	gettimeofday(&endtime_ELLPACK, NULL);	


	M.ELLPACK_P_vreorder();
	gettimeofday(&starttime_ELLPACK_reorder, NULL);
	for(int loop=0;loop<ITER;loop++) {
		unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
		float32xm8_t vtmp, vb;
    	size_t vl;
    	for(int i=0; i<M._rows; i++){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	Val          brd   = M._rowPtr[i][j];
            	unsigned int brd_j = M._coordPtr[i][j];
				unsigned int brd_i = M._rowidx[i];
				Val          *tmp_b  = vin+brd_j*n;
            	Val          *tmp_c  = vout_ELLPACK_reorder+brd_i*n;
				vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
            	for(int k=0; k<n; k+=vl){
                	vl = vsetvli(n-k, RVV_E32, RVV_M8);
                	vb = vlev_float32xm8(tmp_b+k, vl);
					vtmp = vlev_float32xm8(tmp_c+k, vl);
                	vtmp = vfmaccvf_float32xm8(vtmp, brd, vb, vl);
                	vsev_float32xm8(tmp_c+k, vtmp, vl);
            	}
        	}
    	}
	}
	gettimeofday(&endtime_ELLPACK_reorder, NULL);

	gettimeofday(&starttime_EPSW_s, NULL);
	for(int loop=0;loop<ITER;loop++) {
	M.stripe_static(3);
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	unsigned int strips = M._strips;
		float32xm8_t vtmp, vb;
    	size_t vl;
   	 	for(int i=0; i<M._rows; i+=strips){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	for(int k=0; k<strips && i+k<m; k++){
                	if(j<M._nnz_num[i+k]){
                    	Val          brd     = M._rowPtr[i+k][j];
                    	unsigned int brd_j = M._coordPtr[i+k][j];
                    	Val          *tmp_b  = vin+brd_j*n;
                    	Val          *tmp_c  = vout_EPSW_s + M._rowidx[i+k]*n;
						vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<n; l+=vl){
                        	vl = vsetvli(n-l, RVV_E32, RVV_M8);
                			vb = vlev_float32xm8(tmp_b+l, vl);
							vtmp = vlev_float32xm8(tmp_c+l, vl);
                			vtmp = vfmaccvf_float32xm8(vtmp, brd, vb, vl);
                			vsev_float32xm8(tmp_c+l, vtmp, vl);
                    	}
                	}
            	}          
        	}
    	}
	}
	gettimeofday(&endtime_EPSW_s, NULL);


	gettimeofday(&starttime_EPSW_ada, NULL);
	for(int loop=0;loop<ITER;loop++) {
		vector<vector<unsigned int>> stripes;
		stripes = M.stripe_adaptive();
    	unsigned int m=M._rows;
    	unsigned int p=M._cols;
    	unsigned int n=sc;
    	unsigned int strips = M._strips;
		float32xm8_t vtmp, vb;
    	size_t vl;
		int row=0;
    	for(int i=0; i<stripes.size(); i++){
        	for(int j=0; j<M._nnz_num[row]; j++){
            	for(int k=0; k<stripes[i].size() && row+k<m; k++){
                	if(j<M._nnz_num[row+k]){
                    	Val          brd     = M._rowPtr[row+k][j];
                    	unsigned int brd_j = M._coordPtr[row+k][j];
                    	Val          *tmp_b  = vin+brd_j*n;
                    	Val          *tmp_c  = vout_EPSW_ada + M._rowidx[row+k]*n;
						vtmp = vfmvvf_float32xm8(0.0, vsetvli_max(RVV_E32, RVV_M8));
                    	for(int l=0; l<n; l+=vl){
                        	vl = vsetvli(n-l, RVV_E32, RVV_M8);
                			vb = vlev_float32xm8(tmp_b+l, vl);
							vtmp = vlev_float32xm8(tmp_c+l, vl);
                			vtmp = vfmaccvf_float32xm8(vtmp, brd, vb, vl);
                			vsev_float32xm8(tmp_c+l, vtmp, vl);
                    	}
                	}
            	}          
        	}
			row+=stripes[i].size();
    	}
	}
	gettimeofday(&endtime_EPSW_ada, NULL);

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
	
	elapsed[0] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/ITER;
	elapsed[1] = ((endtime_ELLPACK.tv_sec-starttime_ELLPACK.tv_sec)*1000000 + endtime_ELLPACK.tv_usec-starttime_ELLPACK.tv_usec)/ITER;
	elapsed[2] = ((endtime_ELLPACK_reorder.tv_sec-starttime_ELLPACK_reorder.tv_sec)*1000000 + endtime_ELLPACK_reorder.tv_usec-starttime_ELLPACK_reorder.tv_usec)/ITER;
	elapsed[3] = ((endtime_EPSW_s.tv_sec-starttime_EPSW_s.tv_sec)*1000000 + endtime_EPSW_s.tv_usec-starttime_EPSW_s.tv_usec)/ITER;
	elapsed[4] = ((endtime_EPSW_ada.tv_sec-starttime_EPSW_ada.tv_sec)*1000000 + endtime_EPSW_ada.tv_usec-starttime_EPSW_ada.tv_usec)/ITER;

	// for (int i = 0; i < nr*sc; i++)
	// {
	// 	FTYPE p1 = vout[i];
	// 	FTYPE p2 = vout_ELLPACK[i];
	// 	FTYPE p3 = vout_ELLPACK_reorder[i];
	// 	FTYPE p4 = vout_EPSW_s[i];
	// 	FTYPE p5 = vout_EPSW_ada[i];
	// }
	FTYPE p1 = elapsed[0];
	FTYPE p2 = elapsed[1];
	FTYPE p3 = elapsed[2];
	FTYPE p4 = elapsed[3];
	FTYPE p5 = elapsed[4];
	float min = MIN(p1, MIN(p2, MIN(p3, MIN(p4, p5))));
	if(p1==min)
	fprintf(stdout, "%f, %f, %f, %f, %f scalar\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, 
						elapsed[2]/(double)ITER , elapsed[3]/(double)ITER, elapsed[4]/(double)ITER );
	if(p2==min)
				fprintf(stdout, "%f, %f, %f, %f, %f ell\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, 
					elapsed[2]/(double)ITER , elapsed[3]/(double)ITER, elapsed[4]/(double)ITER );
	if(p3==min)	
				fprintf(stdout, "%f, %f, %f, %f, %f order\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, 
					elapsed[2]/(double)ITER , elapsed[3]/(double)ITER, elapsed[4]/(double)ITER );
	if(p4==min)		
				fprintf(stdout, "%f, %f, %f, %f, %f static\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, 
					elapsed[2]/(double)ITER , elapsed[3]/(double)ITER, elapsed[4]/(double)ITER );
	if(p5==min)		
				fprintf(stdout, "%f, %f, %f, %f, %f adaptive\n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, 
					elapsed[2]/(double)ITER , elapsed[3]/(double)ITER, elapsed[4]/(double)ITER );
				
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
