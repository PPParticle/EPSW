#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include "ELLPACK_P.h"
#include <riscv_vector.h>
using namespace std;

#define MAX(a,b) (((a)>(b))?(a):(b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float
#define MFACTOR (32)
#define THRESHOLD (16 * 1)
#define BH (128 * 1)
#define LOG_BH (7)
#define BW (128 * 1)
#define MIN_OCC (BW * 3 / 4)
#define STHRESHOLD (1024 / 2 * 1)
#define NTHREAD (68)
#define SC_SIZE (2048)
#define EXCUTION 2

double vari, avg;
double avg0[NTHREAD];
int sc, nr, nc, ne, gold_ne, npanel, mne, mne_nr;
int nr0;

int *csr_v;
int *csr_e, *csr_e0;
FTYPE *csr_ev, *csr_ev0;
struct ELLPACK_P M;
// int *mcsr_v;
int *mcsr_e; // can be short type
int *mcsr_cnt;
int *mcsr_chk;

int num_dense;

int *special;
int *special2;
int special_p;
char scr_pad_t[SC_SIZE]{0};
double p_elapsed;

void ready(int argc, char **argv)
{
	FILE *fp;
	int *loc;
	char buf[300];
	int nflag, sflag;
	int pre_count = 0, tmp_ne;
	int i;

	fprintf(stdout, "%s, \n", argv[1]);
	// sc = 128;
	sc = 1;
	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);
	printf("%d %d %d\n", nr, nc, ne);
	

	csr_v = new int[nr + 1]{0};	  // csr 行数组
	// csr_e0 = new int[ne]{0};	  // csr 列索引
	// csr_ev0 = new FTYPE[ne]{0.0}; // csr 元素数组
	csr_e = new int[ne]{0};	  // csr 列索引
	csr_ev = new FTYPE[ne]{0.0}; // csr 元素数组
	// csr 赋值
	for(i=0; i<nr+1; i++){
		fscanf(fp, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e[i]);
		// csr_ev[i] = (FTYPE)(rand() % 1048576) / 1048576;
		csr_ev[i] = (FTYPE)i;
	}
	fclose(fp);
	// for(int i=1; i<nr+1; i++){
	// 	// printf("line %d the number of line is %d\n", i, csr_v[i]-csr_v[i-1]);
	// 	for(int j=csr_v[i-1]; j<csr_v[i]; j++)
	// 		printf("%f ", csr_ev[j]);
	// 	printf("\n");
	// }

	// nr0 = nr;
	// nr = CEIL(nr, BH) * BH; // BH的整数倍
	// npanel = CEIL(nr, BH);

	// csr_e = new int[ne];
	// csr_ev = new FTYPE[ne];
	// fprintf(stdout, "%d,%d,%d\n", nr0, nc, ne);
}

void gen()
{
	// spmv rvv test    
	M.ELLPACKinitWithCSR(csr_v, csr_e, csr_ev, nr, nc, ne);
    M.ELLPACK_P_vreorder();
	// M.print_nnz();
	// M.print_rowidx();
	// M.print_ELLPACK_P();
    // M.stripe_adaptive();
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
		// vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
		vin[i] = (FTYPE)i;
	}

	// for(int i=0; i<nr*sc; i++){
	// 	cout << vin[i] << " ";
	// }
	// cout << endl;
	// cout << endl;
	// M.print_rowidx();
	// cout << endl;
	// M.print_ELLPACK_P();

	struct timeval starttime, endtime, starttime1, endtime1;
	#define ITER 128
	////begin
	gettimeofday(&starttime, NULL);
	// M.print_nnz();
	for(int loop=0;loop<ITER;loop++) {
	#if EXCUTION == 2
		vector<vector<unsigned int>> stripes;
		stripes = M.stripe_adaptive();

		// int cnt=0;
		// for(int i=0; i<stripes.size(); i++){
		// 	cnt+=stripes[i].size();
		// 	// for(int j=0; j<stripes[i].size(); j++){
		// 	// 	cout << stripes[i][j] << " ";
		// 	// }
		// 	cout << stripes[i].size() <<endl;
		// }
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
                    	Val          *tmp_b  = vin+brd_idx*n;
                    	Val          *tmp_c  = vout_EPSW + M._rowidx[row+k]*n;
						vfloat32m8_t vans = vfmv_s_f_f32m8(vans, 0.0, vl);
						vfloat32m8_t vtmp = vfmv_s_f_f32m8(vans, 0.0, vl);
                    	for(int l=n; l>0; l-=vl){
	                        vl=vsetvl_e32m8(l);
                        	vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
                        	vtmp=vfmul_vf_f32m8(vb, brd, vl);
                        	vb=vle32_v_f32m8(tmp_c, vl);
                        	vans = vfadd_vv_f32m8(vtmp, vb, vl);
                        	// vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
                        	vse32_v_f32m8(tmp_c, vans, vl);

                        	tmp_c += vl;
                        	tmp_b += vl;
                    	}
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
   	 	for(int i=0; i<M._rows; i+=strips){
        	for(int j=0; j<M._nnz_num[i]; j++){
            	for(int k=0; k<strips && i+k<m; k++){
                	if(j<M._nnz_num[i+k]){
                	// cout << "row " << i+k << endl;
                    	Val          brd     = M._rowPtr[i+k][j];
                    	unsigned int brd_idx = M._coordPtr[i+k][j];
                    	Val          *tmp_b  = vin+brd_idx*n;
                    	Val          *tmp_c  = vout_EPSW + M._rowidx[i+k]*n;
						vfloat32m8_t vans = vfmv_s_f_f32m8(vans, 0.0, vl);
						vfloat32m8_t vtmp = vfmv_s_f_f32m8(vans, 0.0, vl);
                	// cout << "brd " << brd << " brd_idx: "<<brd_idx <<" tmp_c: " << M._rowidx[i+k] << endl;;
                    	for(int l=n; l>0; l-=vl){
                        	vl=vsetvl_e32m8(l);
                        	vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
                        	vtmp=vfmul_vf_f32m8(vb, brd, vl);
                        	vb=vle32_v_f32m8(tmp_c, vl);
                        	vans = vfadd_vv_f32m8(vtmp, vb, vl);
                        	// vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
                        	vse32_v_f32m8(tmp_c, vans, vl);
                        	tmp_c += vl;
                        	tmp_b += vl;
                    	}
                	}
            	}          
        	}
    	}
	#elif EXCUTION == 0
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
            	Val          *tmp_b  = vin+brd_idx*n;
            	Val          *tmp_c  = vout_EPSW+M._rowidx[i]*n;
				vfloat32m8_t vans = vfmv_s_f_f32m8(vans, 0.0, vl);
				vfloat32m8_t vtmp = vfmv_s_f_f32m8(vans, 0.0, vl);
            	for(int k=n; k>0; k-=vl){
                	vl=vsetvl_e32m8(k);
                	vfloat32m8_t vb=vle32_v_f32m8(tmp_b, vl);
                	vtmp=vfmul_vf_f32m8(vb, brd, vl);
                	vb=vle32_v_f32m8(tmp_c, vl);
                	vans = vfadd_vv_f32m8(vtmp, vb, vl);
                	// vtmp=vfadd_vv_f32m8(vb, vtmp, vl);
                	vse32_v_f32m8(tmp_c, vans, vl);

                	tmp_c += vl;
                	tmp_b += vl;
            	}
        	}
    	}
	#endif
	}
	gettimeofday(&endtime, NULL);

		sc = 1;
	gettimeofday(&starttime1, NULL);
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
	
	elapsed[0] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0/ITER;
	elapsed[1] = ((endtime1.tv_sec-starttime1.tv_sec)*1000000 + endtime1.tv_usec-starttime1.tv_usec)/1000000.0/ITER;

	

	for(int i=0; i<5; i++){
		cout << "(" << vout[i] <<"," << vout_EPSW[i]<< ")" << " ";
	}


#define VALIDATE
#if defined VALIDATE
	// validate
	int num_diff = 0;
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
			num_diff++;
		}
		// fprintf(stdout, "%5.2f %5.2f\n", p1, p2);
	}
	fprintf(stdout, "diff : %f\n", (double)num_diff / (nr * sc) * 100);
#endif
// 	fprintf(stdout, "%f,", (double)ne * 2 * 8 * 1 / elapsed[0] / 1000000000);
// 	fprintf(stdout, "%f,%f,", (double)ne * 2 * 32 * 1 / elapsed[1] / 1000000000, (double)ne * 2 * 128 * 1 / elapsed[2] / 1000000000);
// 	fprintf(stdout, "%f\n", p_elapsed / (elapsed[2] / 1));
	fprintf(stdout, "%f,", elapsed[0] / 1000000000);
	fprintf(stdout, "%f,", elapsed[1] / 1000000000);
	cout << endl;
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	
	gen();
	mprocess();
}
