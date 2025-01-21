// 似乎使用 simd 是不可行的
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "ELLPACK_P.h"
#include <riscv_vector.h>
#include <vector>
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float

int nr, nc, ne;
int nr_l, nc_l, ne_l, nr_r, nc_r, ne_r;
int *csr_v;
int *csr_e;
FTYPE *csr_ev;

int *csr_vr;
int *csr_er;
FTYPE *csr_evr;
int col=0;

ELLPACK_P M1;
ELLPACK_P m1;
ELLPACK_P M3;

char *dir_l, *dir_r;
char *app;

void ready(int argc, char **argv)
{
	FILE *fp_l, *fp_r;
	int i;

	app = argv[0];
	dir_l = argv[1];
	dir_r = argv[2];
	fp_l = fopen(dir_l, "r");
	fp_r = fopen(dir_r, "r");

	fscanf(fp_l, "%d, %d, %d", &nr_l, &nc_l, &ne_l);
	fscanf(fp_r, "%d, %d, %d", &nr_r, &nc_r, &ne_r);

	if(nc_l != nr_r) cerr << "wrong shape of matrix" << endl;
	
	nr = nr_l;
	nc = nc_r;
	ne = ne_l;

	csr_v = new int[nr_l + 1]{0};	  
	csr_e = new int[ne_l]{0};	   
	csr_ev = new FTYPE[ne_l]{0.0}; 

	csr_vr = new int[nr_r + 1]{0};	  
	csr_er = new int[ne_r]{0};	   
	csr_evr = new FTYPE[ne_r]{0.0}; 

	for (i = 0; i < nr_l+1; i++){
		fscanf(fp_l, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne_l; i++){
		fscanf(fp_l, "%d", &csr_e[i]);
		// csr_ev[i] = (FTYPE)(rand() & 1048575) / 1048576;
		csr_ev[i] = 1;
	}

	for (i = 0; i < nr_r+1; i++){
		fscanf(fp_r, "%d ", &csr_vr[i]);
	}

	for (i = 0; i < ne_r; i++){
		fscanf(fp_r, "%d", &csr_er[i]);
		// csr_ev[i] = (FTYPE)(rand() & 1048575) / 1048576;
		csr_evr[i] = 1;
	}
	fclose(fp_l);
	fclose(fp_r);
	fprintf(stdout, "%d %d %d \n", ne, nr, nc);
}

void mprocess()
{
	double elapsed[2];
	FTYPE *vout_scalar = new FTYPE[nr*nr]{0};	
	FTYPE *vout_vec = new FTYPE[nr*nr]{0.0};

	struct timeval starttime_scalar, endtime_scalar, starttime_vec, endtime_vec;
	#define ITER 1

	////begin
	gettimeofday(&starttime_vec, NULL);
int cnt=0;
	for(int loop=0; loop<ITER; loop++) {
		size_t vl = vsetvlmax_e32m1();
		vfloat32m1_t va, vb, vtmp;
		vfloat32m1_t vans;
		vuint32m1_t vaidx, vbidx;
		vbool32_t vand;

		for (int i=0; i<nr; i++){ // row i
			int ni = csr_v[i+1] - csr_v[i];

			if(ni>=4){
				
				int itri = ni/vl;
				int remaini = csr_v[i] + vl*itri;

				for(int l=0; l<itri; l++){ // first batch of a line
					vaidx = vle32_v_u32m1((unsigned int*)&csr_e[csr_v[i] + l*vl], vl);
					va = vle32_v_f32m1(&csr_ev[csr_v[i] + l*vl], vl);
	
					for(int j=0; j<nc; j++){ // col j
						int nj = csr_v[j+1] - csr_v[j];

						if(nj>=4){
							int itrj = nj/vl;
							int remainj = csr_v[j] + vl*itrj;
							vans = vfmv_v_f_f32m1(0.0, vsetvlmax_e32m1());
							vtmp = vfmv_v_f_f32m1(0.0, vsetvlmax_e32m1());

							for(int jk=0; jk<itrj; jk++){ // first batch of a col
								vbidx = vle32_v_u32m1((unsigned int*)&csr_e[csr_v[j] + jk*vl], vl);
								vb = vle32_v_f32m1(&csr_ev[csr_v[j] + jk*vl], vl);

								vand = vmseq_vv_u32m1_b32(vaidx, vbidx, vl);
								vtmp = vfmacc_vv_f32m1_m(vand, vtmp, va, vb, vl);  

								cout << i << " " << j << endl;
								float* test = new float[vl];
								vse32_v_f32m1(test, vtmp, vl);
								for(int m=0; m<vl; m++){
									cout << test[m] << " ";
								}
								cout << endl;
							}
							vans = vfredusum_vs_f32m1_f32m1(vans, vtmp, vans, vl);
							vout_vec[i*nc+j] += vfmv_f_s_f32m1_f32(vans);	

							// tail of vec of right
							for(int x = csr_v[i]+l*vl; x<csr_v[i]+(l+1)*vl;){
								for(remainj=csr_v[j] + vl*itrj; remainj<csr_v[j+1]; ){
									if(x>=csr_v[i+1])
										break;
									if(csr_e[x] == csr_e[remainj]){
										vout_vec[i*nc+j] += csr_ev[x]*csr_ev[remainj];
										x++;
										remainj++;
									}
									else if(csr_e[x] > csr_e[remainj]){
										remainj++;
									}
									else x++;
								}
							}
							
						}
						else{
							for(int x = csr_v[i]+l*vl; x<csr_v[i]+(l+1)*vl; x++){
								for(int nj=csr_v[j]; nj<csr_v[j+1];){
									if(x>=csr_v[i+1])
										break;
									if(csr_e[x] == csr_e[nj]){
										vout_vec[i*nc+j] += csr_ev[x]*csr_ev[nj];
										x++;
										nj++;
									}
									else if(csr_e[x] > csr_e[nj]){
										nj++;
									}
									else x++;
								}
							}
						}
						// cout << n_1 << endl;
					}
				}
				// tail  of row
				for(;remaini<csr_v[i+1]; remaini++){
					// cout << csr_e[remaini] << " ";
					for(int j=0; j<nc; j++){
						for(int nj=csr_v[j]; nj<csr_v[j+1]; nj++){
							if(remaini>=csr_v[i+1])
								break;
							if(csr_e[remaini] == csr_e[nj]){
								vout_vec[i*nc+j] += csr_ev[remaini]*csr_ev[nj];
								remaini++;
								nj++;
								}
							else if(csr_e[remaini] > csr_e[nj])
								nj++;
							else remaini++;
						}
					}
				}
			}

			else{
				for(int ni=csr_v[i]; ni<csr_v[i+1]; ni++){
					for(int j=0; j<nc; j++){
						for(int nj=csr_vr[j]; nj<csr_v[j+1]; nj++){
							if(ni>=csr_v[i+1])
								break;
							if(csr_e[ni] == csr_e[nj]){
								vout_vec[i*nc+j] += csr_ev[ni]*csr_ev[nj];
								ni++;
								nj++;
							}
							else if(csr_e[ni] > csr_e[nj]){
								nj++;
							}
							else ni++;
						}
					}
				}
			}
		}
	}
	gettimeofday(&endtime_vec, NULL);

	gettimeofday(&starttime_scalar, NULL);
	for (int loop = 0; loop < ITER; loop++){
		for (int i = 0; i < nr; i++){
				for(int j=0; j<nc; j++){
					int ni = csr_v[i];
					for(int nj=csr_v[j]; nj<csr_v[j+1]; ){
						if(ni>=csr_v[i+1]) {
							break;
						}
						if(csr_e[ni]==csr_e[nj]){
							vout_scalar[i*nc+j] += csr_ev[ni]*csr_ev[nj];
							ni++;
							nj++;
						}
						else if(csr_e[ni]>csr_e[nj])
							nj++;
						else 
							ni++;
					}
				}
		}
	}
	gettimeofday(&endtime_scalar, NULL);
	
	elapsed[0] = ((endtime_vec.tv_sec-starttime_vec.tv_sec)*1000000 + endtime_vec.tv_usec-starttime_vec.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_scalar.tv_sec-starttime_scalar.tv_sec)*1000000 + endtime_scalar.tv_usec-starttime_scalar.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	
	for (int i = 0; i < nr*nc; i++){	

			FTYPE p1 = vout_vec[i]/ITER;
			FTYPE p2 = vout_scalar[i]/ITER;

			if (p1 < 0)
				p1 *= -1;
			if (p2 < 0)
				p2 *= -1;
			FTYPE diff;
			diff = p1 - p2;
			if (diff < 0)
				diff *= -1;
			// if (MAX(p1,p2) !=0 && diff / MAX(p1, p2) > 0.01){
				// fprintf(stdout, "%d %f %f %s %s\n", i, p1, p2, app, dir);
				fprintf(stdout, "%d %f %f \n", i, p1, p2);
				// return;
			// }
	}
	
#endif

	fprintf(stdout, "%f|%f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	mprocess();
}
