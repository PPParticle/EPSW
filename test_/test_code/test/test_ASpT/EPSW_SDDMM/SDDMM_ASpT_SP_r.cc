#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <riscv-vector.h>
using namespace std;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float
#define THRESHOLD (16)
#define BH (32 * 1)
#define BW (32 * 1)
#define MIN_OCC (BW * 3 / 4)
#define NTHREAD (32)
#define SC_SIZE (256)

int sc, nr, nc, ne, npanel;
int nr0;

int *csr_v;
int *csr_e, *csr_e0;
FTYPE *csr_ev, *csr_ev0;

int *mcsr_e; // can be short type
int *mcsr_cnt;
int *mcsr_chk;

int *special;
int *special2;
int special_p;
int scr_pad_t[SC_SIZE]{0}; // 提前定义，应该大于等于列数

int **csr_e1 = new int *[2];
short **coo = new short *[2];
char *dir;

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;
	dir = argv[1];
	sc = atoi(argv[2]);
	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);
		fprintf(stdout, "%d %d %d %d ", ne, nr, nc, sc);

	// printf("%d %d %d\n", nr, nc, ne);
	

	nr0 = nr;
	nr = CEIL(nr, BH) * BH; // BH的整数倍
	npanel = CEIL(nr, BH);

	// 初始化 ne 个元素，每个元素用 struct 进行表示
	csr_v = new int[nr + 1]{0};	  // csr 行数组
	csr_e0 = new int[ne]{0};	  // csr 列索引
	csr_ev0 = new FTYPE[ne]{0.0}; // csr 元素数组

	// csr 赋值
	for (i = 0; i < nr0+1; i++)
	{
		fscanf(fp, "%d ", &csr_v[i]);
	}
	for(;i<nr+1;i++){
		csr_v[i] = csr_v[i-1];
	}
	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e0[i]);
		csr_ev0[i] = (FTYPE)(rand() % 1048576) / 1048576;;
	}

	csr_e = new int[ne]{0};
	csr_ev = new FTYPE[ne+256]{0.0};
	fclose(fp);
}


void gen()
{
	special = new int[ne]{0};
	special2 = new int[ne]{0};
	mcsr_cnt = new int[npanel + 1]{0};
	mcsr_chk = new int[npanel + 1]{0};
	mcsr_e = new int[ne]{0};
	for (int i = 0; i < 2; i++)
	{
		csr_e1[i] = new int[ne]{0};
		coo[i] = new short[ne]{0};
	}


	int bv_size = CEIL(nc, 16);
	unsigned int *bv = new unsigned int [bv_size];

	struct timeval starttime0;
	struct timeval starttime, endtime;

	// 可能存在问题
	for(int row_panel=0; row_panel<nr/BH; row_panel++){
		for(int i=row_panel*BH; i<(row_panel+1)*BH; i++){
			for(int j=csr_v[i]; j<csr_v[i+1]; j++) {
				csr_e1[0][j] = csr_e0[j];
				// fprintf(stdout, "%d ", j);
			}
		}

	}
	int cnt = 0;

	for (int row_panel=0; row_panel<nr/BH; row_panel++){
		int i, j, t_sum = 0;
		memset(scr_pad_t, 0, sizeof(char)*SC_SIZE);
		for(i=row_panel*BH; i<(row_panel+1)*BH; i++){
			for(j=csr_v[i]; j<csr_v[i+1]; j++){
				coo[0][j] = (i & (BH-1));			 // coo[0] 按组存储所有元素的行坐标在当前组中的位置 0~127
				int k = (csr_e0[j] & (SC_SIZE-1));   // 保证元素的列坐标不会越界 (SC_SIZE) 2^10
				// 统计同一列元素个数，如果一列元素个数总和超过 THRESHOLD t_sum++
				if(scr_pad_t[k] < THRESHOLD){       // panel 中行数超过一般就认为稠密列
					if(scr_pad_t[k] == THRESHOLD - 1)
						t_sum++;
					scr_pad_t[k]++;
				}
			}
		}

		if (t_sum < MIN_OCC) {       // 存在 3/4 稠密列，认为稠密块
			mcsr_chk[row_panel] = 1;
			mcsr_cnt[row_panel+1] = 1;
			continue;
		} 

		// sorting(merge sort)
		int flag = 0;
		for(int stride = 1; stride <= BH/2; stride *= 2, flag=1-flag) {
			for(int pivot = row_panel*BH; pivot < (row_panel+1)*BH; pivot += stride*2) {
				int l1, l2;
				for(i = l1 = csr_v[pivot], l2 = csr_v[pivot+stride]; l1 < csr_v[pivot+stride] && l2 < csr_v[pivot+stride*2]; i++) {
					if(csr_e1[flag][l1] <= csr_e1[flag][l2]) {
						coo[1-flag][i] = coo[flag][l1];
						csr_e1[1-flag][i] = csr_e1[flag][l1++];
					}
					else {
						coo[1-flag][i] = coo[flag][l2];
						csr_e1[1-flag][i] = csr_e1[flag][l2++];	
					}
				}
				while(l1 < csr_v[pivot+stride]) {
					coo[1-flag][i] = coo[flag][l1];
					csr_e1[1-flag][i++] = csr_e1[flag][l1++];
				}
				while(l2 < csr_v[pivot+stride*2]) {
					coo[1-flag][i] = coo[flag][l2];
					csr_e1[1-flag][i++] = csr_e1[flag][l2++];
				}
			}
		}		

		int weight=1;
		int cr=0;

		// dense bit extract (and mcsr_e making)
		for(i=csr_v[row_panel*BH]+1; i<csr_v[(row_panel+1)*BH]; i++) {
			if(csr_e1[flag][i-1] == csr_e1[flag][i]) weight++;
			else {
				if(weight >= THRESHOLD) {
					cr++;
				} 				
				weight = 1;
			}
		}
		if(weight >= THRESHOLD) {
			cr++;
		} 		
		mcsr_cnt[row_panel+1] = CEIL(cr,BW)+1;
	}

	for(int i=1; i<=npanel;i++) 
		mcsr_cnt[i] += mcsr_cnt[i-1];
	mcsr_e[BH * mcsr_cnt[npanel]] = ne; // mcsr_e[192] = ne;

	int tid = 0;
	for(int row_panel=0; row_panel<nr/BH; row_panel++, tid++) {
	    if(mcsr_chk[row_panel] == 0) {
			int i, j;
			int flag = 0;
			int cq=0, cr=0;
			for(int stride = 1; stride <= BH/2; stride*=2, flag=1-flag); // 恢复到之前的 flag
			int base = (mcsr_cnt[row_panel]*BH); // 0 96
			int mfactor = mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]; // 3 3
			int weight=1;

		// mcsr_e making
			for(i=csr_v[row_panel*BH]+1; i<csr_v[(row_panel+1)*BH]; i++){
				if(csr_e1[flag][i-1] == csr_e1[flag][i]) weight++;
				else{
					int reminder = (csr_e1[flag][i-1]&15); // csr_e1[flag][i-1] % 16  reminder 0-15
					if(weight >= THRESHOLD){
						cr++;
						bv[csr_e1[flag][i-1]>>4] |= (1<<reminder);  // bv 存储在哪一个2D块中的索引, 第csr_e1[flag][i-1]>>4的索引1<<reminder
						for(j=i-weight; j<=i-1; j++) {
							mcsr_e[base + coo[flag][j] * mfactor + cq + 1]++;
						}
					}else{
						bv[csr_e1[flag][i-1]>>4] &= (0xFFFFFFFF - (1<<reminder)); 
					} 
					if(cr == BW) { cq++; cr=0;}
					weight = 1;
				}
			}	

			int reminder = (csr_e1[flag][i-1]&15);
			// cout << weight << endl;
			if(weight >= THRESHOLD) {
				cr++;
				bv[csr_e1[flag][i-1]>>4] |= (1<<reminder);  // 表示该列是稠密的
				for(j=i-weight; j<=i-1; j++) {
					mcsr_e[base + coo[flag][j] * mfactor + cq + 1]++;
				}
			} else {
				bv[csr_e1[flag][i-1]>>4] &= (0xFFFFFFFF - (1<<reminder)); // 表示该列是稀疏的
			} 
		// reordering
			int delta = mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]; // 3 3
			int base0 = mcsr_cnt[row_panel]*BH; // 0 96
			for(i=row_panel*BH; i<(row_panel+1)*BH; i++) {
				int base = base0+(i-row_panel*BH)*delta;
				int dpnt = mcsr_e[base] = csr_v[i];

				for(int j=1;j<delta;j++) {
					mcsr_e[base+j] += mcsr_e[base+j-1];
				}

				int spnt=mcsr_e[mcsr_cnt[row_panel]*BH + (mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel])*(i - row_panel*BH + 1) - 1];
				for(j=csr_v[i]; j<csr_v[i+1]; j++) {
					int k = csr_e0[j];
					if((bv[k>>4]&(1<<(k&15)))) {
						csr_e[dpnt] = csr_e0[j];
						csr_ev[dpnt++] = csr_ev0[j];
					} else {
						csr_e[spnt] = csr_e0[j];
						csr_ev[spnt++] = csr_ev0[j];
					}
				}
			}
	   	} else {
			int base0 = mcsr_cnt[row_panel]*BH;
			memcpy(&mcsr_e[base0], &csr_v[row_panel*BH], sizeof(int)*BH);
			int bidx = csr_v[row_panel*BH];
			int bseg = csr_v[(row_panel+1)*BH] - bidx;
			memcpy(&csr_e[bidx], &csr_e0[bidx], sizeof(int)*bseg);
			memcpy(&csr_ev[bidx], &csr_ev0[bidx], sizeof(FTYPE)*bseg);
	   	}
	}

	free(bv);
	free(csr_e1);
	free(coo);
}


void mprocess()
{
	double elapsed[3];
  	FTYPE *vin, *vout, *vout_ASpT_vec, *vout_ASpT_scalar, *vout_scalar;
    vin = new FTYPE[nc * sc]{0};
	vout = new FTYPE[nr * sc]{0};
	vout_ASpT_vec = new FTYPE[ne]{0.0};
	vout_ASpT_scalar = new FTYPE[ne]{0.0};
	vout_scalar = new FTYPE[ne]{0.0};

	struct timeval starttime, endtime, 
	starttime_ASpT_vec, endtime_ASpT_vec,
	starttime_ASpT_scalar, endtime_ASpT_scalar;

	for (int i = 0; i < nc * sc; i++)
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
	
	for (int i = 0; i < nr * sc; i++)
		vout[i] = (FTYPE)(rand() % 1048576) / 1048576;
	

#define ITER 10
	gettimeofday(&starttime_ASpT_vec, NULL);
	for(int loop=0;loop<ITER;loop++) {
		for(int row_panel=0; row_panel<nr/BH; row_panel ++){
			int stride;
			for(stride = 0; stride < mcsr_cnt[row_panel+1]-mcsr_cnt[row_panel]-1; stride++){
				for(int i=row_panel*BH; i<(row_panel+1)*BH; i++){	
					int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 					int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
					int interm = loc1 + (((loc2 - loc1)>>3)<<3);
					int j;
					size_t vl=vsetvli_max(RVV_E32, RVV_M2);
					// o*a*b = c
					float32xm2_t va, vb;
					uint32xm2_t vidx, vidx_t;
					vb = vfmvvf_float32xm2(0.0, vl);
					for(j=loc1; j<interm; j+=8){
						vidx = vlev_uint32xm2((unsigned int *)csr_e+j, vl);
						vidx = vmulvx_int32xm2(vidx, sc, vl);
						for(int k=0; k<sc; k++){
							vb = vlev_float32xm2(vout_ASpT_vec+j, vl);
							vidx_t = vaddvx_uint32xm2(vidx, k, vl);
							vidx_t = vsllvx_uint32xm2(vidx_t, 2, vl);
							va = vlxev_float32xm2(vin, vidx_t, vl);
							vb = vfmaccvf_float32xm2(vb, vout[i*sc+k], va, vl);
							vsev_float32xm2(vout_ASpT_vec+j, vb, vl);
						}
						for(int k=0; k<8; k++)
							vout_ASpT_vec[j+k] *= csr_ev[j+k];
					}
					for(; j<loc2; j++){
						for(int k=0; k<sc; k++)
                            vout_ASpT_vec[j] = vout_ASpT_vec[j] + vin[csr_e[j]*sc + k] * vout[i*sc + k];
						vout_ASpT_vec[j] *= csr_ev[j];
					}	
				}
			}
			// sparse
			for(int i=row_panel*BH; i<(row_panel+1)*BH; i++) {
				int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
				int interm = loc1 + (((loc2 - loc1)>>3)<<3);
				int j;
				size_t vl=vsetvli_max(RVV_E32, RVV_M2);
					// o*a*b = c
					float32xm2_t va, vb;
					uint32xm2_t vidx, vidx_t;
					vb = vfmvvf_float32xm2(0.0, vl);
					for(j=loc1; j<interm; j+=8){
						vidx = vlev_uint32xm2((unsigned int *)csr_e+j, vl);
						vidx = vmulvx_int32xm2(vidx, sc, vl);
						for(int k=0; k<sc; k++){
							vb = vlev_float32xm2(vout_ASpT_vec+j, vl);
							vidx_t = vaddvx_uint32xm2(vidx, k, vl);
							vidx_t = vsllvx_uint32xm2(vidx_t, 2, vl);
							va = vlxev_float32xm2(vin, vidx_t, vl);
							vb = vfmaccvf_float32xm2(vb, vout[i*sc+k], va, vl);
							vsev_float32xm2(vout_ASpT_vec+j, vb, vl);
						}
					for(int k=0; k<8; k++) {
						vout_ASpT_vec[j+k] *= csr_ev[j+k];
					}
				}
				for(; j<loc2; j++){
					for(int k=0; k<sc; k++) 
               			vout_ASpT_vec[j] = vout_ASpT_vec[j] + vin[csr_e[j]*sc + k] * vout[i*sc + k];
					vout_ASpT_vec[j] *= csr_ev[j];
				}
			}
		}
		
	}
	gettimeofday(&endtime_ASpT_vec, NULL);
	
	gettimeofday(&starttime_ASpT_scalar, NULL);
	for(int loop=0;loop<ITER;loop++) {
		for(int row_panel=0; row_panel<nr/BH; row_panel ++){
			int stride;
			for(stride = 0; stride < mcsr_cnt[row_panel+1]-mcsr_cnt[row_panel]-1; stride++){
				for(int i=row_panel*BH; i<(row_panel+1)*BH; i++){	
					int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 					int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
					int j;
					for(j=loc1; j<loc2; j++){
						for(int k=0; k<sc; k++)
                        	vout_ASpT_scalar[j] =   vout_ASpT_scalar[j] +   vin[csr_e[j] * sc + k] * vout[i*sc + k];
						vout_ASpT_scalar[j] *= csr_ev[j];
					}
				}
			}
			// sparse
			for(int i=row_panel*BH; i<(row_panel+1)*BH; i++) {
				int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
				int j;
				for(j=loc1; j<loc2; j++){
					for(int k=0; k<sc; k++)
                        vout_ASpT_scalar[j] =   vout_ASpT_scalar[j] +   vin[csr_e[j] * sc + k] * vout[i*sc + k];
					vout_ASpT_scalar[j] *= csr_ev[j];
				}
			}
		}
		
	}
	gettimeofday(&endtime_ASpT_scalar, NULL);

	// matrix (hamard*) (vout * vin)
	// matrix[i,j] = matrix[i, csr_e[j]] = sum[k|0:sc-1](vout[i, k] * vin[csr_e[j]*sc + k])
	gettimeofday(&starttime, NULL);
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
	gettimeofday(&endtime, NULL);




    elapsed[0] = ((endtime_ASpT_vec.tv_sec-starttime_ASpT_vec.tv_sec)*1000000 + endtime_ASpT_vec.tv_usec-starttime_ASpT_vec.tv_usec)/1000000.0;
    elapsed[1] = ((endtime_ASpT_scalar.tv_sec-starttime_ASpT_scalar.tv_sec)*1000000 + endtime_ASpT_scalar.tv_usec-starttime_ASpT_scalar.tv_usec)/1000000.0;
    elapsed[2] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	
#define VALIDATE
#if defined VALIDATE
        //validate
        for(int i=0;i<ne;i++) {
                FTYPE p1 = vout_ASpT_vec[i]; FTYPE p2 = vout_ASpT_scalar[i], p3 = vout_scalar[i];
                if(p1 < 0) p1 *= -1;
                if(p2 < 0) p2 *= -1;
				if(p3 < 0) p3 *= -1;
                FTYPE diff, diff2;
                diff = p1 - p2;
                if(diff < 0) diff *= -1;
				diff2 = p1 - p2;
                if(diff2 < 0) diff2 *= -1;
				diff+=diff2;
                if(MAX(MAX(p1,p2),p3) !=0 && diff / MAX(MAX(p1,p2),p3) > 0.01) {
					fprintf(stdout, "%d %f %f %f %s\n", i, p1, p2, p3, dir);
					return;
                }
				// fprintf(stdout, "%d %f %f %f\n", i, p1, p2, p3);

        }
        
#endif
        fprintf(stdout, " %f,%f,%f\n", elapsed[0]/(double)ITER, elapsed[1]/(double)ITER, elapsed[2]/(double)ITER);
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}


