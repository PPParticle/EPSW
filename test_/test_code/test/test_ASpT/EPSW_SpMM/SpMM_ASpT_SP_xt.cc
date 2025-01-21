#include <stdio.h>
#include<iostream>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <math.h>
#include <riscv_vector.h>
using namespace std;

#define ERR fprintf(stderr, "ERR\n");

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CEIL(a, b) (((a) + (b)-1) / (b))
#define FTYPE float
#define THRESHOLD (64)
#define BH (128 * 1)
#define BW (128 * 1)
#define MIN_OCC (BW * 3 / 4)
#define NTHREAD (68)
#define SC_SIZE (1024)

double avg;
double avg0[NTHREAD];
int sc, nr, nc, ne, npanel;
int nr0;

int *csr_v;
int *csr_e, *csr_e0;
FTYPE *csr_ev, *csr_ev0;

int *mcsr_e; // can be short type
int *mcsr_cnt;
int *mcsr_list;
int *mcsr_chk;


int *special;
int *special2;
int special_p;
int scr_pad_t[SC_SIZE]{0};

void ready(int argc, char **argv)
{
	FILE *fp;
	int i;

	// fprintf(stdout, "%s,\n ", argv[1]);
	sc = atoi(argv[2]);
	fp = fopen(argv[1], "r");
	fscanf(fp, "%d, %d, %d", &nr, &nc, &ne);

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
		csr_ev0[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// csr_ev0[i] = i;
	}

	csr_e = new int[ne];
	csr_ev = new FTYPE[ne];
	fprintf(stdout, "%d,%d,%d,%d ", nr0, nc, ne, sc);
}

void gen()
{
	special = new int[ne]{0};
	special2 = new int[ne]{0};
	mcsr_cnt = new int[npanel + 1]{0};
	mcsr_chk = new int[npanel + 1]{0};
	mcsr_e = new int[ne]{0};

	int **csr_e1 = new int *[2];
	short **coo = new short *[2];
	for (int i = 0; i < 2; i++)
	{
		csr_e1[i] = new int[ne]{0};
		coo[i] = new short[ne]{0};
	}

	int bv_size = CEIL(nc, 32);
	unsigned int **bv = new unsigned int *[NTHREAD];
	for (int i = 0; i < NTHREAD; i++)
		bv[i] = new unsigned int[bv_size]{0};

	

	struct timeval starttime0;
	struct timeval starttime, endtime;

	for (int row_panel = 0; row_panel < nr / BH; row_panel++){
		for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++){
			for (int j = csr_v[i]; j < csr_v[i + 1]; j++){
				csr_e1[0][j] = csr_e0[j];
			}
		}
	}
	int cnt = 0;

	for (int row_panel = 0; row_panel < nr / BH; row_panel++)
	{
		int i, j, t_sum = 0;
		memset(scr_pad_t, 0, sizeof(char) * SC_SIZE);
		for (i = row_panel * BH; i < (row_panel + 1) * BH; i++)
		{
			for (j = csr_v[i]; j < csr_v[i + 1]; j++)
			{
				coo[0][j] = (i & (BH - 1));			 // coo[0] 按组存储所有元素的行坐标 0-127
				int k = (csr_e0[j] & (SC_SIZE - 1)); // 保证元素的列坐标不会越界 (SC_SIZE) 2^11
				// 统计同一列元素个数，如果一列元素个数总和超过 THRESHOLD t_sum++
				if (scr_pad_t[k] < THRESHOLD)
				{
					if (scr_pad_t[k] == THRESHOLD - 1)
						t_sum++;
					scr_pad_t[k]++;
				}
			}
		}

		// printf("%d ", t_sum);
		// 如果 panel 太过稀疏就跳过对 panel 排序的过程，可能认定为稀疏 panel
		if (t_sum < MIN_OCC)
		{ // MIN_OCC = 96
			mcsr_chk[row_panel] = 1;
			mcsr_cnt[row_panel + 1] = 1;
			continue;
		}

		// 将 panel 中每个元素按照列元素的索引大小进行重新排列，coo 记录行坐标，csr_e1 记录列坐标
		int flag = 0;
		for (int stride = 1; stride <= BH / 2; stride *= 2, flag = 1 - flag)
		{
			for (int pivot = row_panel * BH; pivot < (row_panel + 1) * BH; pivot += stride * 2)
			{
				int l1, l2;
				for (i = l1 = csr_v[pivot], l2 = csr_v[pivot + stride]; l1 < csr_v[pivot + stride] && l2 < csr_v[pivot + stride * 2]; i++)
				{
					if (csr_e1[flag][l1] <= csr_e1[flag][l2])
					{
						coo[1 - flag][i] = coo[flag][l1];		  // 第 i 个元素的 row idx
						csr_e1[1 - flag][i] = csr_e1[flag][l1++]; // 第 i 个元素的 col idx
					}
					else
					{
						coo[1 - flag][i] = coo[flag][l2];
						csr_e1[1 - flag][i] = csr_e1[flag][l2++];
					}
				}
				while (l1 < csr_v[pivot + stride])
				{
					coo[1 - flag][i] = coo[flag][l1];
					csr_e1[1 - flag][i++] = csr_e1[flag][l1++];
				}
				while (l2 < csr_v[pivot + stride * 2])
				{
					coo[1 - flag][i] = coo[flag][l2];
					csr_e1[1 - flag][i++] = csr_e1[flag][l2++];
				}
			}
		}

		int weight = 1;
		int cr = 0;

		// dense bit extract (and mcsr_e making)
		for (i = csr_v[row_panel * BH] + 1; i < csr_v[(row_panel + 1) * BH]; i++){
			if (csr_e1[flag][i-1] == csr_e1[flag][i])
				weight++;
			else{
				if (weight >= THRESHOLD)
				{ // 列中数据大于 THRESHOLD 算为稠密列，并计数
					cr++;
				}
				weight = 1;
			}
		}
		if (weight >= THRESHOLD){
			cr++;
		}
		mcsr_cnt[row_panel + 1] = CEIL(cr, BW) + 1; // BW 2D-tile, 记录每个panel中稠密列的个数
	}

	for (int i = 1; i <= npanel; i++)
		mcsr_cnt[i] += mcsr_cnt[i - 1]; // mscr_cnt 记录每个panel中 2D-tile 的个数使用类 csr_v 的存储方式
	mcsr_e[BH * mcsr_cnt[npanel]] = ne; // mcsr_cnt[npanel] = 总 2D-tile 个数？

	// printf("mcsr_cnt  mcsr_chk\n");
	// for(int i=0; i<npanel+1; i++){
	// 	printf("    %d        %d \n", mcsr_cnt[i], mcsr_chk[i]);
	// }
	
	int tid = 0;
	for (int row_panel = 0; row_panel < nr / BH && tid < NTHREAD; row_panel++, tid++){
		if (mcsr_chk[row_panel] == 0)
		{ // 稠密的
			cnt++;
			int i, j;
			int flag = 0;
			int cq = 0, cr = 0;
			int base = mcsr_cnt[row_panel] * BH;						 // 128*128 block
			int mfactor = mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]; // 2D-tile 个数
			int weight = 1;

			// mcsr_e making
			for (i = csr_v[row_panel * BH] + 1; i < csr_v[(row_panel + 1) * BH]; i++){
				if (csr_e1[flag][i - 1] == csr_e1[flag][i])
					weight++;
				else{
					int reminder = (csr_e1[flag][i - 1] & (THRESHOLD-1));
					if (weight >= THRESHOLD){
						cr++;
						bv[tid][csr_e1[flag][i - 1] >> 5] |= (1 << reminder);
						for (j = i - weight; j <= i - 1; j++){
							mcsr_e[base + coo[flag][j] * mfactor + cq + 1]++;
						}
					}
					else{
						bv[tid][csr_e1[flag][i - 1] >> 5] &= (0xFFFFFFFF - (1 << reminder));
					}
					if (cr == BW){
						cq++;
						cr = 0;
					}
					weight = 1;
				}
			}
			int reminder = (csr_e1[flag][i - 1] & 31);
			if (weight >= THRESHOLD)
			{
				cr++;
				bv[tid][csr_e1[flag][i - 1] >> 5] |= (1 << reminder);
				for (j = i - weight; j <= i - 1; j++)
				{
					mcsr_e[base + coo[flag][j] * mfactor + cq + 1]++;
				}
			}
			else
			{
				bv[tid][csr_e1[flag][i - 1] >> 5] &= (0xFFFFFFFF - (1 << reminder));
			}

			// reordering
			int delta = mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel];
			int base0 = mcsr_cnt[row_panel] * BH;
			for (i = row_panel * BH; i < (row_panel + 1) * BH; i++)
			{
				int base = base0 + (i - row_panel * BH) * delta;
				int dpnt = mcsr_e[base] = csr_v[i];
				for (int j = 1; j < delta; j++)
				{
					mcsr_e[base + j] += mcsr_e[base + j - 1];
				}
				int spnt = mcsr_e[base0 + delta * (i - row_panel * BH + 1) - 1];
				avg0[tid] += csr_v[i + 1] - spnt;
				for (j = csr_v[i]; j < csr_v[i + 1]; j++)
				{
					int k = csr_e0[j];
					if ((bv[tid][k >> 5] & (1 << (k & 31))))
					{
						csr_e[dpnt] = csr_e0[j];
						csr_ev[dpnt++] = csr_ev0[j];
					}
					else
					{
						csr_e[spnt] = csr_e0[j];
						csr_ev[spnt++] = csr_ev0[j];
					}
				}
			}
		}
		else
		{
			int base0 = mcsr_cnt[row_panel] * BH;
			memcpy(&mcsr_e[base0], &csr_v[row_panel * BH], sizeof(int) * BH);
			avg0[tid] += csr_v[(row_panel + 1) * BH] - csr_v[row_panel * BH];
			int bidx = csr_v[row_panel * BH]; // 0
			int bseg = csr_v[(row_panel + 1) * BH] - bidx;
			memcpy(&csr_e[bidx], &csr_e0[bidx], sizeof(int) * bseg);
			memcpy(&csr_ev[bidx], &csr_ev0[bidx], sizeof(FTYPE) * bseg);
		}
	}
	
	for (int i = 0; i < NTHREAD; i++)
		free(bv[i]);
	for (int i = 0; i < 2; i++)
	{
		free(csr_e1[i]);
		free(coo[i]);
	}
	free(bv);
	free(csr_e1);
	free(coo);
}

void mprocess()
{
	double elapsed[3];
	FTYPE *vin, *vout, *vout_ASpT_vec, *vout_ASpT_sca;
	vin = new FTYPE[nc * sc]{0};
	vout = new FTYPE[nr * nc]{0};
	vout_ASpT_vec = new FTYPE[nr * nc]{0};
	vout_ASpT_sca = new FTYPE[nr * nc]{0};

	struct timeval starttime, endtime, 
	starttime_ASpT_sca, endtime_ASpT_sca,
	starttime_ASpT_vec, endtime_ASpT_vec;
	for (int i = 0; i < nc * sc; i++)
	{
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
		// vin[i] = i;
	}

	#define ITER 1
	
	gettimeofday(&starttime_ASpT_vec, NULL);
	for (int loop = 0; loop < ITER; loop++){
		for (int row_panel=0; row_panel<nr/BH; row_panel++){
			// dense
			int stride;
			for(stride = 0; stride < mcsr_cnt[row_panel+1]-mcsr_cnt[row_panel]-1; stride++){
				for(int i=row_panel*BH; i<(row_panel+1)*BH; i++){	
					int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 					int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
					int interm = loc1 + (((loc2 - loc1)>>3)<<3);
					int j;

					vfloat32m2_t va, vb, vtmp;
					vuint32m2_t vidx, vidx_t;
					vfloat32m1_t vans;
					size_t vl = vsetvlmax_e32m2();
					for(j=loc1; j<interm; j+=8){
						va = vle32_v_f32m2(csr_ev+j, vl);
						vidx = vle32_v_u32m2((unsigned int*)csr_e+j, vl);
						vidx = vmul_vx_u32m2(vidx, sc, vl);
						for(int k=0; k<sc; k++){
							vans = vfmv_v_f_f32m1(0.0, vl);
							vtmp = vfmv_v_f_f32m2(0.0, vl);
							vidx_t = vadd_vx_u32m2(vidx, k, vl);
							vidx_t = vsll_vx_u32m2(vidx_t, 2, vl);
							vb = vloxei32_v_f32m2(vin, vidx_t, vl);
							vtmp = vfmul_vv_f32m2(va, vb, vl);
							vans = vfredusum_vs_f32m2_f32m1(vans, vtmp, vans, vl);
							vout_ASpT_vec[i*sc+k] += vfmv_f_s_f32m1_f32 (vans);
						}
					}
					
					for(; j<loc2; j++) {
						for(int k=0; k<sc; k++) {
							vout_ASpT_vec[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
						}
					}

				}
			}

			// sparse
			for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++){
				int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
 				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];
				int interm = loc1 + (((loc2 - loc1)>>3)<<3);
				int j;

				vfloat32m2_t va, vb, vtmp;
				vuint32m2_t vidx, vidx_t;
				vfloat32m1_t vans;
				size_t vl = vsetvlmax_e32m2();
				for(j=loc1; j<interm; j+=8){
					va = vle32_v_f32m2(csr_ev+j, vl);
					vidx = vle32_v_u32m2((unsigned int*)csr_e+j, vl);
					vidx = vmul_vx_u32m2(vidx, sc, vl);
					for(int k=0; k<sc; k++){
						vans = vfmv_v_f_f32m1(0.0, 1);
						vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());
				    	vidx_t = vadd_vx_u32m2(vidx, k, vl);
						vidx_t = vsll_vx_u32m2(vidx_t, 2, vl);
						vb = vloxei32_v_f32m2(vin, vidx_t, vl);
						vtmp = vfmul_vv_f32m2(va, vb, vl);
						vans = vfredusum_vs_f32m2_f32m1(vans, vtmp, vans, vl);
						vout_ASpT_vec[i*sc+k] += vfmv_f_s_f32m1_f32 (vans);
					}
				}
				for(; j<loc2; j++) {
					for(int k=0; k<sc; k++) {
						vout_ASpT_vec[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
					}
				}
			}
		}
		////end
	}
	gettimeofday(&endtime_ASpT_vec, NULL);

	gettimeofday(&starttime_ASpT_sca, NULL);
	for (int loop = 0; loop < ITER; loop++){
		for (int row_panel=0; row_panel<nr/BH; row_panel++){
			// dense
			int stride;
			for (stride=0; stride<mcsr_cnt[row_panel+1]-mcsr_cnt[row_panel]-1; stride++){
				for (int i=row_panel*BH; i<(row_panel+1)*BH; i++){
					int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride; // i&(BH-1)  = I%BH 表示 i 在panel中的相对位置																								// mcsr_cnt 以 panel 为单位
					int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy + 1];
					// cout << loc1 << " " << loc2 << endl;
					int interm = loc1 + (((loc2 - loc1) >> 3) << 3);
					// cout << interm << endl;
					int j;
					for(j=loc1; j<loc2; j++) 
						for(int k=0; k<sc; k++) 
							vout_ASpT_sca[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
				}
			}

			// sparse
			for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++){
				int dummy = mcsr_cnt[row_panel] * BH + (i & (BH - 1)) * (mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]) + stride;
				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy + 1];
				int interm = loc1 + (((loc2 - loc1) >> 3) << 3);
				int j;
				for(j=loc1; j<loc2; j++) 
					for(int k=0; k<sc; k++) 
						vout_ASpT_sca[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
			}
		}
	}
	gettimeofday(&endtime_ASpT_sca, NULL);

	gettimeofday(&starttime, NULL);	
	for (int loop = 0; loop < ITER; loop++)
	{
		for (int i = 0; i < nr; i++){
			for (int j = csr_v[i]; j < csr_v[i + 1]; j++){
				for (int k = 0; k < sc; k++){
					vout[i*sc + k] += csr_ev[j] * vin[csr_e[j] * sc + k];
			   	}
			}
		}
	}
	gettimeofday(&endtime, NULL);

	elapsed[0] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_ASpT_sca.tv_sec-starttime_ASpT_sca.tv_sec)*1000000 + endtime_ASpT_sca.tv_usec-starttime_ASpT_sca.tv_usec)/1000000.0;
	elapsed[2] = ((endtime_ASpT_vec.tv_sec-starttime_ASpT_vec.tv_sec)*1000000 + endtime_ASpT_vec.tv_usec-starttime_ASpT_vec.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	for (int i = 0; i < nr * sc; i++)
	{
		FTYPE p1 = vout[i];
		FTYPE p2 = vout_ASpT_sca[i];
		FTYPE p3 = vout_ASpT_vec[i];

		if (p1 < 0)
			p1 *= -1;
		if (p2 < 0)
			p2 *= -1;
		if (p3 < 0)
			p3 *= -1;
		FTYPE diff1, diff2, diff3;
		diff1 = p1 - p2;
		diff2 = p1 - p2;
		diff3 = p2 - p3;
		if (diff1 < 0)
			diff1 *= -1;
		if (diff2 < 0)
			diff2 *= -1;
		if (diff3 < 0)
			diff3 *= -1;
		diff1+=diff2+diff3;
		if (diff1 > 0.01)
		{	
			cout << i << " " << p1 << " " << p2 << " " <<p3<<endl;
			break;
		}
	}
#endif
	fprintf(stdout, "%f,%f,%f \n", elapsed[0]/(double)ITER , elapsed[1]/(double)ITER, elapsed[2]/(double)ITER );
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
