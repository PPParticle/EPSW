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

#define MFACTOR (32)
#define LOG_MFACTOR (5)
#define BSIZE (1024 / 1)
#define BF (BSIZE / 32)
#define INIT_GRP (10000000)
#define INIT_LIST (-1)
#define THRESHOLD (16 * 1)
#define BH (128 * 1)
#define LOG_BH (7)
#define BW (128 * 1)
#define MIN_OCC (BW * 3 / 4)
#define SBSIZE (128)
#define SBF (SBSIZE / 32)
#define DBSIZE (1024)
#define DBF (DBSIZE / 32)
#define SPBSIZE (256)
#define SPBF (SPBSIZE / 32)
#define STHRESHOLD (1024 / 2 * 1)
#define SSTRIDE (STHRESHOLD / SPBF)
#define NTHREAD (68)
#define SC_SIZE (2048)

double vari, avg;
double avg0[NTHREAD];
struct v_struct *temp_v, *gold_temp_v;
int sc, nr, nc, ne, gold_ne, npanel, mne, mne_nr;
int nr0;

int *csr_v;
int *csr_e, *csr_e0;
FTYPE *csr_ev, *csr_ev0;
// int *mcsr_v;
int *mcsr_e; // can be short type
int *mcsr_cnt;
int *mcsr_list;
int *mcsr_chk;

int *baddr, *saddr;
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

	fprintf(stdout, "%s\n, ", argv[1]);
	sc = 128;
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
	for (i = 0; i < nr + 1; i++)
	{
		fscanf(fp, "%d ", &csr_v[i]);
	}
	for (i = 0; i < ne; i++)
	{
		fscanf(fp, "%d ", &csr_e0[i]);
		csr_ev0[i] = (FTYPE)(rand() % 1048576) / 1048576;
	}

	csr_e = new int[ne];
	csr_ev = new FTYPE[ne];
	fprintf(stdout, "%d,%d,%d\n", nr0, nc, ne);
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

	for (int row_panel = 0; row_panel < nr / BH; row_panel++)
	{
		for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++)
		{
			for (int j = csr_v[i]; j < csr_v[i + 1]; j++)
			{
				csr_e1[0][j] = csr_e0[j];
			}
		}
	}

	gettimeofday(&starttime, NULL);

	for (int row_panel = 0; row_panel < nr / BH; row_panel++)
	{
		int i, j, t_sum = 0;
		for (i = row_panel * BH; i < (row_panel + 1) * BH; i++)
		{
			memset(scr_pad_t, 0, sizeof(char) * SC_SIZE);
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

		// for(int l= row_panel*BH; l < (row_panel+1)*BH; l++){
		// 	printf("(%d, %d) ",csr_e1[flag][l], csr_e0[l]);
		// }
		// printf("\n"); // 数据集太小

		int weight = 1;

		int cr = 0;

		// dense bit extract (and mcsr_e making)
		for (i = csr_v[row_panel * BH] + 1; i < csr_v[(row_panel + 1) * BH]; i++)
		{
			if (csr_e1[flag][i - 1] == csr_e1[flag][i])
				weight++;
			else
			{
				if (weight >= THRESHOLD)
				{ // 列中数据大于 THRESHOLD 算为稠密列，并计数
					cr++;
				}
				weight = 1;
			}
		}
		if (weight >= THRESHOLD)
		{
			cr++;
		}
		mcsr_cnt[row_panel + 1] = CEIL(cr, BW) + 1; // BW 2D-tile, 记录每个panel中稠密列的个数
	}

	// printf("mcsr_cnt  mcsr_chk\n");
	// for(int i=0; i<npanel; i++){
	// 	printf("%d %d \n", mcsr_cnt[i], mcsr_chk[i]);
	// }

	for (int i = 1; i <= npanel; i++)
		mcsr_cnt[i] += mcsr_cnt[i - 1]; // mscr_cnt 记录每个panel中 2D-tile 的个数使用类 csr_v 的存储方式
	mcsr_e[BH * mcsr_cnt[npanel]] = ne; // mcsr_cnt[npanel] = 总 2D-tile 个数？

	// printf("mcsr_cnt  mcsr_chk\n");
	// for(int i=0; i<npanel; i++){
	// 	printf("%d %d \n", mcsr_cnt[i], mcsr_chk[i]);
	// }

	for(int i=0;i<npanel;i++){
		cout << "("<< mcsr_cnt[i] << ", " << mcsr_chk[i] << ") ";
	}
	cout << endl;

	int tid = 0;
	for (int row_panel = 0; row_panel < nr / BH && tid < NTHREAD; row_panel++, tid++)
	{
		if (mcsr_chk[row_panel] == 0)
		{ // 稠密的
			int i, j;
			int flag = 0;
			int cq = 0, cr = 0;
			int base = mcsr_cnt[row_panel] * BH;						 // 128*128 block
			int mfactor = mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]; // 2D-tile 个数
			int weight = 1;

			// mcsr_e making
			for (i = csr_v[row_panel * BH] + 1; i < csr_v[(row_panel + 1) * BH]; i++)
			{
				if (csr_e1[flag][i - 1] == csr_e1[flag][i])
					weight++;
				else
				{
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
					if (cr == BW)
					{
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
			int bidx = csr_v[row_panel * BH];
			int bseg = csr_v[(row_panel + 1) * BH] - bidx;
			memcpy(&csr_e[bidx], &csr_e0[bidx], sizeof(int) * bseg);
			memcpy(&csr_ev[bidx], &csr_ev0[bidx], sizeof(FTYPE) * bseg);
		}
	}

	for (int i = 0; i < NTHREAD; i++)
	{
		avg += avg0[i];
		// fprintf(stdout, "%5.2f ", avg0[i]);
	}

	avg /= (double)nr;

	for (int i = 0; i < nr; i++)
	{

		int idx = (mcsr_cnt[i >> LOG_BH]) * BH + (mcsr_cnt[(i >> LOG_BH) + 1] - mcsr_cnt[i >> LOG_BH]) * ((i & (BH - 1)) + 1);
		int diff = csr_v[i + 1] - mcsr_e[idx - 1];
		double r = ((double)diff - avg);
		vari += r * r;

		if (diff >= STHRESHOLD)
		{
			int pp = (diff) / STHRESHOLD;
			// printf("%d %d\n", pp, special_p);
			for (int j = 0; j < pp && special_p <= ne; j++)
			{
				special[special_p] = i;
				special2[special_p] = j * STHRESHOLD;
				special_p++;
			}
		}
	}
	// fprintf(stdout, "%5.2lf\n", vari);
	vari /= (double)nr;

	gettimeofday(&endtime, NULL);

	double elapsed0 = ((starttime.tv_sec - starttime0.tv_sec) * 1000000 + starttime.tv_usec - starttime0.tv_usec) / 1000000.0;
	p_elapsed = ((endtime.tv_sec - starttime.tv_sec) * 1000000 + endtime.tv_usec - starttime.tv_usec) / 1000000.0;
	// fprintf(stdout, "%f,%f\n", elapsed0*1000, p_elapsed*1000);

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
	// FILE *fpo = fopen("SpMM_KNL_SP.out", "a");
	// FILE *fpo2 = fopen("SpMM_KNL_SP_preprocessing.out", "a");

	double elapsed[3];
	FTYPE *vin, *vout, *vout_ASpT;
	vin = new FTYPE[nc * sc]{0};
	vout = new FTYPE[nr * sc]{0};
	vout_ASpT = new FTYPE[nr * sc]{0};

	// for(sc=8; sc<=128; sc*=4) {
	sc = 128;
	struct timeval starttime, endtime;
	for (int i = 0; i < nc * sc; i++)
	{
		vin[i] = (FTYPE)(rand() % 1048576) / 1048576;
	}
	// fin matrix * vin matrix  = vout matrix

	// #define ITER (128*1/128)
	// if(vari < 5000*1/1*1) {

	gettimeofday(&starttime, NULL);

	////begin
	// for(int loop=0;loop<ITER;loop++) {
	// printf("npanel = %d\n", npanel);
	// for(int i=0; i<npanel; i++){
	// 	printf("%d ", mcsr_cnt[i]);
	// }
	// printf("\n");

	for (int loop = 0; loop < 1; loop++)
	{
		for (int row_panel = 0; row_panel < nr / BH; row_panel++)
		{
			// dense
			int stride;
			for (stride = 0; stride < mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]; stride++)
			{
				for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++)
				{
					int dummy = mcsr_cnt[row_panel] * BH + (i & (BH - 1)) * (mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]) + stride; // i&(BH-1)  = I%BH 表示 i 在panel中的相对位置																								// mcsr_cnt 以 panel 为单位
					int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy + 1];

					int interm = loc1 + (((loc2 - loc1) >> 3) << 3);
					int j;

					for (j = loc1; j < interm; j += 8)
					{
						for (int k = 0; k < sc; k++)
						{
							vout[i * sc + k] = vout[i * sc + k] + csr_ev[j] * vin[csr_e[j] * sc + k] + csr_ev[j + 1] * vin[csr_e[j + 1] * sc + k] + csr_ev[j + 2] * vin[csr_e[j + 2] * sc + k] + csr_ev[j + 3] * vin[csr_e[j + 3] * sc + k] + csr_ev[j + 4] * vin[csr_e[j + 4] * sc + k] + csr_ev[j + 5] * vin[csr_e[j + 5] * sc + k] + csr_ev[j + 6] * vin[csr_e[j + 6] * sc + k] + csr_ev[j + 7] * vin[csr_e[j + 7] * sc + k];
						}
					}
					for (; j < loc2; j++)
					{
						for (int k = 0; k < sc; k++)
						{
							vout[i * sc + k] += csr_ev[j] * vin[csr_e[j] * sc + k];
						}
					}

					// vfloat32m2_t va, vb, vtmp;
					// vint32m2_t vidx;
					// vfloat32m1_t vans = vfmv_v_f_f32m1(0.0, 1);
					// size_t vl = vsetvlmax_e32m2();
					// vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());
					// fprintf(stdout, "rvv \n");
					// fflush(stdout);
					// for(j=loc1; j<interm; j+=8){
					// 	va = vle32_v_f32m2(csr_ev+j, vl);
					// 	fprintf(stdout, "va load \n");
					// 	fflush(stdout);
					// 	vidx = vle32_v_i32m2(csr_e+j, vl);
					// 	vidx = vmul_vx_i32m2(vidx, sc, vl);
					// 	vidx = vsll_vx_i32m2(vidx, 2, vl);
					// 	fprintf(stdout, "vidx load \n");
					// 	fflush(stdout);

					// 	for(int k=0; k<sc; k++){
					// 		vidx = vadd_vx_i32m2(vidx, k, vl);
					// 		vb = vluxei32_v_f32m2(vin, vidx, vl);
					// 		vtmp = vfmul_vv_f32m2(va, vb, vl);
					// 		vans = vfredusum_vs_f32m2_f32m1(vans, vtmp, vans, vl);
					// 		vout_rvv[i*sc+k] += vfmv_f_s_f32m1_f32 (vans);
					// 	}
					// }
					// for(; j<loc2; j++) {
					// 	for(int k=0; k<sc; k++) {
					// 		vout_rvv[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
					// 	}
					// }

					// for(long i=0; i<nc*sc; i++){
					// 	if(vout[i]!=vout_rvv[i]){
					// 		fprintf(stdout, "%d %d %d\n", i, vout[i], vout_rvv[i]);
					// 		break;
					// 	}
					// }
					// fflush(stdout);
				}
			}

			// sparse
			for (int i = row_panel * BH; i < (row_panel + 1) * BH; i++)
			{

				int dummy = mcsr_cnt[row_panel] * BH + (i & (BH - 1)) * (mcsr_cnt[row_panel + 1] - mcsr_cnt[row_panel]) + stride;
				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy + 1];

				// printf("(%d %d %d %d %d)\n", i, csr_v[i], loc1, csr_v[i+1], loc2);
				// printf("%d %d %d %d %d %d %d\n", i, dummy, stride, csr_v[i], loc1, csr_v[i+1], loc2);

				int interm = loc1 + (((loc2 - loc1) >> 3) << 3);
				int j;
				for (j = loc1; j < interm; j += 8)
				{
					for (int k = 0; k < sc; k++)
					{
						vout[i * sc + k] = vout[i * sc + k] + csr_ev[j] * vin[csr_e[j] * sc + k] + csr_ev[j + 1] * vin[csr_e[j + 1] * sc + k] + csr_ev[j + 2] * vin[csr_e[j + 2] * sc + k] + csr_ev[j + 3] * vin[csr_e[j + 3] * sc + k] + csr_ev[j + 4] * vin[csr_e[j + 4] * sc + k] + csr_ev[j + 5] * vin[csr_e[j + 5] * sc + k] + csr_ev[j + 6] * vin[csr_e[j + 6] * sc + k] + csr_ev[j + 7] * vin[csr_e[j + 7] * sc + k];
					}
				}
				for (; j < loc2; j++)
				{
					for (int k = 0; k < sc; k++)
					{
						vout[i * sc + k] += csr_ev[j] * vin[csr_e[j] * sc + k];
					}
				}

				// vfloat32m2_t va, vb, vtmp;
				// vint32m2_t vidx;
				// vfloat32m1_t vans = vfmv_v_f_f32m1(0.0, 1);
				// size_t vl = vsetvlmax_e32m2();
				// vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());

				// for(j=loc1; j<interm; j+=8){
				// 	va = vle32_v_f32m2(csr_ev+j, vl);
				// 	vidx = vle32_v_i32m2(csr_e+j, vl);
				// 	vidx = vmul_vx_i32m2(vidx, sc, vl);
				// 	vidx = vsll_vx_i32m2(vidx, 2, vl);
				// 	for(int k=0; k<sc; k++){
				// 		vidx = vadd_vx_i32m2(vidx, k, vl);
				// 		vb = vluxei32_v_f32m2(vin, vidx, vl);
				// 		vtmp = vfmul_vv_f32m2(va, vb, vl);
				// 		vans = vfredusum_vs_f32m2_f32m1(vans, vtmp, vans, vl);
				// 		vout_rvv[i*sc+k] += vfmv_f_s_f32m1_f32 (vans);
				// 	}
				// }
				// for(; j<loc2; j++) {
				// 	for(int k=0; k<sc; k++) {
				// 		vout_rvv[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
				// 	}
				// }
			}
		}
		////end
	}
	gettimeofday(&endtime, NULL);

	// } else { // big var

	// 	gettimeofday(&starttime, NULL);
	// ////begin
	// for(int loop=0;loop<ITER;loop++) {
	// 	for(int row_panel=0; row_panel<nr/BH; row_panel ++) {
	// 		//dense
	// 		int stride;
	// 		for(stride = 0; stride < mcsr_cnt[row_panel+1]-mcsr_cnt[row_panel]-1; stride++) {

	// 			for(int i=row_panel*BH; i<(row_panel+1)*BH; i++) {
	// 		        int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
	//  				int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];

	// 				int interm = loc1 + (((loc2 - loc1)>>3)<<3);
	// 				int j;
	// 				for(j=loc1; j<interm; j+=8) {
	// 					for(int k=0; k<sc; k++) {
	// 						vout[i*sc+k] = vout[i*sc+k] + csr_ev[j] * vin[csr_e[j]*sc + k]
	// 						+ csr_ev[j+1] * vin[csr_e[j+1]*sc + k]
	// 						+ csr_ev[j+2] * vin[csr_e[j+2]*sc + k]
	// 						+ csr_ev[j+3] * vin[csr_e[j+3]*sc + k]
	// 						+ csr_ev[j+4] * vin[csr_e[j+4]*sc + k]
	// 						+ csr_ev[j+5] * vin[csr_e[j+5]*sc + k]
	// 						+ csr_ev[j+6] * vin[csr_e[j+6]*sc + k]
	// 						+ csr_ev[j+7] * vin[csr_e[j+7]*sc + k];
	// 					}
	// 				}
	// 				for(; j<loc2; j++) {
	// 					for(int k=0; k<sc; k++) {
	// 						vout[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
	// 					}
	// 				}

	// 				vfloat32m2_t va, vb, vtmp;
	// 				vint32m2_t vidx;
	// 				vfloat32m2_t vans = vfmv_v_f_f32m2(0.0, 1);
	// 				size_t vl = vsetvlmax_e32m2();
	// 				vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());

	// 				for(j=loc1; j<interm; j+=8){
	// 					va = vle32_v_f32m2(csr_ev+j, vl);
	// 					vidx = vle32_v_i32m2(csr_e+j, vl);
	// 					vidx = vmul_vx_i32m2(vidx, sc, vl);
	// 					vidx = vsll_vx_i32m2(vidx, 2, vl);
	// 					for(int k=0; k<sc; k++){
	// 						vidx = vadd_vx_i32m2(vidx, k, vl);
	// 						vb = vluxei32_v_f32m2(vin, vidx, vl);
	// 						vtmp = vfmul_vv_f32m2(va, vb, vl);
	// 						vans = vfredusum_vs_f32m2_f32m2(vans, vtmp, vans, vl);
	// 						vout_rvv[i*sc+k] += vfmv_f_s_f32m2_f32 (vans);
	// 					}
	// 				}
	// 				for(; j<loc2; j++) {
	// 					for(int k=0; k<sc; k++) {
	// 						vout_rvv[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
	// 					}
	// 				}
	// 			}

	// 		}
	// 		//sparse
	// 		for(int i=row_panel*BH; i<(row_panel+1)*BH; i++) {

	// 			int dummy = mcsr_cnt[row_panel]*BH + (i&(BH-1))*(mcsr_cnt[row_panel+1] - mcsr_cnt[row_panel]) + stride;
	//  			int loc1 = mcsr_e[dummy], loc2 = mcsr_e[dummy+1];

	// 			loc1 += ((loc2 - loc1)/STHRESHOLD)*STHRESHOLD;

	// 			int interm = loc1 + (((loc2 - loc1)>>3)<<3);
	// 			int j;
	// 			for(j=loc1; j<interm; j+=8) {
	// 				for(int k=0; k<sc; k++) {
	// 					vout[i*sc+k] = vout[i*sc+k] + csr_ev[j] * vin[csr_e[j]*sc + k]
	// 					+ csr_ev[j+1] * vin[csr_e[j+1]*sc + k]
	// 					+ csr_ev[j+2] * vin[csr_e[j+2]*sc + k]
	// 					+ csr_ev[j+3] * vin[csr_e[j+3]*sc + k]
	// 					+ csr_ev[j+4] * vin[csr_e[j+4]*sc + k]
	// 					+ csr_ev[j+5] * vin[csr_e[j+5]*sc + k]
	// 					+ csr_ev[j+6] * vin[csr_e[j+6]*sc + k]
	// 					+ csr_ev[j+7] * vin[csr_e[j+7]*sc + k];
	// 				}
	// 			}
	// 			for(; j<loc2; j++) {
	// 				for(int k=0; k<sc; k++) {
	// 					vout[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
	// 				}
	// 			}

	// 				vfloat32m2_t va, vb, vtmp;
	// 				vint32m2_t vidx;
	// 				vfloat32m2_t vans = vfmv_v_f_f32m2(0.0, 1);
	// 				size_t vl = vsetvlmax_e32m2();
	// 				vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());

	// 				for(j=loc1; j<interm; j+=8){
	// 					va = vle32_v_f32m2(csr_ev+j, vl);
	// 					vidx = vle32_v_i32m2(csr_e+j, vl);
	// 					vidx = vmul_vx_i32m2(vidx, sc, vl);
	// 					vidx = vsll_vx_i32m2(vidx, 2, vl);
	// 					for(int k=0; k<sc; k++){
	// 						vidx = vadd_vx_i32m2(vidx, k, vl);
	// 						vb = vluxei32_v_f32m2(vin, vidx, vl);
	// 						vtmp = vfmul_vv_f32m2(va, vb, vl);
	// 						vans = vfredusum_vs_f32m2_f32m2(vans, vtmp, vans, vl);
	// 						vout_rvv[i*sc+k] += vfmv_f_s_f32m2_f32 (vans);
	// 					}
	// 				}
	// 				for(; j<loc2; j++) {
	// 					for(int k=0; k<sc; k++) {
	// 						vout_rvv[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
	// 					}
	// 				}
	// 		}
	// 	}
	// 	for(int row_panel=0; row_panel<special_p;row_panel ++) {
	// 		int i=special[row_panel];

	// 		int dummy = mcsr_cnt[i>>LOG_BH]*BH + ((i&(BH-1))+1)*(mcsr_cnt[(i>>LOG_BH)+1] - mcsr_cnt[i>>LOG_BH]);

	// 		int loc1 = mcsr_e[dummy-1] + special2[row_panel];
	// 		int loc2 = loc1 + STHRESHOLD;

	// 		//int interm = loc1 + (((loc2 - loc1)>>3)<<3);
	// 		int j;
	// //assume to 128
	// 		FTYPE temp_r[128]={0,};
	// 		//for(int e=0;e<128;e++) {
	// 		//	temp_r[e] = 0.0f;
	// 		//}

	// 		for(j=loc1; j<loc2; j+=8) {
	// 			for(int k=0; k<sc; k++) {
	// 				temp_r[k] = temp_r[k] + csr_ev[j] * vin[csr_e[j]*sc + k]
	// 				+ csr_ev[j+1] * vin[csr_e[j+1]*sc + k]
	// 				+ csr_ev[j+2] * vin[csr_e[j+2]*sc + k]
	// 				+ csr_ev[j+3] * vin[csr_e[j+3]*sc + k]
	// 				+ csr_ev[j+4] * vin[csr_e[j+4]*sc + k]
	// 				+ csr_ev[j+5] * vin[csr_e[j+5]*sc + k]
	// 				+ csr_ev[j+6] * vin[csr_e[j+6]*sc + k]
	// 				+ csr_ev[j+7] * vin[csr_e[j+7]*sc + k];
	// 			}
	// 		}
	// 		for(int k=0; k<sc; k++) {
	// 			vout[i*sc+k] += temp_r[k];
	// 		}

	// 		vfloat32m2_t va, vb, vtmp;
	// 				vint32m2_t vidx;
	// 				vfloat32m2_t vans = vfmv_v_f_f32m2(0.0, 1);
	// 				size_t vl = vsetvlmax_e32m2();
	// 				vtmp = vfmv_v_f_f32m2(0.0, vsetvlmax_e32m2());

	// 				for(j=loc1; j<interm; j+=8){
	// 					va = vle32_v_f32m2(csr_ev+j, vl);
	// 					vidx = vle32_v_i32m2(csr_e+j, vl);
	// 					vidx = vmul_vx_i32m2(vidx, sc, vl);
	// 					vidx = vsll_vx_i32m2(vidx, 2, vl);
	// 					for(int k=0; k<sc; k++){
	// 						vidx = vadd_vx_i32m2(vidx, k, vl);
	// 						vb = vluxei32_v_f32m2(vin, vidx, vl);
	// 						vtmp = vfmul_vv_f32m2(va, vb, vl);
	// 						vans = vfredusum_vs_f32m2_f32m2(vans, vtmp, vans, vl);
	// 						vout_rvv[i*sc+k] += vfmv_f_s_f32m2_f32 (vans);
	// 					}
	// 				}
	// 				for(; j<loc2; j++) {
	// 					for(int k=0; k<sc; k++) {
	// 						vout_rvv[i*sc + k] += csr_ev[j] * vin[csr_e[j]*sc + k];
	// 					}
	// 				}
	// 	}
	// } // end loop
	// 	gettimeofday(&endtime, NULL);

	// }

	// double elapsed = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;

	// if(sc == 8) elapsed[0] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	// else if(sc == 32) elapsed[1] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;
	// else elapsed[2] = ((endtime.tv_sec-starttime.tv_sec)*1000000 + endtime.tv_usec-starttime.tv_usec)/1000000.0;

	//  }

#define VALIDATE
#if defined VALIDATE
	// validate
	for (int i = 0; i < nr * sc; i++)
	{
		vout[i] = 0.0f;
	}

	for (int i = 0; i < nr; i++)
	{
		for (int j = 0; j < sc; j++)
		{
			for (int k = csr_v[i]; k < csr_v[i + 1]; k++)
			{
				vout[i * sc + j] += csr_ev[k] * vin[csr_e[k] * sc + j];
			}
		}
	}

	int num_diff = 0;
	for (int i = 0; i < nr * sc; i++)
	{
		FTYPE p1 = vout[i];
		FTYPE p2 = vout_ASpT[i];

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
			// if(num_diff < 20*1*1) fprintf(stdout, "%d %f %f\n", i, vout[i], vout[i]*ITER);
			// if(vout[i] < vout[i]) fprintf(stdout, "%d %f %f\n", i, vout[i], vout[i]);

			num_diff++;
		}
		// fprintf(stdout, "%5.2f %5.2f\n", p1, p2);
	}
	//      fprintf(stdout, "num_diff : %d\n", num_diff);
	fprintf(stdout, "diff : %f\n", (double)num_diff / (nr * sc) * 100);
//      fprintf(stdout, "ne : %d\n", gold_ne);
#endif
	// fprintf(stdout, "%f,%f\n", (double)elapsed*1000/ITER, (double)ne*2*sc*ITER/elapsed/1000000000);
	//  fprintf(stdout, "%f,", (double)ne*2*8*ITER/elapsed[0]/1000000000);
	//  fprintf(stdout, "%f,%f,", (double)ne*2*32*ITER/elapsed[1]/1000000000, (double)ne*2*128*ITER/elapsed[2]/1000000000);
	//  fprintf(stdout, "%f\n", p_elapsed/(elapsed[2]/ITER));
	fprintf(stdout, "%f,", (double)ne * 2 * 8 * 1 / elapsed[0] / 1000000000);
	fprintf(stdout, "%f,%f,", (double)ne * 2 * 32 * 1 / elapsed[1] / 1000000000, (double)ne * 2 * 128 * 1 / elapsed[2] / 1000000000);
	fprintf(stdout, "%f\n", p_elapsed / (elapsed[2] / 1));

	// fprintf(fpo, "%f,%f,", (double)ne*2*32*ITER/elapsed[1]/1000000000, (double)ne*2*128*ITER/elapsed[2]/1000000000);
	// fprintf(fpo, "%f,", (double)num_diff/(nr*sc)*100);
	// fprintf(fpo2, "%f", p_elapsed/(elapsed[2]/ITER));
	// fclose(fpo); fclose(fpo2);
}

int main(int argc, char **argv)
{
	ready(argc, argv);
	gen();
	mprocess();
}
