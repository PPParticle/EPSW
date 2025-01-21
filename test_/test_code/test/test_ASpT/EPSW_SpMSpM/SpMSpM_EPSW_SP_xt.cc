#include <stdio.h>
#include <sys/time.h>
#include <cstring>
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
ELLPACK_P M2;
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
		csr_ev[i] = (FTYPE)(rand() & 1048575) / 1048576;
		// csr_ev[i] = i;
	}

	for (i = 0; i < nr_r+1; i++){
		fscanf(fp_r, "%d ", &csr_vr[i]);
	}

	for (i = 0; i < ne_r; i++){
		fscanf(fp_r, "%d", &csr_er[i]);
	}

	fclose(fp_l);
	fclose(fp_r);
	fprintf(stdout, "%d %d %d ", ne, nr, nc);
}
void gen(){

	vector<vector<FTYPE>> csr_evt(nc_l);
	for(int i=0; i<nc_l; i++){
		vector<FTYPE> tmp;
		csr_evt.push_back(tmp);
	}

	for(int i=0; i<nr_l; i++){
		for(int j=csr_v[i]; j<csr_v[i+1]; j++){
			csr_evt[csr_e[j]].push_back(csr_ev[j]);
		}
	}

	int cnt=0;
	for(int i=0; i<nc_l+1; i++){
		for(int j=0; j<csr_evt[i].size(); j++){
			// cout << cnt << " " << csr_evt[i][j] << endl;  正确
			csr_evr[cnt++] = csr_evt[i][j];
			
		}
	}
	M1.ELLPACKinitWithCSR(csr_v, csr_e, csr_ev, nr_l, nc_l, ne_l);
	M2.ELLPACKinitWithCSR(csr_vr, csr_er, csr_evr, nr_r, nc_r, ne_r);
	for(int i=0; i<ne_r; i++){
		col = MAX(col, csr_er[i]);
	}
	M3._rows = M1._rows;
	M3._cols = M2._cols;

	M3._rowidx=new unsigned int [M3._rows];
	M3._nnz_num=new unsigned int [M3._rows];
	M3._nnz_old_num=new unsigned int [M3._rows];
	M3._rowPtr=new FTYPE*[M3._rows];
	M3._coordPtr=new unsigned int*[M3._rows];

	for(unsigned int i=0; i<M2._rows; i++){
        col = MAX(col, M2._coordPtr[i][M2._nnz_num[i]-1]);
    }
	col++;

    for(int i=0; i<M3._rows; i++){
		M3._rowidx[i] = M1._rowidx[i];
        M3._nnz_num[i] = M3._nnz_old_num[i] = col;
        M3._rowPtr[i] = new Val [col]{0.0f};
        M3._coordPtr[i] = new unsigned int [col]{0};
    }
}

void mprocess()
{	
	double elapsed[2];
	FTYPE *vout_scalar = new FTYPE[nr*nc]{0};	
	FTYPE *vout_EPSW = new FTYPE[nr*nc]{0.0};

	struct timeval starttime_scalar, endtime_scalar, starttime_EPSW, endtime_EPSW;
	#define ITER 10

	////begin
	gettimeofday(&starttime_EPSW, NULL);

for(int loop=0; loop<ITER; loop++) {
	size_t vl;
    vfloat32m8_t vb, vc;
    vuint32m8_t vbidx, vcidx, vo; 
	vbool4_t vmask;


    for(unsigned int i=0; i<M1._rows; i++){        
        for(unsigned int j=0; j<M1._nnz_num[i]; j++){
            
            Val a_v = M1._rowPtr[i][j];
			unsigned int a_col = M1._coordPtr[i][j];
	        unsigned int a_row = M1._rowidx[i];
		
            for(unsigned int l=0; l<M2._nnz_num[a_col]; l+=vl){
				vc = vfmv_v_f_f32m8(0.0, vsetvlmax_e32m8());
                vl = vsetvl_e32m8(M2._nnz_num[a_col]-l);
                vb = vle32_v_f32m8(M2._rowPtr[a_col]+l, vl);
                vbidx = vle32_v_u32m8(M2._coordPtr[a_col]+l, vl);
                vcidx = vmv_v_v_u32m8(vbidx, vl);
                vcidx = vsll_vx_u32m8(vbidx, 2, vl);
                vc = vloxei32_v_f32m8(M3._rowPtr[i], vcidx, vl);			
                vc = vfmacc_vf_f32m8(vc, a_v, vb, vl);  
                vsuxei32_v_u32m8(M3._coordPtr[i], vcidx, vbidx, vl); // base, index, value, vl
				vsuxei32_v_f32m8(M3._rowPtr[i], vcidx, vc, vl);
            }
        }
    }
}
	gettimeofday(&endtime_EPSW, NULL);

	gettimeofday(&starttime_scalar, NULL);
	for (int loop = 0; loop < ITER; loop++){
		for (int i = 0; i < nr; i++){
				for(int j=0; j<nc; j++){
					int ni = csr_v[i];
					for(int nj=csr_v[j]; nj<csr_v[j+1]; ){
						if(ni>=csr_v[i+1]) {
							// cout << "b" <<ni << " " << csr_v[i+1]<< endl;
							break;
						}
						// cout << i << " " << j << " " << ni << " " << nj << " "<< csr_e[ni] << " " << csr_e[nj] << endl;
							
						if(csr_e[ni]==csr_e[nj]){
							vout_scalar[i*nc+j] += csr_ev[ni]*csr_ev[nj];
							ni++;
							nj++;
						}
						else if(csr_e[ni]>csr_e[nj]){
							nj++;
						}
							else ni++;
					}
				}
		}
	}
	
	gettimeofday(&endtime_scalar, NULL);
	
	elapsed[0] = ((endtime_EPSW.tv_sec-starttime_EPSW.tv_sec)*1000000 + endtime_EPSW.tv_usec-starttime_EPSW.tv_usec)/1000000.0;
	elapsed[1] = ((endtime_scalar.tv_sec-starttime_scalar.tv_sec)*1000000 + endtime_scalar.tv_usec-starttime_scalar.tv_usec)/1000000.0;

#define VALIDATE
#if defined VALIDATE
	// M3.print_coordPtr();
	// M3.print_value();

	for(unsigned int i=0; i<M3._rows; i++){
        unsigned int nnz=0;
        for(unsigned int j=0; j<col; j++){
            if(M3._rowPtr[i][j]!=0.0){
				// cout<< M3._rowPtr[i][j] << " ";
                nnz++;
            }
        }
		// cout << endl;
        M3._nnz_old_num[i] = nnz;
        M3._nnz_num[i] = M3._nnz_old_num[i];

        int tmp=0;
        FTYPE *tmp_v = new FTYPE [nnz]{0.0f};
        unsigned int *tmp_col = new unsigned int [nnz]{0};
		
        for(unsigned int j=0; j<col; j++){
            if(M3._rowPtr[i][j]!=0){
                tmp_v[tmp]=M3._rowPtr[i][j]; 
                tmp_col[tmp++]=M3._coordPtr[i][j]; 
            }
        } 

		delete [] M3._rowPtr[i];
        delete [] M3._coordPtr[i];
        M3._rowPtr[i] = tmp_v;
        M3._coordPtr[i] = tmp_col;
	}
	
	int *sca_csr_v, *sca_csr_e;
	FTYPE *sca_csr_ev;

	sca_csr_v = new int[nr+1];
    int cnt=0;
    for(int i=0; i<nr; i++){
        for(int j=0; j<nc; j++){
            if(vout_scalar[i*nc+j]!=0){
                cnt++;
            }
        }
        sca_csr_v[i+1] = cnt;
    }
    sca_csr_e = new int [cnt];
    sca_csr_ev = new Val [cnt];
	// cout << cnt << endl;
	int cnts=0;
    for(int i=0; i<nr; i++){
        for(int j=0; j<nc; j++){
            if(vout_scalar[i*nc+j]!=0){
                sca_csr_e[cnts] = j;
				// cout << sca_csr_e[cnt] << endl;
                sca_csr_ev[cnts++] = vout_scalar[i*nc+j];
            }
        }
    }

	FILE* fout = fopen("easy.ans", "a");
	for (int i=0; i < cnt; i++){	
		fprintf(fout, "%5.2f ", sca_csr_ev[i]);	
		// fprintf(fout, "%d ", sca_csr_e[i]);	

	}
	fprintf(fout, "\n");
	
	FTYPE *vec_csr_ev = new FTYPE [cnt];
	FTYPE *tmp_vec = vec_csr_ev;
	for (int i=0; i < nr; i++){	
		memcpy(tmp_vec, M3._rowPtr[i], M3._nnz_num[i]*4);
		tmp_vec+=M3._nnz_num[i];
		// delete [] M3._rowPtr[i];
		for(int j=0; j<M3._nnz_num[i]; j++){
			fprintf(fout, "%5.2f ", M3._rowPtr[i][j]);	
			// fprintf(fout, "%d ", M3._coordPtr[i][j]);	
		}
	}
	fprintf(fout, "\n");
	for(int i=0; i<cnt; i++){
		if(abs(sca_csr_ev[i]-vec_csr_ev[i])>0.01){
			cout << i << " " << sca_csr_ev[i] << " " << vec_csr_ev[i] << endl;
			break;	
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
 