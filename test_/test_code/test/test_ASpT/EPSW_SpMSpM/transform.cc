#include<iostream>
#include<cstring>
#include<vector>
#define FTYPE float
using namespace std;
int nr, nc, ne;
int *csr_v;
int *csr_e;
FTYPE *csr_ev;

int *csr_vt;
int *csr_et;
FTYPE *csr_evt;

char *dir_r;
char *app;

void ready(int argc, char **argv)
{
	FILE *fp_r;
	int i;

	app = argv[0];
	dir_r = argv[1];
	fp_r = fopen(dir_r, "r");
	fscanf(fp_r, "%d, %d, %d", &nr, &nc, &ne);

	csr_v = new int[nr + 1]{0};	  
	csr_e = new int[ne]{0};	   
	csr_ev = new FTYPE[ne]{0.0}; 

	csr_vt = new int[nc + 1]{0};	  
	csr_et = new int[ne]{0};	   
	csr_evt = new FTYPE[ne]{0.0}; 

	for (i = 0; i < nr+1; i++){
		fscanf(fp_r, "%d ", &csr_v[i]);
	}

	for (i = 0; i < ne; i++){
		fscanf(fp_r, "%d", &csr_e[i]);
	}
	fclose(fp_r);
	// fprintf(stdout, "%d %d %d \n", ne, nr, nc);
}
void gen(){
	for(int j=0; j<ne; j++){
		csr_vt[csr_e[j]+1]++;
	}
	
	for(int i=0; i<nc+1; i++){	
		// cout << csr_vt[i] << " ";
		csr_vt[i+1] += csr_vt[i];
	}
	// cout << endl;
	vector<vector<int>> csr_etv(nc);
	for(int i=0; i<nc; i++){
		vector<int> tmp;
		csr_etv.push_back(tmp);
	}
	vector<vector<FTYPE>> csr_evtv(nc);
	for(int i=0; i<nc; i++){
		vector<FTYPE> tmp;
		csr_evtv.push_back(tmp);
	}


	for(int i=0; i<nr; i++){
		for(int j=csr_v[i]; j<csr_v[i+1]; j++){
			csr_etv[csr_e[j]].push_back(i);
			csr_evtv[csr_e[j]].push_back(csr_ev[j]);
		}
	}

	int cnt=0;
	for(int i=0; i<nc+1; i++){
		for(int j=0; j<csr_etv[i].size(); j++){
			csr_et[cnt] = csr_etv[i][j];
			csr_evt[cnt++] = csr_evtv[i][j];
		}
	}
    // string dir(dir_r);
    int len = strlen(dir_r);
    char* dir = new char[len+6];
    strcpy(dir, dir_r);
    strcpy(dir+len, ".trans");
    FILE* fp_w = fopen(dir, "a");
    fprintf(fp_w, "%d,%d,%d \n", nc, nr, ne);
    for(int i=0; i<nc+1; i++){
        fprintf(fp_w, "%d ", csr_vt[i]);
    }
    fprintf(fp_w, "\n");
    for(int i=0; i<ne; i++){
        fprintf(fp_w, "%d ", csr_et[i]);
    }
}
int main(int argc, char **argv){
    ready(argc, argv);
    gen();
}
