#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<cstring>
using namespace std;

int main(int argc, char **argv){
    fstream fin;
    fin.open(argv[1]);
	if(!fin)
	{
		cout<<"can not open" << argv[1] <<endl;
		exit(1);
	}
    string line;
    while(getline(fin, line)){
        // string s(argv[2]);
        char *sub = strstr(const_cast<char*>(line.c_str()), const_cast<char*>(argv[2]));
        if(sub)
        cout << sub << endl;
    }


    return 0;
}