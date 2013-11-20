#if !defined(FILEIO_UTILS_HPP)
#define FILEIO_UTILS_HPP

#include <sstream>
#include <string.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

template <class Tp> 
void parser(vector<Tp>& output,string str, string delims=" \t") {
	char * str_c = new char[str.size()+1];
	str_c[str.size()] = '\0';
	memcpy(str_c,str.c_str(),str.size());

	char * pch = strtok(str_c,delims.c_str());
	while(pch !=NULL) {
		istringstream iss(pch);
		Tp tok;
		if(!(iss>>tok)) {
			cout << "# Error parsing string!" << endl;
			exit(-1);
		}
		output.push_back(tok);
		pch = strtok(NULL,delims.c_str());
	}
	delete[] str_c;
}

#endif
