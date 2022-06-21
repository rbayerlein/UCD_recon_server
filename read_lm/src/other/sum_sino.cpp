// Sums sinograms from each of 8 EXPLORER nodes into one sinogram 


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
// #include <time.h>
// #include <omp.h>

#define PI 3.141592653589793

using namespace std;



int main(int argc, char** argv) {

	
	if (argc < 3) {
		cout << "not enough input arguments, quit\n"; 
		return 1; 
	}

	cout << argc << "\n"; 

	string infile_fullpath;
	string outfile_fullpath; 
	string outfolder;
	string str_temp;

	
	int din = 0; 


	
	// open first sinogram
	infile_fullpath = argv[2]; 

	ifstream sino_in1; 
	sino_in1.open(infile_fullpath.c_str(), ios::in | ios::binary); 

	ifstream sino_in2;
	ifstream sino_in3;
	ifstream sino_in4;
	ifstream sino_in5;
	ifstream sino_in6;
	ifstream sino_in7;
	ifstream sino_in8; 


	sino_in1.seekg(0, sino_in1.end);
	long file_size = sino_in1.tellg(); 
	long sino_size = file_size / 4; 
	sino_in1.seekg(0, sino_in1.beg); 


	if (argc > 3) {
		infile_fullpath = ""; 
		infile_fullpath = argv[3]; 
		sino_in2.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in2.seekg(0, sino_in2.end);
		if (sino_in2.tellg() != file_size) {
			cout << "Invalid sinogram 2\n"; 
			return 1; 
		}
		sino_in2.seekg(0, sino_in2.beg); 
		
	}
	if (argc > 4) {
		infile_fullpath = ""; 
		infile_fullpath = argv[4]; 
		sino_in3.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in3.seekg(0, sino_in3.end);
		if (sino_in3.tellg() != file_size) {
			cout << "Invalid sinogram 3\n"; 
			return 1; 
		}
		sino_in3.seekg(0, sino_in3.beg);  
	}
	if (argc > 5) {
		infile_fullpath = ""; 
		infile_fullpath = argv[5]; 
		sino_in4.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in4.seekg(0, sino_in4.end);
		if (sino_in4.tellg() != file_size) {
			cout << "Invalid sinogram 4\n"; 
			return 1; 
		}
		sino_in4.seekg(0, sino_in4.beg);  
	}
	if (argc > 6) {
		infile_fullpath = ""; 
		infile_fullpath = argv[6]; 
		sino_in5.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in5.seekg(0, sino_in5.end);
		if (sino_in5.tellg() != file_size) {
			cout << "Invalid sinogram 5\n"; 
			return 1; 
		}
		sino_in5.seekg(0, sino_in5.beg);  
	}
	if (argc > 7) {
		infile_fullpath = ""; 
		infile_fullpath = argv[7]; 
		sino_in6.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in6.seekg(0, sino_in6.end);
		if (sino_in6.tellg() != file_size) {
			cout << "Invalid sinogram 6\n"; 
			return 1; 
		}
		sino_in6.seekg(0, sino_in6.beg);  
	}
	if (argc > 8) {
		infile_fullpath = ""; 
		infile_fullpath = argv[8]; 
		sino_in7.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in7.seekg(0, sino_in7.end);
		if (sino_in7.tellg() != file_size) {
			cout << "Invalid sinogram 7\n"; 
			return 1; 
		}
		sino_in7.seekg(0, sino_in7.beg);  
	}
	if (argc > 9) {
		infile_fullpath = ""; 
		infile_fullpath = argv[9]; 
		sino_in8.open(infile_fullpath.c_str(), ios::in | ios::binary);
		sino_in8.seekg(0, sino_in8.end);
		if (sino_in8.tellg() != file_size) {
			cout << "Invalid sinogram 8\n"; 
			return 1; 
		}
		sino_in8.seekg(0, sino_in8.beg);  
	}


	// open output sinogram
	outfile_fullpath = argv[1]; 
	//int str_len = outfile_fullpath.length(); 
	//outfile_fullpath.erase(str_len-6, 6);
	//outfile_fullpath.append(".raw"); 

	cout << "sino out name = " << outfile_fullpath << "\n"; 

	vector<int> sino_sum(sino_size); 

	for (long i = 0; i < sino_size; i++) {
		sino_sum[i] = 0; 
	}





	// Main run program
	for (long k = 0; k < sino_size; k++) {
		sino_in1.read(reinterpret_cast<char *>(&din), sizeof(int));
		sino_sum[k] = sino_sum[k] + din; 
		din = 0; 

		if (argc > 2) {
			sino_in2.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 3) {
			sino_in3.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 4) {
			sino_in4.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 5) {
			sino_in5.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 6) {
			sino_in6.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 7) {
			sino_in7.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}
		if (argc > 8) {
			sino_in8.read(reinterpret_cast<char *>(&din), sizeof(int));
			sino_sum[k] = sino_sum[k] + din; 
			din = 0; 
		}




	}


	// write summed sinogram
	ofstream sino_out; 
	sino_out.open(outfile_fullpath.c_str(), ios::out | ios::binary); 

	for (long j = 0; j < sino_size; j++) {
		sino_out.write(reinterpret_cast<const char *>(&sino_sum[j]), sizeof(int));
	}
	sino_out.close(); 
	
	
	sino_in1.close(); 
	
	if (argc > 2) {
		sino_in2.close();  
	}
	if (argc > 3) {
		sino_in3.close();  
	}
	if (argc > 4) {
		sino_in4.close();
	}
	if (argc > 5) {
		sino_in5.close(); 
	}
	if (argc > 6) {
		sino_in6.close();
	}
	if (argc > 7) {
		sino_in7.close(); 
	}
	if (argc > 8) {
		sino_in8.close();  
	}







return 1; 

}
