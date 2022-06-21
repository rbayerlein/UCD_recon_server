// Inserts coincidence data from one dataset (e.g. point source) into another dataset
// Eric Berg; July 2019

// inputs are:
// 1) full path of the new listmode file
// 2) full path of the primary listmode file
// 3) full path of the secondary listmode file to be inserted



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <chrono>
// #include <omp.h>

#define PI 3.141592653589793

using namespace std;



int main(int argc, char** argv) {
	//cout << "You have entered " << argc << " arguments:" << "\n";

 	if (argc < 2) {
 		cout << "Not enough input arguments, try again!";
 		return 0;
 	}
 	if (argc > 10) {
 		cout << "Too many input arguments, try again!";
 		return 0;
 	}

	// variables that are set by user in GUI
	string outname;
	string outname_p; 
	string outname_s; 
	string outname_m; 

    string inname_p1;
    string inname_p2;
    string inname_p3;
    string inname_p4;
    string inname_p5;
    string inname_p6;
    string inname_p7;
    string inname_p8;

    string inname_s1;
    string inname_s2;
    string inname_s3;
    string inname_s4;
    string inname_s5;
    string inname_s6;
    string inname_s7;
    string inname_s8;

    string inname_m1;
    string inname_m2;
    string inname_m3;
    string inname_m4;
    string inname_m5;
    string inname_m6;
    string inname_m7;
    string inname_m8;


	outname = argv[1];
	inname_p1 = argv[2];
	inname_p2 = argv[3];
	inname_p3 = argv[4];
	inname_p4 = argv[5];
	inname_p5 = argv[6];
	inname_p6 = argv[7];
	inname_p7 = argv[8];
	inname_p8 = argv[9];

	// int num_outfiles = argv[10];  
	int num_outfiles = 300; 

	inname_s1 = inname_p1; 
	inname_s1.erase(inname_p1.length()-2, 2); 
	inname_s1.append("add_fac"); 

	inname_s2 = inname_p2; 
	inname_s2.erase(inname_p2.length()-2, 2); 
	inname_s2.append("add_fac"); 

	inname_s3 = inname_p3; 
	inname_s3.erase(inname_p3.length()-2, 2); 
	inname_s3.append("add_fac"); 

	inname_s4 = inname_p4; 
	inname_s4.erase(inname_p4.length()-2, 2); 
	inname_s4.append("add_fac"); 

	inname_s5 = inname_p5; 
	inname_s5.erase(inname_p5.length()-2, 2); 
	inname_s5.append("add_fac"); 

	inname_s6 = inname_p6; 
	inname_s6.erase(inname_p6.length()-2, 2); 
	inname_s6.append("add_fac"); 

	inname_s7 = inname_p7; 
	inname_s7.erase(inname_p7.length()-2, 2); 
	inname_s7.append("add_fac"); 

	inname_s8 = inname_p8; 
	inname_s8.erase(inname_p8.length()-2, 2); 
	inname_s8.append("add_fac"); 




	inname_m1 = inname_p1; 
	inname_m1.erase(inname_p1.length()-2, 2); 
	inname_m1.append("mul_fac"); 

	inname_m2 = inname_p2; 
	inname_m2.erase(inname_p2.length()-2, 2); 
	inname_m2.append("mul_fac"); 

	inname_m3 = inname_p3; 
	inname_m3.erase(inname_p3.length()-2, 2); 
	inname_m3.append("mul_fac"); 

	inname_m4 = inname_p4; 
	inname_m4.erase(inname_p4.length()-2, 2); 
	inname_m4.append("mul_fac"); 

	inname_m5 = inname_p5; 
	inname_m5.erase(inname_p5.length()-2, 2); 
	inname_m5.append("mul_fac"); 

	inname_m6 = inname_p6; 
	inname_m6.erase(inname_p6.length()-2, 2); 
	inname_m6.append("mul_fac"); 

	inname_m7 = inname_p7; 
	inname_m7.erase(inname_p7.length()-2, 2); 
	inname_m7.append("mul_fac"); 

	inname_m8 = inname_p8; 
	inname_m8.erase(inname_p8.length()-2, 2); 
	inname_m8.append("mul_fac"); 



	//outname_p = outname; 

	//outname_s = outname; 
	//outname_s.erase(outname.length()-2, 2); 
	//outname_s.append("add_fac"); 

	//outname_m = outname; 
	//outname_m.erase(outname.length()-2, 2); 
	//outname_m.append("mul_fac"); 


	


	////////////////////////////////////////////////////////////////////////////////
	//  ********************* Scanner parameters and global variables *********  //

	

	//////////////////////////////////////////////////////////////////////////////

	// ********************  Read Listmode  *********************************** //


	// open listmode files
	ifstream infile_p1;
	infile_p1.open(inname_p1.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p2;
	infile_p2.open(inname_p2.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p3;
	infile_p3.open(inname_p3.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p4;
	infile_p4.open(inname_p4.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p5;
	infile_p5.open(inname_p5.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p6;
	infile_p6.open(inname_p6.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p7;
	infile_p7.open(inname_p7.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_p8;
	infile_p8.open(inname_p8.c_str(), ios::in | ios::binary); //open primary list mode file


	if (!infile_p1 || !infile_p2 || !infile_p3 || !infile_p4 || !infile_p5 || !infile_p6 || !infile_p7 || !infile_p8 ) {
		infile_p1.close();
		cout << inname_p1 << "\n";
		cout << "Cannot open primary listmode file\nCheck folder names and locations\n";

		return 0;
	}



	// open listmode files
	ifstream infile_s1;
	infile_s1.open(inname_s1.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s2;
	infile_s2.open(inname_s2.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s3;
	infile_s3.open(inname_s3.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s4;
	infile_s4.open(inname_s4.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s5;
	infile_s5.open(inname_s5.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s6;
	infile_s6.open(inname_s6.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s7;
	infile_s7.open(inname_s7.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_s8;
	infile_s8.open(inname_s8.c_str(), ios::in | ios::binary); //open primary list mode file







	// open listmode files
	ifstream infile_m1;
	infile_m1.open(inname_m1.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m2;
	infile_m2.open(inname_m2.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m3;
	infile_m3.open(inname_m3.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m4;
	infile_m4.open(inname_m4.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m5;
	infile_m5.open(inname_m5.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m6;
	infile_m6.open(inname_m6.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m7;
	infile_m7.open(inname_m7.c_str(), ios::in | ios::binary); //open primary list mode file

	// open listmode files
	ifstream infile_m8;
	infile_m8.open(inname_m8.c_str(), ios::in | ios::binary); //open primary list mode file



	// open listmode files
	//ofstream outfile;
	ofstream outfile[300]; 
	ofstream outfile_s[300]; 
	ofstream outfile_m[300];


	 

	string outname_sub; 
	string outname_sub_s; 
	string outname_sub_m; 



	for (int fi = 0; fi < num_outfiles; fi++)  {
		stringstream subf;  
		stringstream subf_m; 
		stringstream subf_s; 

		subf << (fi + 1) << "_prompts.lm"; 
		subf_m << (fi + 1) << "_prompts.mul_fac"; 
		subf_s  << (fi + 1) << "_prompts.add_fac"; 

		outname_sub = ""; 
		outname_sub = outname; 
		outname_sub.append(subf.str()); 

		outname_sub_s = ""; 
		outname_sub_s = outname; 
		outname_sub_s.append(subf_s.str()); 

		outname_sub_m = ""; 
		outname_sub_m = outname; 
		outname_sub_m.append(subf_m.str()); 

		subf << ""; 
		subf.clear(); 
		subf_m << ""; 
		subf_m.clear(); 
		subf_s << ""; 
		subf_s.clear(); 

		

		outfile[fi].open(outname_sub.c_str(), ios::out  | ios::binary);
		outfile_m[fi].open(outname_sub_m.c_str(), ios::out  | ios::binary);
		outfile_s[fi].open(outname_sub_s.c_str(), ios::out  | ios::binary); 
	}


	//outfile.open(outname.c_str(), ios::out | ios::binary); //open primary list mode file


	//ofstream outfile_s;
	//outfile_s.open(outname_s.c_str(), ios::out | ios::binary); //open primary list mode file



	//ofstream outfile_m;
	//outfile_m.open(outname_m.c_str(), ios::out | ios::binary); //open primary list mode file

	

	long long max_events_file = 0; 

	infile_p1.seekg(0, infile_p1.end);
	long long file_size_p1 = infile_p1.tellg(); //get size (events) of list mode file
	long long num_events_p1 = file_size_p1 / 10;
	max_events_file = num_events_p1; 
	infile_p1.seekg(0, infile_p1.beg);

	infile_p2.seekg(0, infile_p2.end);
	long long file_size_p2 = infile_p2.tellg(); //get size (events) of list mode file
	long long num_events_p2 = file_size_p2 / 10;
	if (num_events_p2 > max_events_file) { 
		max_events_file = num_events_p2; 
	}
	infile_p2.seekg(0, infile_p2.beg);

	infile_p3.seekg(0, infile_p3.end);
	long long file_size_p3 = infile_p3.tellg(); //get size (events) of list mode file
	long long num_events_p3 = file_size_p3 / 10;
	if (num_events_p3 > max_events_file) { 
		max_events_file = num_events_p3; 
	}
	infile_p3.seekg(0, infile_p3.beg);

	infile_p4.seekg(0, infile_p4.end);
	long long file_size_p4 = infile_p4.tellg(); //get size (events) of list mode file
	long long num_events_p4 = file_size_p4 / 10;
	if (num_events_p4 > max_events_file) { 
		max_events_file = num_events_p4; 
	}
	infile_p4.seekg(0, infile_p4.beg);

	infile_p5.seekg(0, infile_p5.end);
	long long file_size_p5 = infile_p5.tellg(); //get size (events) of list mode file
	long long num_events_p5 = file_size_p5 / 10;
	if (num_events_p5 > max_events_file) { 
		max_events_file = num_events_p5; 
	}
	infile_p5.seekg(0, infile_p5.beg);

	infile_p6.seekg(0, infile_p6.end);
	long long file_size_p6 = infile_p6.tellg(); //get size (events) of list mode file
	long long num_events_p6 = file_size_p6 / 10;
	if (num_events_p6 > max_events_file) { 
		max_events_file = num_events_p6; 
	}
	infile_p6.seekg(0, infile_p6.beg);

	infile_p7.seekg(0, infile_p7.end);
	long long file_size_p7 = infile_p7.tellg(); //get size (events) of list mode file
	long long num_events_p7 = file_size_p7 / 10;
	if (num_events_p7 > max_events_file) { 
		max_events_file = num_events_p7; 
	}
	infile_p7.seekg(0, infile_p7.beg);

	infile_p8.seekg(0, infile_p8.end);
	long long file_size_p8 = infile_p8.tellg(); //get size (events) of list mode file
	long long num_events_p8 = file_size_p8 / 10;
	if (num_events_p8 > max_events_file) { 
		max_events_file = num_events_p8; 
	}
	infile_p8.seekg(0, infile_p8.beg);




	long long num_events_all = num_events_p1 + num_events_p2 + num_events_p3 + num_events_p4 + num_events_p6 + num_events_p7 + num_events_p8; 

	double p1_prob = (double)num_events_p1 / (double)max_events_file; 
	double p2_prob = (double)num_events_p2 / (double)max_events_file; 
	double p3_prob = (double)num_events_p3 / (double)max_events_file; 
	double p4_prob = (double)num_events_p4 / (double)max_events_file; 
	double p5_prob = (double)num_events_p5 / (double)max_events_file; 
	double p6_prob = (double)num_events_p6 / (double)max_events_file; 
	double p7_prob = (double)num_events_p7 / (double)max_events_file; 
	double p8_prob = (double)num_events_p8 / (double)max_events_file; 


	



	//////////////////////////////////////////////////////////////////////////////

	// **************		Main Run Program		*********************//
// *****************************************************************************************************************************************************************************


	//cout << "\n\n*** Processing Dataset # ... *** \n\n";

	// for reading lm file
	short din = 0; 
	float din_f = 0.0; 
int f_towrite = 0; 
    double r_num, r_file = 0.0; 
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	default_random_engine generator (seed);
	//default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0,1.0);



	//unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count(); 
	//default_random_engine generator_f (seed2); 


	//uniform_real_distribution2<double> distribution(0.0,1.0);


	bool notdone = true; 
	bool notdone1 = true; 
	bool notdone2 = true; 
	bool notdone3 = true; 
	bool notdone4 = true; 
	bool notdone5 = true; 
	bool notdone6 = true; 
	bool notdone7 = true; 
	bool notdone8 = true; 

	while (notdone) {
	//for (long long nn = 0; nn < num_events_p; nn++) {

		r_num = distribution(generator);

		if (r_num < p1_prob && notdone1) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p1.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s1.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m1.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 
			 
			if (infile_p1.tellg() > (file_size_p1-10)) {
				notdone1 = false; 
				infile_p1.close(); 
				infile_s1.close(); 
				infile_m1.close(); 
			}
		}

		if (r_num < p2_prob && notdone2) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p2.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s2.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m2.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p2.tellg() > (file_size_p2-10)) {
				notdone2 = false; 
				infile_p2.close(); 
				infile_s2.close(); 
				infile_m2.close(); 
			}
		}

		if (r_num < p3_prob && notdone3) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p3.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s3.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
                        din_f = din_f / (float)num_outfiles; 
    			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m3.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p3.tellg() > (file_size_p3-10)) {
				notdone3 = false;
				infile_p3.close(); 
				infile_s3.close(); 
				infile_m3.close();  
			}
		}

		if (r_num < p4_prob && notdone4) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p4.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s4.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m4.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p4.tellg() > (file_size_p4-10)) {
				notdone4 = false; 
				infile_p4.close();
				infile_s4.close(); 
				infile_m4.close();  
			}
		}

		if (r_num < p5_prob && notdone5) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p5.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s5.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
                        din_f = din_f / (float)num_outfiles; 
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m5.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p5.tellg() > (file_size_p5-10)) {
				notdone5 = false; 
				infile_p5.close();
				infile_s5.close(); 
				infile_m5.close();  
			}
		}

		if (r_num < p6_prob && notdone6) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p6.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}
			infile_s6.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m6.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p6.tellg() > (file_size_p6-10)) {
				notdone6 = false; 
				infile_p6.close(); 
				infile_s6.close(); 
				infile_m6.close(); 
			}
		}

		if (r_num < p7_prob && notdone7) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p7.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}
			infile_s7.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m7.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p7.tellg() > (file_size_p7-10)) {
				notdone7 = false; 
				infile_p7.close(); 
				infile_s7.close(); 
				infile_m7.close(); 
			}
		}

		if (r_num < p8_prob && notdone8) {
			r_file = (double)num_outfiles * distribution(generator);
			f_towrite  = floor(r_file); 
			for (int k = 0; k < 5; k++) {
				infile_p8.read(reinterpret_cast<char *>(&din), sizeof(short));
				outfile[f_towrite].write(reinterpret_cast<const char *>(&din), sizeof(short)); 
			}

			infile_s8.read(reinterpret_cast<char *>(&din_f), sizeof(float));
                        din_f = din_f / (float)num_outfiles;  
			outfile_s[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			infile_m8.read(reinterpret_cast<char *>(&din_f), sizeof(float)); 
			outfile_m[f_towrite].write(reinterpret_cast<const char *>(&din_f), sizeof(float)); 

			if (infile_p8.tellg() > (file_size_p8-10)) {
				notdone8 = false; 
				infile_p8.close(); 
				infile_s8.close(); 
				infile_m8.close(); 
			}
		}


		if (!notdone1 && !notdone2 && !notdone3 && !notdone4 && !notdone5 && !notdone6 && !notdone7 && !notdone8) {
			notdone = false; 
		}




	}



	for (int fi = 0; fi < num_outfiles; fi++)  {
		outfile[fi].close();
		outfile_m[fi].close();
		outfile_s[fi].close(); 
	}

	//outfile.close();
	//outfile_s.close(); 
	//outfile_m.close(); 


	cout << "Listmode combine finished!\n";


	return 0;



}
