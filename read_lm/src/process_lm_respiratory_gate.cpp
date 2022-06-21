

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
// #include <time.h>
// #include <omp.h>

#define PI 3.141592653589793

using namespace std;



int main(int argc, char** argv) {


	/////////////////////////////////////////////////////////////////////////////////////

	// ******************* File info ********************************************* ///
	
	cout << "\n"; 
	
	// variables that are set by user in GUI

	// file names, directories
	//string data_base_dir = "c:/Documents/Primate scanner/monkey/";
	string infile_fullpath;
	string outfolder;
	string str_temp;

	// dynamic framing
	string static_dynamic;
	bool dyn = false;
	vector<int> frames;
	vector<float> frame_t;
	bool ft1 = true;
	int vec_len = 0;

	bool write_lmfile = false;
	bool write_sino = false;
	bool write_sino_block = true; 
	bool write_histo_img = true; 
	
	bool rand_sub = false; 

	// Read config file
	//string scan_details_fullpath = "/run/media/meduser/data/read_lm/Reconstruction_Parameters_1";
	string scan_details_fullpath;
	scan_details_fullpath = argv[1];  

	//string scan_details_fullpath = "c:/Documents/Software/miniEXPLORER/process_lm_sino/Reconstruction_Parameters_1.txt";

	ifstream detes;
	detes.open(scan_details_fullpath.c_str());
	if (!detes) {
		cout << "could not open read listmode config file\n";
		return 1;
	}

	int subint = 0;
	float subtime = 0.0;

	// path to output folder
	getline(detes, str_temp);
	// cout << str_temp << "\n";
	outfolder = str_temp;
	//cout << outfolder << "\n"; 
	//outfolder.append("/");
	str_temp = "";

	//string scan_details_fullpath_out = outfolder; 
	//scan_details_fullpath_out.append("Reconstruction_Parameters_1"); 
	//ofstream detesout; 
	//detesout.open(scan_details_fullpath_out.c_str()); 
	//if (!detesout) {
	//	cout << "Could not create output scan log\n";
    //}
	//detesout << outfolder << "\n"; 


	// path to .lm data
	getline(detes, str_temp);
	//detesout << str_temp << "\n";
	infile_fullpath = str_temp;
	
	str_temp = "";

	// get raw data number (1 to 8)
	string str_find = ".raw"; 
	size_t ext_ind = infile_fullpath.find(str_find);
	char raw_num_char = infile_fullpath[ext_ind - 1]; 
	string raw_num(1, raw_num_char); 
	cout << "File number = " << raw_num << "\n"; 

	// Static or dynamic (1,2)
	getline(detes, str_temp);
	//detesout << str_temp << "\n";
	stringstream convert2(str_temp);
	int subint2;
	convert2 >> subint2;
	bool run_recon = false; 
	if (subint2 == 1) {
		run_recon = true; 
	}
	
	
	convert2.str("");
	convert2.clear();
	

	str_temp = "";

	// dynamic framing (1, 100, 2, 200, etc)
	getline(detes, str_temp);
	//detesout << str_temp << "\n";
	stringstream convert3(str_temp);
	int subint3;
	convert3 >> subint3;
	static_dynamic = str_temp; 

	if (subint3 == 0) {
		dyn = false;
	}

	else {
		dyn = true;
	}
	stringstream sd(static_dynamic);
	while (sd.good()) {
		string substr;

		getline(sd, substr, ',');
		//cout << substr << "\n";
		stringstream convert3(substr);

		if (ft1) {
			convert3 >> subint;
			convert3.str("");
			convert3.clear();
			frames.push_back((int)subint);
			vec_len++;
		}
		else {
			convert3 >> subtime;
			convert3.str("");
			convert3.clear();
			frame_t.push_back((double)subtime);
		}
		ft1 = !ft1;
	}
	str_temp = "";

	// write processed listmode data (0 or 1)
	getline(detes, str_temp);
	//detesout << str_temp << "\n"; 
	stringstream convert4(str_temp);
	convert4 >> subint;
	convert4.str("");
	convert4.clear();
	if (subint == 1) {
		write_lmfile = true;
	}
	else {
		write_lmfile = false;
	}
	subint = 0;
	str_temp = ""; 

	// write 3D sinograms (0 or 1)
	getline(detes, str_temp);
	//detesout << str_temp << "\n";
	stringstream convert5(str_temp);
	convert5 >> subint;
	convert5.str("");
	convert5.clear();
	if (subint == 1) {
		write_sino = true;
		cout << "write sino TRUE\n"; 
	}
	else {
		write_sino = false;
		cout << "write sino FALSE\n"; 
	}
	subint = 0;
	str_temp = ""; 
	
	
	// write block sinograms (0 or 1)
	getline(detes, str_temp);
	//detesout << str_temp << "\n";
	stringstream convert6(str_temp);
	convert6 >> subint;
	convert6.str("");
	convert6.clear();
	if (subint == 1) {
		write_sino_block = true;
		cout << "write block sino TRUE\n"; 
	}
	else {
		write_sino_block = false;
		cout << "write block sino FALSE\n"; 
	}
	subint = 0;
	str_temp = "";

	// write TOF histo image
	getline(detes, str_temp);
	stringstream convert7(str_temp);
	convert7 >> subint;
	convert7.str("");
	convert7.clear();
	if (subint == 1) {
		write_histo_img = true;
		cout << "write TOF HISTO IMG = TRUE\n"; 
	}
	else {
		write_histo_img = false;
		cout << "write TOF HISTO IMG = FALSE\n"; 
	}
	subint = 0;
	str_temp = "";


	//detesout.close(); 
	detes.close(); 
	

	// dynamic framing 
	int tot_frames = 0;
	for (int i = 0; i < vec_len; i++) {
		tot_frames = tot_frames + frames[i];
	}
	if (dyn) {
		//cout << "Dynamic imaging\n"; 
		cout << "total frames = " << tot_frames << "\n";
	}
	else {
		cout << "Static imaging\n"; 
	}
	vector<int> frame_length(tot_frames);
	vector<double> frame_length_s(tot_frames);
	int frame_count = 0; 
	int cur_pos = 0;
	for (int i = 0; i < tot_frames; i++) {
		if (i == (frames[cur_pos] + frame_count)) {
			frame_count = frame_count + frames[cur_pos];
			cur_pos++;

		}
		frame_length[i] = (int)(1000*frame_t[cur_pos]);
		frame_length_s[i] = (frame_t[cur_pos]);
		
		if (dyn) {
			cout << "Frame " << i << " = " << frame_length_s[i] << " seconds\n";

		}
	}


	
	////////////////////////////////////////////////////////////////////////////////
	//  ********************* Scanner parameters and global variables *********  //

	// scanner parameters
	int num_bins_sino = 549*420; 
	int num_bins_sino2 = num_bins_sino / 6; 
	int num_bins_sino_block = 91*60; 
	int num_slice = 1343; 
	int num_slice_block = 223; 
	int num_crystals = 840;
	int num_crystals_ax_wgap = 679; 
	int num_crystals_all = 564480;

	double tof_res = 450.0;
	
	//int num_crystals_all = 94080; 
	
	// for reading lm file	
	
	unsigned int din1, din2, dout1, dout2 = 0;
	unsigned int a32 = 2147483648;
  
	int p_d, coinc = 0; 
	bool coinc_tag = true; 
	bool blockcount_tag = false; 
	bool time_tag = false; 
	bool time_tag0 = false; 
	int bk1, bk2, bk3 = 0; 

	
	
	bool corrupted = false; 
	signed short dtemp = 0; 
	short transcA, transcB, crysaxA, crysaxB, dt = 0; 
	short dout[5];
	unsigned short eout [2];  

	
	// time tag variables
	int year1, month1, day1, hour1, minute1, second1, milli1 = 0;
	int year0, month0, day0, hour0, minute0, second0, milli0 = 0;
	int year00, month00, day00, hour00, minute00, second00, milli00 = 0; 
	int ss0, ss1, ssdiff = 0; 

	
	
	// dynamic 
	int frame_start = 0.0;
	int frame_end = 0.0;
	int frame_num = 0; 
	int time_s, time_s0 = 0; 
	
	frame_end = frame_start + frame_length[0];
	
	
	// to run 
	bool run = true; 
	bool write_lm = false; 


	//////////////////////////////////////////////////////////////////////////////
	
	// ********************  Read Listmode  *********************************** // 
	
	
	// open listmode file
	ifstream infile1;
	infile1.open(infile_fullpath.c_str(), ios::in | ios::binary); //open list mode file
	

	if (!infile1) {
		infile1.close();
		cout << infile_fullpath << "\n"; 
		cout << "Cannot open listmode file\nCheck folder names and locations\n";
	}
	
	
	long long file_size = 0; 
		
	infile1.seekg(0, infile1.end);
	long long file_size1 = infile1.tellg(); //get size (events) of list mode file
	long long num_events = file_size1 / 8;
	cout << num_events << " total events\n"; 
	infile1.seekg(0, infile1.beg);
	


	string fname_out;
	string outfile_fullpath = "";

	stringstream ss;
	ss << "lm_reorder_f" << frame_num;
	fname_out = ss.str();
	outfile_fullpath = outfolder;
	outfile_fullpath.append(fname_out);
	
	outfile_fullpath.append("_prompts."); 
	outfile_fullpath.append(raw_num); 
	outfile_fullpath.append(".raw");
	
	ss << "";
	ss.clear();
	str_temp = "";
	fname_out = "";

	ofstream outfile;	//prompts output
	outfile.open(outfile_fullpath.c_str(), ios::out | ios::binary); //create binary file containing new crystal + time data


	cout << "\nFinding start time...\n"; 
	
	bool run1 = true;  

	while (run1) {
		
		
		
		infile1.read(reinterpret_cast<char *>(&din1), sizeof(uint));
		infile1.read(reinterpret_cast<char *>(&din2), sizeof(uint));
		
		
		
		if (!(din1<=a32) && !(din2>=a32)) {
			// shift by one byte at a time
			infile1.read(reinterpret_cast<char *>(&dtemp), 1);
		}
		else {
			din2 = din2 - a32; 
				
			coinc = din1 >> 4;
			coinc = coinc % 2; 
			if (coinc == 1) {
				coinc_tag = true; 
			}
			else {
				coinc_tag = false;
			}
			// check for blk count tag
			bk1 = din1 >> 4; 
			bk2 = din1 >> 5; 
			bk3 = din1 >> 6; 
			
			bk1 = bk1 % 2; 
			bk2 = bk2 % 2; 
			bk3 = bk3 % 2; 
			
			// time tag
			if (bk1 == 0 && bk2 == 0 && bk3 == 0) {
					
					 
				day1 = din1 >> 19; 
				din1 = din1 - (day1<<19);
				if (day1 > 512) {
					day1 = day1 - 512; 
				}
					 
				year1 = din1 >> 7; 
				din1 = din1 - (year1<<7);
				month1 = din1; 
					
				milli1 = din2 >> 19; 
				din2 = din2 - (milli1 << 19);
				second1 = din2 >> 13; 
				din2 = din2 - (second1 << 13); 
				minute1 = din2 >> 7; 
				din2 = din2 - (minute1<<7); 
				hour1 = din2;  
				
				if (1) {
					year00 = year1; 
					month00 = month1; 
					day00 = day1; 
					hour00 = hour1; 
					minute00 = minute1; 
					second00 = second1; 
					milli00 = milli1;
				}
				else if (year1 != year00 || month1 != month00 || day1 != day00 || hour1 != hour00 || minute1 != minute00 || second1 != second00 || milli1 != milli00) {
					cout << "listmode files not aligned in time!\n"; 
				}
				
				run1 = false; 
			}
		}
	
	}
	
	cout << "Scan start: " << year1 << "/" << month1 << "/" << day1 << " " << hour1 << ":" << minute1 << ":" << second1 << "." << milli1 << "\n";
	
	ss0 = milli00; 
	infile1.seekg(-8, infile1.cur);
	
	
	cout << "\n\n*** Processing Frame 0... *** \n\n"; 
	
 	bool next_frame = false; 
 	time_s = time_s0; 
 
	while (run) {
		
		infile1.read(reinterpret_cast<char *>(&din1), sizeof(unsigned int));
		infile1.read(reinterpret_cast<char *>(&din2), sizeof(unsigned int));
		
		dout1 = din1; 
		dout2 = din2; 
		if (infile1.tellg() > (file_size1-8*10)) {
			run = false; 
		}
		
				
		if (!(din1<=a32) && !(din2>=a32)) {
			// shift by one byte at a time
			infile1.read(reinterpret_cast<char *>(&dtemp), 1);				
		}

		else {
			din2 = din2 - a32; 
				
			coinc = din1 >> 4;
			coinc = coinc % 2; 
			if (coinc == 1) {
				coinc_tag = true; 
			}
			else {
				coinc_tag = false;
			}
			// check for blk count tag
			bk1 = din1 >> 4; 
			bk2 = din1 >> 5; 
			bk3 = din1 >> 6; 
			
			bk1 = bk1 % 2; 
			bk2 = bk2 % 2; 
			bk3 = bk3 % 2; 
			
				
			// time tag
			if (next_frame) {
				bk1 = 0; 
				bk2 = 0; 
				bk3 = 0; 
			}
			if (bk1 == 0 && bk2 == 0 && bk3 == 0) {
					
					 
				day1 = din1 >> 19; 
				din1 = din1 - (day1<<19);
				if (day1 > 512) {
					day1 = day1 - 512; 
				}
					 
				year1 = din1 >> 7; 
				din1 = din1 - (year1<<7);
				month1 = din1; 
					
				milli1 = din2 >> 19; 
				din2 = din2 - (milli1 << 19);
				second1 = din2 >> 13; 
				din2 = din2 - (second1 << 13); 
				minute1 = din2 >> 7; 
				din2 = din2 - (minute1<<7); 
				hour1 = din2;  
	
				
				
				if (year1 == year00) {
					ss1 = milli1;  
					//ss1 = second1; 
					ssdiff = ss1 - ss0;
					if (ssdiff<0) {
						//ssdiff = ssdiff + 60;
						ssdiff = ssdiff + 1000; 
					}
					time_s = time_s + ssdiff; 
					ss0 = ss1; 
				
				}
				 
					
				if (time_s >= frame_end || next_frame) {
					//fnum = fnum + 1; 
					 
					ss0 = milli00; 
					ssdiff = 0; 
					time_s = time_s0; 
					next_frame = false; 
					cout << "Cur time: " << hour1 << ":" << minute1 << ":" << second1 << "." << milli1 << "\n";
					cout << infile1.tellg() << "\n"; 

					milli00 = milli1; 
					ss0 = milli00;
					time_s0 = frame_end; 
					time_s = time_s0; 
					
					frame_num++;	
					frame_start = frame_end;
					frame_end = frame_start + frame_length[frame_num];
					if (frame_num % 2 == 1)	{				
						write_lm = true; 
					}
					else {
						write_lm = false;
					} 	
				}
				
			}
	
				
			if (coinc_tag && write_lm) {
				// write output data
				outfile.write(reinterpret_cast<const char *>(&dout1), sizeof(dout1));
				outfile.write(reinterpret_cast<const char *>(&dout2), sizeof(dout2));
			}
			if (!coinc_tag) {
				outfile.write(reinterpret_cast<const char *>(&dout1), sizeof(dout1));
				outfile.write(reinterpret_cast<const char *>(&dout2), sizeof(dout2));
			}

		}
	}


	infile1.close(); 
	outfile.close(); 


	return 0; 

}







































	
