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
	bool write_cod = true; 
	
	bool rand_sub = false; 

	// Read config file
	//string scan_details_fullpath = "/run/media/meduser/data/read_lm/Reconstruction_Parameters_1";
	string scan_details_fullpath;
	scan_details_fullpath = argv[1];  

	string roi_image_fullpath; 
	roi_image_fullpath = argv[2]; 


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
			//cout << "Frame " << i << " = " << frame_length_s[i] << " seconds\n";

		}
	}




	// Read roi mask image 
	// Image must be 256 x 256 x 830 (2.344 isotropic)
	ifstream roi_image;
	roi_image.open(roi_image_fullpath.c_str(), ios::in | ios::binary); //open list mode file
	

	if (!roi_image) {
		roi_image.close();
		cout << roi_image_fullpath << "\n"; 
		cout << "Cannot open ROI image\nCheck folder names and locations\n";
		
		return 1;
	} 
		
	roi_image.seekg(0, roi_image.end);
	int img_size = roi_image.tellg() / sizeof(int); //get image size
	
	if (img_size != 256*256*830) {
		cout << "Invalid ROI image size\nSize must be 256 x 256 x 830\n"; 
		return 1; 
	}
	
	
	roi_image.seekg(0, roi_image.beg);

	vector<int> img_mask(256*256*830); 
	int num_roi = 0; 
	for (int r = 0; r<256*256*830; r++) {
		roi_image.read(reinterpret_cast<char *>(&img_mask[r]), sizeof(int));
		if (img_mask[r] > num_roi) {
			num_roi = img_mask[r]; 
		}
	}


	if (num_roi == 0) {
		cout << "mask image needs at least 1 roi\n"; 
		return 1; 
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

	
	
	//int num_crystals_all = 94080; 
	
	// for reading lm file	
	
	uint din1, din2 = 0;
	uint a32 = 2147483648;
	
	int energyA, energyB = 0; 
	int ttA, ttB, tt = 0;
	double ttd, t_v = 0.0;  
	int p_d, coinc = 0; 
	bool coinc_tag = true; 
	bool blockcount_tag = false; 
	bool time_tag = false; 
	bool time_tag0 = false; 
	//bool found_lut = false; 
	int bk1, bk2, bk3 = 0; 
	bool prompt = true; 
	int ua0, ua1, ua2, ub0, ub1, ub2, ua, ub = 0; 
	int crys1, crys2 = 0; 
	int transA, transB, axA, axB = 0; 
	int bankA, bankB, bankk = 0; 
	int modA, modB = 0; 
	int ui0, ui1, ui2, ui3, uid = 0; 
	int unitA, unitB, unit_diff = 0;
	int aa, bb = 0; 
	
	 
	bool corrupted = false; 
	signed short dtemp = 0; 
	short transcA, transcB, crysaxA, crysaxB, dt = 0; 
	short dout[5];
	unsigned short eout [2];  
	
	int e_low = 0; 
	int e_high = 200; 
	
	// time tag variables
	int year1, month1, day1, hour1, minute1, second1, milli1 = 0;
	int year0, month0, day0, hour0, minute0, second0, milli0 = 0;
	int year00, month00, day00, hour00, minute00, second00, milli00 = 0; 
	int ss0, ss1, ssdiff = 0; 
	
	// block count rate variables
	vector<double> block_rate(120*14*8); 
	int blkX, blkXa, blkXb, blkY, blkYa, blkYb, blk_abs, blk_absA, blk_absB = 0; 
	int br_temp = 0; 
	double singles_c1, singles_c2 = 0.0; 
	double rstemp = 0.0; 
	for (int br = 0; br<(120*14*8); br++) {
		block_rate[br] = 0.0; 
	}
	 
	
	// sinogram variables
	long ind, ind2, indind, indind2, ind_block, ind_block2 = 0; 
	int nv, nu, nvtemp, nutemp, ntemp = 0;
	int sino_ax_span = 10000; 
	int block_ax_span = 2; 

	 
	vector<int> sino(num_bins_sino*num_slice); 
	vector<int> sino_r(num_bins_sino*num_slice);

	vector<int> sino_block(num_bins_sino_block*num_slice_block); 
	vector<int> sino_block_r(num_bins_sino_block*num_slice_block);

	// HistoTOF image variables
	int histo_img_size = 128 * 128 * 324; 
	vector<float> histo_image(histo_img_size); 
	double c_xx, c_yy, c_zz, v_xx, v_yy, v_zz, l_v = 0.0; 
	int p_xx, p_yy, p_zz, ind_histo = 0; 
	double img_counts = 0.0; 

	// tof calibration variables
	vector<double> tof_store(num_crystals_all); 
	vector<double> counts_store(num_crystals_all); 
	double phiA, phiB, phiphi, fan_angle; 
	double max_fan_angle = PI / 4.0; 

	// COD motion detect variables
	double tof_res = 450.0;
	double frame_counter = 0.0; 
	vector<double> x_mean(num_roi, 0.0); 
	vector<double> y_mean(num_roi, 0.0);  
	vector<double> z_mean(num_roi, 0.0); 
	vector<double> xyz_mean_counts(num_roi, 0.0); 
	//vector<double> y_mean_counts(8, 0.0);  
	//vector<double> z_mean_counts(8, 0.0); 

	double tof_x, tof_y, tof_z = 0.0; 
	int unit_z = 0; 
	int pix_x, pix_y, pix_z = 0; 
	int img_index = 0; 
	int roi_val = 0; 
	double pix_xd, pix_yd, pix_zd = 0.0; 

	//vector<int> sino(num_bins_sino2*num_slice); 
	//vector<int> sino_r(num_bins_sino2*num_slice);
	
	// randoms
	float rtemp = 0.0;
	int sec1, sec2 = 0;
	int ct1, ct2, rt1, rt2 = 0;
	float lorsum1, lorsum2, lorsumall = 0.0;
	float t_window = 1.0; 
	
	
	
	// dynamic 
	int frame_start = 0.0;
	int frame_end = 0.0;
	int frame_num = 0; 
	int time_s, time_s0 = 0; 
	
	frame_end = frame_start + frame_length[0];

	
	
	// counters
	double cou = 100000000.0;
	double num_prompts, num_randoms, num_coinc = 0;
	
	// to run 
	bool run = true; 
	bool write_lm = false; 
	string runrecon_str = ""; 
	string removeraws_str = ""; 


	//////////////////////////////////////////////////////////////////////////////
	
	// ********************  Read Listmode  *********************************** // 
	
	// open all 8 .raw files
	//string infile_fullpath_cut = infile_fullpath; 
	//string toErase = "1.raw"; 
	//size_t pps = infile_fullpath.find(toErase); 
	//if (pps != string::npos) {
	//	infile_fullpath_cut.erase(pps, toErase.length()); 
	//}
	//cout << infile_fullpath_cut << "\n"; 
	string infile_fullpath1 = infile_fullpath;
	
		
	//infile_fullpath1.append("1.raw");
 
	
	// open listmode file
	ifstream infile1;
	infile1.open(infile_fullpath1.c_str(), ios::in | ios::binary); //open list mode file
	

	if (!infile1) {
		infile1.close();
		cout << infile_fullpath1 << "\n"; 
		cout << "Cannot open listmode file\nCheck folder names and locations\n";
		
		//return 1;
	}
	
	
	
	long long file_size = 0; 
		
	infile1.seekg(0, infile1.end);
	long long file_size1 = infile1.tellg(); //get size (events) of list mode file
	long long num_events = file_size1 / 8;
	cout << num_events << " total events\n"; 
	infile1.seekg(0, infile1.beg);
	

	string fname_out;
	string outfile_fullpath_p = "";
	string outfile_fullpath_r = "";
	string outfile_fullpath_s = ""; 
	string outfile_fullpath_s_nontof = ""; 
	

	//cout << outfolder <<"\n"; 

	stringstream ss;
	//ss << "lm_reorder_f" << frame_num << "_" << frame_start << "s_" << frame_end << "s";
	ss << "lm_reorder_f" << frame_num;
	fname_out = ss.str();
	outfile_fullpath_p = outfolder;
	outfile_fullpath_r = outfolder;
	outfile_fullpath_s = outfolder; 
	
	outfile_fullpath_p.append(fname_out);
	outfile_fullpath_r.append(fname_out);
	outfile_fullpath_s.append(fname_out); 
	
	outfile_fullpath_p.append("_prompts."); 
	outfile_fullpath_p.append(raw_num); 
	outfile_fullpath_p.append(".lm");

	outfile_fullpath_r.append("_randoms.");
	outfile_fullpath_r.append(raw_num); 
	outfile_fullpath_r.append(".lm"); 

	outfile_fullpath_s.append("_sub."); 
	outfile_fullpath_s.append(raw_num); 
	outfile_fullpath_s.append(".lm"); 
	
	ss << "";
	ss.clear();
	str_temp = "";
	fname_out = "";

	//ofstream outfile_p;	//prompts output
	//outfile_p.open(outfile_fullpath_p.c_str(), ios::out | ios::binary); //create binary file containing new crystal + time data

	//ofstream outfile_r;	//randoms lm output
	//outfile_r.open(outfile_fullpath_r.c_str(), ios::out | ios::binary); //create binary file containing new crystal + time data
	
	//ofstream outfile_s; 
	//outfile_s.open(outfile_fullpath_s.c_str(), ios::out | ios::binary); 
	
	
	//outfile_fullpath_p=""; 
	//outfile_fullpath_r="";
	//outfile_fullpath_s=""; 
	
	
	// ssrb sino names
	//string fname_sino;
	//string sino_fullpath_p;
	//string sino_fullpath_r; 


	// tof names

	/*
	string fname_tof_out; 
	string fname_tofcounts_out; 
	string tof_fullpath; 
	string tofcounts_fullpath; 

	tof_fullpath = outfolder; 
	tofcounts_fullpath = outfolder; 
	fname_tof_out = "t_off."; 
	fname_tofcounts_out = "t_off_counts."; 
	tof_fullpath.append(fname_tof_out); 
	tofcounts_fullpath.append(fname_tofcounts_out); 
	tof_fullpath.append(raw_num); 
	tofcounts_fullpath.append(raw_num); 
	tof_fullpath.append(".swt");  
	tofcounts_fullpath.append(".swc"); 

	fname_tof_out = ""; 
	fname_tofcounts_out = ""; 

	ofstream outfile_tof; 
	outfile_tof.open(tof_fullpath.c_str(), ios::out | ios::binary);

	ofstream outfile_tofcounts; 
	outfile_tofcounts.open(tofcounts_fullpath.c_str(), ios::out | ios::binary); 
	*/


	// cod xyz names
	string fname_cod_xyz; 
	string fname_cod_xyz_counts; 
	string outfile_cod_xyz_fullpath;
	string outfile_cod_xyz_counts_fullpath; 

	outfile_cod_xyz_fullpath = outfolder;
	//outfile_cod_xyz_counts_fullpath = outfolder; 
	fname_cod_xyz = "cod_xyz."; 
	//fname_cod_xyz = "cod_xyz_counts."; 
	
	outfile_cod_xyz_fullpath.append(fname_cod_xyz); 
	outfile_cod_xyz_fullpath.append(raw_num);
	outfile_cod_xyz_fullpath.append(".cg");

	//outfile_cod_xyz_counts_fullpath.append(fname_cod_xyz_counts); 
	//outfile_cod_xyz_counts_fullpath.append(raw_num);
	//outfile_cod_xyz_counts_fullpath.append(".cgc");


	ofstream outfile_cod_xyz; 
	outfile_cod_xyz.open(outfile_cod_xyz_fullpath.c_str(), ios::out | ios::binary); 


	//ofstream outfile_cod_xyz_counts; 
	//outfile_cod_xyz_counts.open(outfile_cod_xyz_counts_fullpath.c_str(), ios::out | ios::binary); 


	/*
	stringstream sinoss;
	sinoss << "sinogram_ssrb_f" << (frame_num);
	fname_sino = sinoss.str();
	
	sino_fullpath_p = outfolder;
	sino_fullpath_r = outfolder; 
	
	sino_fullpath_p.append(fname_sino);
	sino_fullpath_p.append("_prompts."); 
	sino_fullpath_p.append(raw_num); 
	sino_fullpath_p.append(".raw");
	
	sino_fullpath_r.append(fname_sino); 
	sino_fullpath_r.append("_randoms."); 
	sino_fullpath_r.append(raw_num); 
	sino_fullpath_r.append(".raw"); 
	sinoss << "";
	sinoss.clear();
	str_temp = "";
	fname_sino = "";
	ofstream sinogram;
	sinogram.open(sino_fullpath_p.c_str(), ios::out | ios::binary); //create binary file containing new	
	
	ofstream sinogram_r;
	sinogram_r.open(sino_fullpath_r.c_str(), ios::out | ios::binary); //create binary file containing new	
	
	sino_fullpath_p = "";
	sino_fullpath_r = ""; 
	*/
	
	// block sino names
	//string fname_sino_block;
	//string sino_fullpath_block_p;
	//string sino_fullpath_block_r; 

	string fname_histo_img; 
	string histo_img_fullpath; 

	//string fname_tof_cal; 
	//string tof_cal_fullpath; 
	//string tof_cal_counts_fullpath;


	
	/*
	stringstream sinobss;
	sinobss << "sinogram_block_f" << (frame_num);
	fname_sino_block = sinobss.str();
	
	sino_fullpath_block_p = outfolder;
	sino_fullpath_block_r = outfolder; 
	
	sino_fullpath_block_p.append(fname_sino_block);
	sino_fullpath_block_p.append("_prompts."); 
	sino_fullpath_block_p.append(raw_num); 
	sino_fullpath_block_p.append(".raw");
	
	sino_fullpath_block_r.append(fname_sino_block); 
	sino_fullpath_block_r.append("_randoms."); 
	sino_fullpath_block_r.append(raw_num); 
	sino_fullpath_block_r.append(".raw"); 
	sinobss << "";
	sinobss.clear();
	str_temp = "";
	fname_sino_block = "";
	ofstream sinogram_block;
	sinogram_block.open(sino_fullpath_block_p.c_str(), ios::out | ios::binary); //create binary file containing new	
	
	ofstream sinogram_block_r;
	sinogram_block_r.open(sino_fullpath_block_r.c_str(), ios::out | ios::binary); //create binary file containing new	
	
	sino_fullpath_block_p = "";
	sino_fullpath_block_r = ""; 
	
	*/


	//string frameinfo_fullpath;
	
	
	//cout << outfolder <<"\n";


	// ************		Load LUTs  ****************//
	string fdir_code = "../read_lm/lut/"; 
	
	// open bank lut
	int lutsum = 0;
	string bank_lut_path = fdir_code; 
	bank_lut_path.append("bank_lut");
	//cout << bank_lut_path << "\n";
	ifstream lut_read; 
	lut_read.open(bank_lut_path.c_str(), ios::in | ios::binary); 
	if (!lut_read) {
		cout << "could not open bank LUT file\n";
		return 1;  
	}
	int bank_lut[54][2]; 
	for (int ii=0; ii<108; ii++) {
		if (ii<54) {
			lut_read.read(reinterpret_cast<char *>(&bank_lut[ii][0]), sizeof(int));
			lutsum = lutsum + bank_lut[ii][0];
		}
		else {
			lut_read.read(reinterpret_cast<char *>(&bank_lut[ii-54][1]), sizeof(int));
			lutsum = lutsum + bank_lut[ii-54][1]; 
		}
	}
	
	if (lutsum != 594) {
		cout << "Bank LUT incorrect!\n";
		return 1; 
	} 
	

	// crystal index - sinogram bin LUTs
	vector<int> index_crystalpairs_transaxial_int16_1(num_bins_sino);
	vector<int> index_crystalpairs_transaxial_int16_2(num_bins_sino);
	
	string fpath_crysidx = fdir_code; 
	fpath_crysidx.append("index_crystalpairs_transaxial_2x549x420_int16");  
	ifstream crysidx_file;
	crysidx_file.open(fpath_crysidx.c_str(), ios::in | ios::binary); //open list mode file

	if (!crysidx_file) {
		crysidx_file.close();
		cout << "Cannot open crys_idx file\nQUIT\n\n";
		
		return 1;
	}

	for (int ns = 0; ns < num_bins_sino; ns++) {

		crysidx_file.read(reinterpret_cast<char *>(&index_crystalpairs_transaxial_int16_1[ns]), sizeof(signed short));
		crysidx_file.read(reinterpret_cast<char *>(&index_crystalpairs_transaxial_int16_2[ns]), sizeof(signed short));
	}

	

	int noindex_crystalpairs_transaxial_int16[840][840];

	for (int jj = 0; jj < num_crystals; jj++) {
		for (int kk = 0; kk < num_crystals; kk++) {
			noindex_crystalpairs_transaxial_int16[jj][kk] = -1;
		}
	}

	

	for (int nss = 0; nss < num_bins_sino; nss++) {
		nv = index_crystalpairs_transaxial_int16_1[nss];
		nu = index_crystalpairs_transaxial_int16_2[nss];
		//cout << nv << " " << nu << "\n"; 
		noindex_crystalpairs_transaxial_int16[nv][nu] = nss;
		noindex_crystalpairs_transaxial_int16[nu][nv] = nss;
	}





	// block index - block sinogram bin LUTs
	vector<int> index_blockpairs_transaxial_int16_1(num_bins_sino_block);
	vector<int> index_blockpairs_transaxial_int16_2(num_bins_sino_block);
	
	string fpath_blockidx = fdir_code; 
	fpath_blockidx.append("index_blockpairs_transaxial_2x91x60_int16");  
	ifstream blockidx_file;
	blockidx_file.open(fpath_blockidx.c_str(), ios::in | ios::binary); //open list mode file

	if (!blockidx_file) {
		blockidx_file.close();
		cout << "Cannot open block_idx file\nQUIT\n\n";
		
		return 1;
	}

	for (int nbs = 0; nbs < num_bins_sino_block; nbs++) {

		blockidx_file.read(reinterpret_cast<char *>(&index_blockpairs_transaxial_int16_1[nbs]), sizeof(signed short));
		blockidx_file.read(reinterpret_cast<char *>(&index_blockpairs_transaxial_int16_2[nbs]), sizeof(signed short));
	}

	

	int noindex_blockpairs_transaxial_int16[120][120];

	for (int jbj = 0; jbj < 120; jbj++) {
		for (int kbk = 0; kbk < 120; kbk++) {
			noindex_blockpairs_transaxial_int16[jbj][kbk] = -1;
		}
	}

	

	for (int nbss = 0; nbss < num_bins_sino_block; nbss++) {
		nv = index_blockpairs_transaxial_int16_1[nbss];
		nu = index_blockpairs_transaxial_int16_2[nbss];
		//cout << nv << " " << nu << "\n"; 
		noindex_blockpairs_transaxial_int16[nv][nu] = nbss;
		noindex_blockpairs_transaxial_int16[nu][nv] = nbss;
	}



	// crystal position x,y,z
	string rx_lut = fdir_code;
	string ry_lut = fdir_code;
	string zz_lut = fdir_code; 
	rx_lut.append("rx");
	ry_lut.append("ry");
	zz_lut.append("zz"); 

	ifstream rx_read; 
	ifstream ry_read; 
	ifstream zz_read; 

	rx_read.open(rx_lut.c_str(), ios::in | ios::binary); 
	ry_read.open(ry_lut.c_str(), ios::in | ios::binary); 
	zz_read.open(zz_lut.c_str(), ios::in | ios::binary); 

	vector<double> rx(num_crystals); 
	vector<double> ry(num_crystals); 
	vector<double> zz(num_crystals_ax_wgap); 

	for (int rr = 0; rr < num_crystals; rr++) { 
		rx_read.read(reinterpret_cast<char *>(&rx[rr]), sizeof(double));
		ry_read.read(reinterpret_cast<char *>(&ry[rr]), sizeof(double));
	}

	for (int zzz = 0; zzz < num_crystals_ax_wgap; zzz++) { 
		zz_read.read(reinterpret_cast<char *>(&zz[zzz]), sizeof(double));
	}	
	rx_read.close(); 
	ry_read.close(); 
	zz_read.close(); 
	
	
	
	// timing LUTs
	string tLUT_temp = infile_fullpath1; 
	int str_len = tLUT_temp.length(); 
	tLUT_temp.erase(str_len-3, 3); 
	
	string tLUT_hw = tLUT_temp; 
	string tLUT_sw = tLUT_temp; 
	//string tLUT_tdc = tLUT_temp;
	
	tLUT_hw.append("hw"); 
	tLUT_sw.append("swt"); 
	//tLUT_tdc.append("tdc");  
	
	
	vector<double> t_sw(num_crystals_all); 
	vector<int> t_hw(num_crystals_all); 
	//vector<double> t_tdc(168*560*128); 
	
	ifstream t_sw_read; 
	ifstream t_hw_read; 
	t_sw_read.open(tLUT_sw.c_str(), ios::in | ios::binary); 
	t_hw_read.open(tLUT_hw.c_str(), ios::in | ios::binary); 
	
	if (!t_sw_read) {
		cout << "could not open timing offset files \n"; 
	}
	
	for (int tt1 = 0; tt1<(num_crystals_all); tt1++) {
		t_sw_read.read(reinterpret_cast<char *>(&t_sw[tt1]), sizeof(double));
		t_hw_read.read(reinterpret_cast<char *>(&t_hw[tt1]), sizeof(int)); 
		
		//cout << t_hw[tt1] << "\n"; 
	}

	//cout << outfolder <<"\n";
	
	
	//ifstream t_tdc_read; 
	//t_tdc_read.open(tLUT_tdc.c_str(), ios::in | ios::binary); 
	
	//if (!t_tdc_read) {
		//cout << "could not open TDC file\n"; 
	//}
	
	//for (int tt2 = 0; tt2<(168*560*128); tt2++) { 
	//	t_tdc_read.read(reinterpret_cast<char *>(&t_tdc[tt2]), sizeof(double)); 
	//}
	
	
	// timing offset file
	//string TO_SW_path = fdir_code; 
	//TO_SW_path.append("TO_SW");
	// //cout << TO_SW_path << "\n";
	//ifstream TOSW_read; 
	//TOSW_read.open(TO_SW_path.c_str(), ios::in | ios::binary);
	//vector<double> toff_sw(num_crystals_all); 
	//if (!TOSW_read) {
	//	cout << "could not open timing offset file \n"; 
	//}
	//for (int tt=0; tt<num_crystals_all; tt++) {
	//	TOSW_read.read(reinterpret_cast<char *>(&toff_sw[tt]), sizeof(double));
	//}
	
	vector<float> coinc_window(8); 
	coinc_window[0] = 4500.0; 
	coinc_window[1] = 4800.0;
	coinc_window[2] = 5400.0;
	coinc_window[3] = 6100.0;
	coinc_window[4] = 6900.0;
	coinc_window[5] = 6900.0;
	coinc_window[6] = 6900.0;
	coinc_window[7] = 6900.0;
	coinc_window[8] = 6900.0;
	
	
	
	//////////////////////////////////////////////////////////////////////////////

	// **************		Main Run Program		*********************//
// *****************************************************************************************************************************************************************************
	
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
	
	//return 0; 
	
	
	
	
	
	
	cout << "\n\n*** Processing Frame 0... *** \n\n"; 
	
 	bool next_frame = false; 
 	time_s = time_s0; 
 
	while (run) {
		//event_counter++;
		//tag = false; 
		
		prompt = true;
		
		dout[0] = 0; 
		dout[1] = 0; 
		dout[2] = 0; 
		dout[3] = 0; 
		dout[4] = 0; 
		
		
		
		
		infile1.read(reinterpret_cast<char *>(&din1), sizeof(unsigned int));
		infile1.read(reinterpret_cast<char *>(&din2), sizeof(unsigned int));
		if (infile1.tellg() > (file_size1-100)) {
			//fnum = fnum + 1; 
			next_frame = true; 
		}
		
		
		//infile1.read(reinterpret_cast<char *>(&din1), sizeof(uint));	//read 64 bits of data
		//infile1.read(reinterpret_cast<char *>(&din2), sizeof(uint));

		
		
		//if (infile8.tellg() > (file_size-100)) {	 
			//run = false; 
			//cout << "\nReached end of files\nnum coinc = " << num_coinc << "\n\n"; 	
		//}
		
				
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
					write_lm = true; 		
					
				}					
			}
	
				
			if (coinc_tag) {
				
				if (num_coinc > cou) {
					cout << "num events = " << num_coinc << "\n"; 
					cou = cou + 100000000.0;
					 
				}
				num_coinc = num_coinc+1.0; 
					
				p_d = din2 >> 4; 
				p_d = p_d % 2; 
				if (p_d == 1) { 
					prompt = false;
				}
				else {
					prompt = true; 
				}
					
				ua2 = din1 >> 30; 
				ub2 = din2 >> 30; 
					
				din1 = din1 - (ua2<<30); 
				din2 = din2 - (ub2<<30); 
				
				energyA = din1 >> 23; 
				energyB = din2 >> 23;
				
				din1 = din1 - (energyA<<23);					
				din2 = din2 - (energyB<<23);
				
				ttA = (din1 >> 19);
				ttB = (din2 >> 19); 
					
				tt = ttA + (16*ttB); 
					 
				// two's complement for delta T
				if (tt > 127) {
					tt = tt - 256; 
				}
				
				din1 = din1 - (ttA<<19);
				din2 = din2 - (ttB<<19);
				
				ua0 = din1 >> 3; 
				ub0 = din2 >> 3; 
				ua1 = din1 >> 18; 
				ub1 = din2 >> 18; 
				
				ua0 = ua0 % 2; 
				ub0 = ub0 % 2;
				ua1 = ua1 % 2; 
				ub1 = ub1 % 2; 
				
				ua = ua0 + (2*ua1) + (4*ua2);
				ub = ub0 + (2*ub1) + (4*ub2);
				
				//ua = ua0 + (2*ua1); 
				//ub = ub0 + (2*ub1); 
					
				unitA = ua; 
				unitB = ub; 
				
				//unitA = unitid[ua]; 
				//unitB = unitid[ub]; 
				
				din1 = din1 - (ua1<<18);
				din2 = din2 - (ub1<<18);
				
				crys1 = (din1 >> 5); 
				crys2 = (din2 >> 5); 
				
				din1 = din1 - (crys1<<5);
				din2 = din2 - (crys2<<5); 
				    
				din1 = din1 - (coinc<<4); 
				din2 = din2 - (p_d<<4); 
				
				din1 = din1 - (ua0<<3); 
				din2 = din2 - (ub0<<3);
				
				bankA = din1;
				bankB = din2; 
				
				bankk = bankA + (8*bankB); 
				bankk = bankk-1; 
				
				transA = crys1 % 70; 
				transB = crys2 % 70; 
				modA = bank_lut[bankk][0];
				modB = bank_lut[bankk][1]; 
					
				transA = transA + (modA*70); 
				transB = transB + (modB*70); 
				
				axA = floor(crys1 / 70); 
				axB = floor(crys2 / 70);
				
				axA = axA + (unitA*84); 
				axB = axB + (unitB*84); 
				
				//blkXa = floor(transA / 7); 
				//blkXb = floor(transB / 7); 
				//blkYa = floor(axA / 6); 
				//blkYb = floor(axB / 6); 
				//blk_absA = blkXa + (120*blkYa); 
				//blk_absB = blkXb + (120*blkYb); 
				
				transcA = (signed short)transA; 
				transcB = (signed short)transB; 
				crysaxA = (signed short)axA; 
				crysaxB = (signed short)axB;
					
				 
				// add timing offset to dt; 
				
				crys1 = transcA + (840*crysaxA); 
				crys2 = transcB + (840*crysaxB);
				
				//crys1 = crysaxA + (transcA*672); 
				//crys2 = crysaxB + (transcB*672); 
				
				//tt = tt + t_hw[crys1] - t_hw[crys2]; 
				
				ttd = (double)tt; 
				ttd = ttd*39.0625; 
				 
				//ttd = ttd - toff_sw[crys1] + toff_sw[crys2];
				ttd = ttd - t_sw[crys1] + t_sw[crys2]; 
				t_v = ttd; 
				ttd = round(ttd / 39.0625); 
					
				dt = (signed short)ttd; 
				
				dt = (signed short)tt; 
				
				crysaxA = crysaxA + (short)unitA; 
				crysaxB = crysaxB + (short)unitB; 
				
				dout[0] = transcA; 
				dout[1] = crysaxA; 
				dout[2] = transcB; 
				dout[3] = crysaxB; 
				dout[4] = dt; 
					
				eout[0] = (unsigned short)energyA; 
				eout[1] = (unsigned short)energyB;		
				
				unit_diff = abs(unitA - unitB);
				//t_window = coinc_window[unit_diff]; 
					
					
					
				ind = noindex_crystalpairs_transaxial_int16[transcA][transcB];
				//indind = (ind % 549) + (549 * floor(ind/(549*6)));  

				//ind_block = noindex_blockpairs_transaxial_int16[blkXa][blkXb]; 
				
				
					 
				if (ind>0) {
					//ind2 = ind + (num_bins_sino*(axA + axB)); 
					//indind2 = indind + (num_bins_sino2*(axA + axB)); 

					//ind_block2 = ind_block + (num_bins_sino_block*(blkYa + blkYb)); 
				    
					if (prompt) {

						
						//Put event in sinogram			
						//if (abs(crysaxA - crysaxB) <= sino_ax_span) { 
							//sino[ind2]++;
							//sino[indind2]++; 
						//}
						//if ( (ind_block > 0) && (abs(blkYa - blkYb) <= block_ax_span) ) {
							//sino_block[ind_block2]++; 
						//}
						if (write_histo_img) {
							c_xx = ((double)rx[transcA] + (double)rx[transcB]) / 2.0; 
							c_yy = ((double)ry[transcA] + (double)ry[transcB]) / 2.0;
							c_zz = ((double)zz[crysaxA] + (double)zz[crysaxB]) / 2.0;
							
							v_xx = ((double)rx[transcB] - (double)rx[transcA]);
							v_yy = ((double)ry[transcB] - (double)ry[transcA]);
							v_zz = ((double)zz[crysaxB] - (double)zz[crysaxA]); 
							l_v = sqrt(v_xx*v_xx + v_yy*v_yy + v_zz*v_zz); 
							v_xx = v_xx / l_v; 
							v_yy = v_yy / l_v; 
							v_zz = v_zz / l_v;

							c_xx = c_xx + (t_v * v_xx);
							c_yy = c_yy + (t_v * v_yy);
							c_zz = c_zz + (t_v * v_zz); 

							p_xx = round(c_xx / 4.688); 
							p_yy = round(c_yy / 4.688); 
							p_zz = round(c_zz / 6.0); 
							
							if ( (abs(p_xx) < 64) && (abs(p_yy) < 64) && (p_zz < 324) && (p_zz > 0) ) {
								p_xx = p_xx + 64; 
								p_yy = p_yy + 64;
								//p_zz = p_zz + 324; 
								ind_histo = p_yy + (128*p_xx) + (128*128*p_zz); 								
								histo_image[ind_histo] = histo_image[ind_histo] + 1.0; 
								img_counts = img_counts + 1.0; 
							}
							
						}
						if (write_cod) {
							c_xx = ((double)rx[transcA] + (double)rx[transcB]) / 2.0; 
							c_yy = ((double)ry[transcA] + (double)ry[transcB]) / 2.0;
							c_zz = ((double)zz[crysaxA] + (double)zz[crysaxB]) / 2.0;
							
							v_xx = ((double)rx[transcB] - (double)rx[transcA]);
							v_yy = ((double)ry[transcB] - (double)ry[transcA]);
							v_zz = ((double)zz[crysaxB] - (double)zz[crysaxA]); 
							l_v = sqrt(v_xx*v_xx + v_yy*v_yy + v_zz*v_zz); 
							v_xx = v_xx / l_v; 
							v_yy = v_yy / l_v; 
							v_zz = v_zz / l_v;

							tof_x = (0.3 * t_v / 2.0) * v_xx;
							tof_y = (0.3 * t_v / 2.0) * v_yy;
							tof_z = (0.3 * t_v / 2.0) * v_zz; 

							c_xx = c_xx + tof_x;
							c_yy = c_yy + tof_y;
							c_zz = c_zz + tof_z; 

							// check if event is inside mask image
							pix_xd = 128.0 + (c_xx / 2.344); 
							pix_yd = 128.0 + (-1.0 * c_yy / 2.344);
							pix_zd = (c_zz / 2.344);

							if (pix_xd < 0) {
								pix_xd  = 0.0; 
							}
							if (pix_yd < 0) {
								pix_yd  = 0.0; 
							}
							if (pix_zd < 0) {
								pix_zd  = 0.0; 
							}

							if (pix_xd > 255.0) {
								pix_xd  = 255.0; 
							}
							if (pix_yd > 255.0) {
								pix_yd  = 255.0; 
							}
							if (pix_zd > 829.0) {
								pix_zd  = 829.0; 
							}

							pix_x = round(pix_xd); 
							pix_y = round(pix_yd); 
							pix_z = round(pix_zd); 

							img_index = pix_y + (256*pix_x) + (256*256*pix_z); 

							roi_val = img_mask[img_index]; 


							//unit_z = floor(c_zz / 243.0); 
							//if (unit_z > 7) { 
							//	unit_z = 7; 
							//}
							//if (unit_z < 0) { 
							//	unit_z = 0; 
							//}
							
							unit_z = roi_val - 1; 

							if (roi_val > 0) {
								x_mean[unit_z] = x_mean[unit_z] + c_xx;
								y_mean[unit_z] = y_mean[unit_z] + c_yy;
								z_mean[unit_z] = z_mean[unit_z] + c_zz; 

								xyz_mean_counts[unit_z] = xyz_mean_counts[unit_z] + 1.0;
							}
							
							//y_mean_counts[unit_z] = y_mean_counts[unit_z] + 1.0;
							//z_mean_counts[unit_z] = z_mean_counts[unit_z] + 1.0;  




						}

						if (write_lmfile) {
							//for (int yy=0; yy<5; yy++) {
							//	outfile_p.write(reinterpret_cast<const char *>(&dout[yy]), sizeof(dout[yy]));
							//}
						    
							//singles_c1 = block_rate[blk_absA] / 42.0; 
							//singles_c2 = block_rate[blk_absB] / 42.0; 
							//rstemp = singles_c1 * singles_c2 * frame_length_s[frame_num]; 
							//if (isnan(rstemp) || isinf(rstemp)) {
							//	rstemp = 0.0; 
							//}
							//rtemp = (float)rstemp; 
							//rtemp = rtemp * (39.0625E-12); 
							//outfile_s.write(reinterpret_cast<const char *>(&rtemp), sizeof(float));
							
						}		
						num_prompts = num_prompts+1.0; 
					}
					
					 
					if (!prompt) {
						//Put event in sinogram
						//if (abs(crysaxA - crysaxB) <= sino_ax_span) { 
							//sino_r[ind2]++;
							//sino_r[indind2]++; 
						//}
						//if ( (ind_block > 0) && (abs(blkYa - blkYb) <= block_ax_span) ) {
							//sino_block_r[ind_block2]++; 
						//}
						//if (write_lmfile) {
						//	for (int yy=0; yy<5; yy++) {
						//		outfile_r.write(reinterpret_cast<const char *>(&dout[yy]), sizeof(dout[yy]));
						//	}
						//}			
						num_randoms = num_randoms+1.0; 
					}
				}		
			}  // coincidence
				
		}  // data good
			
			
		if (write_lm == true || run == false || infile1.eof()) {
			
			time_tag = false; 
			
			/*
			frameinfo_fullpath = ""; 
			fname_out = ""; 
			stringstream ss;
			ss << "lm_info_f" << (frame_num);
			fname_out = ss.str(); 
			frameinfo_fullpath = outfolder;
			//cout << outfolder << "\n";
			frameinfo_fullpath.append(fname_out);
			ofstream frameinfo_out;
			frameinfo_out.open(frameinfo_fullpath.c_str());
			//cout << frameinfo_fullpath << "\n";
			if (!frameinfo_out) {
				cout << "Could not create frame info file\n";
			}
			ss << "";
			ss.clear();
			fname_out = "";

			//stringstream ss;
			ss << "startDateTime=" << year00 << "," << month00 << "," << day00 << "," << hour00 << "," << minute00 << "," << second00 << "\n" << "endDateTime=" << year1 << "," << month1 << "," << day1 << "," << hour1 << "," << minute1 << "," << second1 << "\n" << "frame_start=" << (frame_start) << "\n" << "frame_length=" << (frame_end - frame_start) << "\n" << infile1.tellg() << "\n"; // << infile2.tellg() << "\n" << infile3.tellg() << "\n" << infile4.tellg() << "\n" << infile5.tellg() << "\n" << infile6.tellg() << "\n" << infile7.tellg() << "\n" << infile8.tellg() << "\n";
			//ssinfo_out << "frame_length=" << frame_length_s[frame_num] << "\n" << "frame_start=" << (frame_start/100.0) << "\n";
			str_temp = ss.str();
			frameinfo_out << str_temp;
			ss << "";
			ss.clear();
			str_temp = "";
			frameinfo_out.close();
			*/




			//for (int ti = 0; ti < num_crystals_all; ti++) { 
			//	outfile_tof.write(reinterpret_cast<const char *>(&tof_store[ti]), sizeof(double)); 
			//	outfile_tofcounts.write(reinterpret_cast<const char *>(&counts_store[ti]), sizeof(double)); 
			//}
			//outfile_tof.close(); 
			//outfile_tofcounts.close(); 

			for (int ci = 0; ci < num_roi; ci++) {
				outfile_cod_xyz.write(reinterpret_cast<const char *>(&x_mean[ci]), sizeof(double));
				outfile_cod_xyz.write(reinterpret_cast<const char *>(&y_mean[ci]), sizeof(double));
				outfile_cod_xyz.write(reinterpret_cast<const char *>(&z_mean[ci]), sizeof(double));

				outfile_cod_xyz.write(reinterpret_cast<const char *>(&xyz_mean_counts[ci]), sizeof(double));
				//outfile_cod_xyz_counts.write(reinterpret_cast<const char *>(&y_mean_counts[ci]), sizeof(double));
				//outfile_cod_xyz_counts.write(reinterpret_cast<const char *>(&z_mean_counts[ci]), sizeof(double));

				x_mean[ci] = 0.0;
				y_mean[ci] = 0.0;
				z_mean[ci] = 0.0; 

				xyz_mean_counts[ci] = 0.0;
				//y_mean_counts[ci] = 0.0;
				//z_mean_counts[ci] = 0.0; 
			}
			

			


			/* 
			fname_sino = ""; 
			sino_fullpath_p = "";
			sino_fullpath_r = "";
					
			stringstream sinoss;
			sinoss << "sinogram_ssrb_f" << (frame_num);
			fname_sino = sinoss.str();
			sino_fullpath_p = outfolder;
			sino_fullpath_r = outfolder; 
					
			sino_fullpath_p.append(fname_sino);
			sino_fullpath_p.append("_prompts."); 
			sino_fullpath_p.append(raw_num); 
			sino_fullpath_p.append(".raw");

			sino_fullpath_r.append(fname_sino); 
			sino_fullpath_r.append("_randoms."); 
			sino_fullpath_r.append(raw_num); 
			sino_fullpath_r.append(".raw"); 

			sinoss << "";
			sinoss.clear();
			str_temp = "";
			fname_sino = "";

			ofstream sinogram; 			
			sinogram.open(sino_fullpath_p.c_str(), ios::out | ios::binary); //create binary file containing new	
	
			ofstream sinogram_r; 
			sinogram_r.open(sino_fullpath_r.c_str(), ios::out | ios::binary); //create binary file containing new	
	
			sino_fullpath_p = "";
			sino_fullpath_r = ""; 
			
			
			
			// write sinograms
			if (write_sino) {
				for (int si = 0; si < (num_bins_sino * num_slice); si++) {
				
					//sino[si] = sino[si] - sino_r[si]; 
				
				
					sinogram.write(reinterpret_cast<const char *>(&sino[si]), sizeof(sino[si]));
					sinogram_r.write(reinterpret_cast<const char *>(&sino_r[si]), sizeof(sino_r[si]));
				
				
					sino[si] = 0; 
					sino_r[si] = 0; 
				}
			}
			sinogram.close();
			sinogram_r.close(); 



			// write block sinograms
			// block sino names
			fname_sino_block = ""; 
			sino_fullpath_block_p = ""; 
			sino_fullpath_block_r = ""; 

			stringstream sinobss;
			sinobss << "sinogram_block_f" << (frame_num);
			fname_sino_block = sinobss.str();
	
			sino_fullpath_block_p = outfolder;
			sino_fullpath_block_r = outfolder; 
	
			sino_fullpath_block_p.append(fname_sino_block);
			sino_fullpath_block_p.append("_prompts."); 
			sino_fullpath_block_p.append(raw_num); 
			sino_fullpath_block_p.append(".raw");
	

			sino_fullpath_block_r.append(fname_sino_block); 
			sino_fullpath_block_r.append("_randoms."); 
			sino_fullpath_block_r.append(raw_num); 
			sino_fullpath_block_r.append(".raw"); 

			sinobss << "";
			sinobss.clear();
			str_temp = "";
			fname_sino_block = "";
			
			
			ofstream sinogram_block;
			sinogram_block.open(sino_fullpath_block_p.c_str(), ios::out | ios::binary); 

			ofstream sinogram_block_r;
			sinogram_block_r.open(sino_fullpath_block_r.c_str(), ios::out | ios::binary); 	
	
			sino_fullpath_block_p = "";
			sino_fullpath_block_r = ""; 
			if (write_sino_block) {
				for (int sbi = 0; sbi < (num_bins_sino_block * num_slice_block); sbi++) {
				
					//sino[si] = sino[si] - sino_r[si]; 
				
				
					sinogram_block.write(reinterpret_cast<const char *>(&sino_block[sbi]), sizeof(sino_block[sbi]));
					sinogram_block_r.write(reinterpret_cast<const char *>(&sino_block_r[sbi]), sizeof(sino_block_r[sbi]));
					sino_block[sbi] = 0; 
					sino_block_r[sbi] = 0;
				}
				
				 
			}
			sinogram_block.close();
			sinogram_block_r.close();

			histo_img_fullpath = "";
			fname_histo_img = "";
			stringstream histo_img_ss;
			histo_img_ss << "histo_tof_img_f" << (frame_num) << ".";
			fname_histo_img = histo_img_ss.str();
			histo_img_fullpath = outfolder; 
			histo_img_fullpath.append(fname_histo_img); 
			histo_img_fullpath.append(raw_num); 
			histo_img_fullpath.append(".img"); 
			histo_img_ss << ""; 
			histo_img_ss.clear(); 
			
			ofstream histo_img_write; 
			histo_img_write.open(histo_img_fullpath.c_str(), ios::out | ios::binary); 
			
			if (write_histo_img) { 
				for (int hii = 0; hii < histo_img_size; hii++) {
					histo_img_write.write(reinterpret_cast<const char *>(&histo_image[hii]), sizeof(float));
					histo_image[hii] = 0.0; 
				}
			}
			histo_img_write.close();  

			cout << "histo image counts = " << img_counts << "\n"; 

			*/

			// close output lm files
			//outfile_p.close();
			//outfile_r.close();
			//outfile_s.close(); 
			
			

			
				
			

			// Create files for next frame, or quit if end frame is reached

			frame_num++;
			if (frame_num == tot_frames || infile1.eof()) {

				run = false;
				//infile1.close();
				cout << "Scan time: " << (time_s/1000) << " seconds\n";
				cout << "\nFinish listmode processing and sinogram binning!\n";
				//return 1;
			}

			else {
				cout << "\n\n*** Processing Frame " << frame_num << "... ***\n\n";
				frame_start = frame_end;
				frame_end = frame_start + frame_length[frame_num];
				if (write_lmfile || write_sino) {

					/*					
					outfile_fullpath_p = "";
					outfile_fullpath_r = "";
					outfile_fullpath_s = ""; 
					

					stringstream ss;
					//ss << "lm_reorder_frame" << frame_num << "_" << frame_start << "s_" << frame_end << "s";
					ss << "lm_reorder_f" << frame_num;
					fname_out = ss.str();
					outfile_fullpath_p = outfolder;
					outfile_fullpath_r = outfolder;
					outfile_fullpath_s = outfolder; 
					
					outfile_fullpath_p.append(fname_out);
					outfile_fullpath_r.append(fname_out);
					outfile_fullpath_s.append(fname_out); 
					
					outfile_fullpath_p.append("_prompts."); 
					outfile_fullpath_p.append(raw_num); 
					outfile_fullpath_p.append(".lm");

					outfile_fullpath_r.append("_randoms."); 
					outfile_fullpath_r.append(raw_num); 
					outfile_fullpath_r.append(".lm");

					outfile_fullpath_s.append("_sub."); 
					outfile_fullpath_s.append(raw_num); 
					outfile_fullpath_s.append(".lm"); 
					
					ss << "";
					ss.clear();
					str_temp = "";
					fname_out = "";

					outfile_p.open(outfile_fullpath_p.c_str(), ios::out | ios::binary);
					outfile_r.open(outfile_fullpath_r.c_str(), ios::out | ios::binary);
					outfile_s.open(outfile_fullpath_s.c_str(), ios::out | ios::binary); 
					
					*/


					

				}
			
				num_prompts = 0;
				num_randoms = 0;
				write_lm = false;
			}
		}
		
	}

	
		
	infile1.close();
	
	
	//if (write_lmfile) {
		//outfile_p.close();
		//outfile_r.close();
		//outfile_s.close(); 
		
	//}

	outfile_cod_xyz.close(); 
	//outfile_cod_xyz_counts.close(); 
	
	
	//cout << "Scan time: " << (time_s) << " seconds\n";
	
	//cout << "\nListmode pre-processing completed \n";

	
	return 1;



}
