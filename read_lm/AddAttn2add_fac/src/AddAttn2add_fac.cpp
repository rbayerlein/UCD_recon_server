#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

#define BUFFER_SIZE 65536 // = elements per buffer


using namespace std;

int num_ax_crys_per_unit_w_gap = 85;
int num_ax_crys_per_block = 6;
int num_tx_crys_per_block = 7;
int num_crystals_all = 570360; // 679*840
int num_plane_efficiencies = 461041; // 679*679

struct Lut {
	short nv,nu;
}; // block lut

struct ids {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

int main(int argc, char **argv) {

	if (argc != 6 ) {
		cout << "usage: " << argv[0] << " [lm_file_name] [attn_fp_filename] [user_name] [crys_eff] [plane_eff]" << endl; // instead of nc file read in crys eff and plane eff
		cerr << "wrong number of input parameters" << endl;
		exit(1);
	}

	cout << "Starting AddAttn2add_fac executable..." << endl;
	
	// get folder
	string infile_lm = argv[1];
	string slash = "/";
	size_t found = infile_lm.find_last_of(slash);

	string infile_fullpath = infile_lm.substr(0,found);
	cout << "-> opening list-mode files from folder " << infile_fullpath << endl;

	// cut file extention
	size_t found_ext = infile_lm.find_last_of(".lm");
	string fname_raw = infile_lm.substr(0,found_ext-2);
	cout << "found file name without extension: " << fname_raw << endl;

	// get add fac and mul fac files
	stringstream ss_add_fac, ss_mul_fac, ss_add_fac_out, ss_mul_fac_out;
	ss_add_fac << fname_raw << ".add_fac";
	ss_add_fac_out << fname_raw << ".add_fac_out";	// for writing output
	ss_mul_fac << fname_raw << ".mul_fac";
	ss_mul_fac_out << fname_raw << ".mul_fac_out"; // for writing output

	string infile_add_fac = ss_add_fac.str();
	string infile_mul_fac = ss_mul_fac.str();
	string outfile_add_fac = ss_add_fac_out.str();
	string outfile_mul_fac = ss_mul_fac_out.str();

	cout << "-> input and output add_fac and mul_fac files:\n\t" <<
	infile_lm << "\n\t" << infile_add_fac << "\n\t" << infile_mul_fac << "\n\t" << outfile_add_fac << "\n\t" << outfile_mul_fac << endl;

	// get attenuation file
	string infile_attn = argv[2];

	// always use original add_fac file coming from the read_lm executable 
	// check if original file exists
	cout << "-> checking if original add_fac file exists" << endl;
	string infile_add_fac_orig = infile_add_fac;
	infile_add_fac_orig.append(".original");
	ifstream add_fac_tmp(infile_add_fac_orig.c_str());
	if( ! add_fac_tmp.good()){
		cout << "----> creating original by copying *.add_fac file to *.add_fac.original" << endl;
		ifstream  src(infile_add_fac.c_str(), std::ios::binary);
    	ofstream  dst(infile_add_fac_orig.c_str(), std::ios::binary);
	    dst << src.rdbuf();
	    cout << "----> done." << endl;
	}else cout << "----> exists." << endl;
	add_fac_tmp.close();

	// always use original mul_fac file coming from the read_lm executable 
	// check if original file exists
	cout << "-> checking if original mul_fac file exists" << endl;
	string infile_mul_fac_orig = infile_mul_fac;
	infile_mul_fac_orig.append(".original");
	ifstream mul_fac_tmp(infile_mul_fac_orig.c_str());
	if( ! mul_fac_tmp.good()){
		cout << "----> creating original by copying *.mul_fac file to *.mul_fac.original" << endl;
		ifstream  src(infile_mul_fac.c_str(), std::ios::binary);
    	ofstream  dst(infile_mul_fac_orig.c_str(), std::ios::binary);
	    dst << src.rdbuf();
	    cout << "----> done." << endl;
	}else cout << "----> exists." << endl;
	add_fac_tmp.close();


	// open input files
	FILE *pInputFile_lm =  fopen(infile_lm.c_str(), "rb");
	if (pInputFile_lm == NULL) { // check return value of file pointer
		cerr << infile_lm << " cannot be opened." << endl;
		exit(1);
	}	
	// always read from original add_fac file, which only contains randoms corrections
	FILE *pInputFile_add_fac =  fopen(infile_add_fac_orig.c_str(), "rb");
	if (pInputFile_add_fac == NULL) { // check return value of file pointer
		cerr << infile_add_fac_orig << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pInputFile_mul_fac =  fopen(infile_mul_fac_orig.c_str(), "rb");
	if (pInputFile_mul_fac == NULL) { // check return value of file pointer
		cerr << infile_mul_fac_orig << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pOutputFile_add_fac =  fopen(outfile_add_fac.c_str(), "wb");
	if (pOutputFile_add_fac == NULL) { // check return value of file pointer
		cerr << outfile_add_fac << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pOutputFile_mul_fac =  fopen(outfile_mul_fac.c_str(), "wb");
	if (pOutputFile_mul_fac == NULL) { // check return value of file pointer
		cerr << outfile_mul_fac << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pInputFile_attn = fopen(infile_attn.c_str(), "rb");
	if (pInputFile_attn == NULL ) {
		cerr << infile_attn << " cannot be opened." << endl;
		exit(1);
	}

	//here read in crys eff and plane eff instead of nc file.  
	cout <<"-> getting crystal efficiencies from file: " << endl;
	string crys_fullpath = argv[4];
	cout << "\t" << crys_fullpath << endl;
	ifstream crys_read;
	crys_read.open(crys_fullpath.c_str(), ios::in | ios::binary);

	cout <<"-> getting plane efficiencies from file: " << endl;
	string plane_fullpath = argv[5];
	cout << "\t" << plane_fullpath << endl;
	ifstream plane_read;
	plane_read.open(plane_fullpath.c_str(), ios::in | ios::binary);

	vector<float> nc_crys(num_crystals_all); 
	vector<float> nc_plane(num_plane_efficiencies);

	bool use_plane_default = false;
	if (!plane_read) {
		cout << "Could not open plane efficiencies file, use default values (1) \n";
		use_plane_default = true;
		for (int nci = 0; nci < num_plane_efficiencies; ++nci) {
			nc_plane[nci] = 1.0;
		}
	}

	bool use_crys_default = false;
	if (!crys_read) {
		cout << "Could not open crystal efficiencies file, use default values (1) \n";
		use_crys_default = true;
		for (int nci = 0; nci < num_crystals_all; ++nci) {
			nc_crys[nci] = 1.0;
		}
	}

	if (!use_plane_default) {
		// get plane efficiencies
	//	nc_read.seekg( 2079453 * 4, nc_read.beg);
		for (int p = 0;  p < num_plane_efficiencies; ++p) {
			plane_read.read(reinterpret_cast<char*>(&nc_plane[p]), sizeof(float));
		}
		plane_read.close();
	}
	if (!use_crys_default){
		//get crystal efficiencies
	//	nc_read.seekg(3095517 * 4, nc_read.beg);
		for (int c = 0; c < num_crystals_all; c++) {
			crys_read.read(reinterpret_cast<char*>(&nc_crys[c]), sizeof(float));
		}
		crys_read.close();
	}

/*
	// open michelogram LUT
	cout << "-> getting michelogram LUT" << endl;
	string mich_fullpath = "michel_lut_672x672";
	ifstream mich_read;
	mich_read.open(mich_fullpath.c_str(), ios::in | ios::binary);
	int mich[672][672];

	if (!mich_read){
		cout << "Could not open michelogram LUT." << endl;
		exit(1);
	}else{
		for(int m_x = 0; m_x < 672; m_x++){
			for(int m_y = 0; m_y < 672; m_y++){
				mich_read.read(reinterpret_cast<char*>(&mich[m_x][m_y]), sizeof(int));
			}
		}
	}
*/

// ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo ==== main program start

	// read in events from lm file
	cout << "-> reading from list mode file...(may take a while)" << endl;

	// get buffers 
	ids * pids_in = new ids[BUFFER_SIZE];
	float *add_fac_out = new float[BUFFER_SIZE];
	float *add_fac = new float[BUFFER_SIZE];
	float *mul_fac = new float[BUFFER_SIZE];
	float *mul_fac_out = new float[BUFFER_SIZE];
	float *attn = new float[BUFFER_SIZE];

	int num_buffers_read = 1;
	int buffer_indx = 0;
	
	while(!feof(pInputFile_lm)){
		//cout << "----> reading buffer number " << num_buffers_read << endl;

		int read_ct = fread(pids_in, sizeof(ids), BUFFER_SIZE, pInputFile_lm);
		fread(add_fac, sizeof(float), BUFFER_SIZE, pInputFile_add_fac);
		fread(mul_fac, sizeof(float), BUFFER_SIZE, pInputFile_mul_fac);
		fread(attn, sizeof(float), BUFFER_SIZE, pInputFile_attn);
		//cout << "Events in this buffer (read_ct): "<< read_ct << endl;

		for (int i = 0; i < read_ct; ++i)
		{	
			// get crystal IDs
			short txCrysA = pids_in[i].txIDA;
			short txCrysB = pids_in[i].txIDB;
//			short axCrysA = pids_in[i].axIDA - (pids_in[i].axIDA / num_ax_crys_per_unit_w_gap);
//			short axCrysB = pids_in[i].axIDB - (pids_in[i].axIDB / num_ax_crys_per_unit_w_gap);

			short axCrysA = pids_in[i].axIDA;
			short axCrysB = pids_in[i].axIDB;

			// convert to block pids_in
/*			short txBiA = txCrysA / num_tx_crys_per_block;
			short axBiA = axCrysA / num_ax_crys_per_block;
			short txBiB = txCrysB / num_tx_crys_per_block;
			short axBiB = axCrysB / num_ax_crys_per_block;
*/
			// get add_fac from file (i.e. randoms)
			float rtemp = add_fac[i];

			//divide add_fac by mul_fac (apply dead time and decay correction)
			if(mul_fac[i] != 0){
				rtemp /= mul_fac[i];
			}else{
				rtemp = 0;
			}
			//divide add_fac by crystal normalization i.e. multiply by UIH-defined correction factors
			rtemp *= (nc_crys[axCrysA + 679*txCrysA] * nc_crys[axCrysB + 679*txCrysB]);
			//divide add_fac by plane normalization i.e. multiply by UIH-defined correction factors
			rtemp *= (nc_plane[axCrysA + 679*axCrysB]);
			//divide add_fac by attn
			if (attn[i] !=0){ // catch dividing by zero. Should theoretically never happen as attenuation is never zero, even in air
				rtemp /= attn[i];
			}else{
				rtemp = 0;
			}

			add_fac_out[i]=rtemp;
			buffer_indx++;

			// REMOVE THIS LINE!! IT IS JUST FOR DEBUGGING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			//add_fac_out[i] = 0;
			// REMOVE THIS LINE!! IT IS JUST FOR DEBUGGING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			// fill new mul fac file with ones as they are now incorporated into add_fac files.
			mul_fac_out[i] = 1;

			if(buffer_indx == BUFFER_SIZE){
				fwrite(add_fac_out, sizeof(float), BUFFER_SIZE, pOutputFile_add_fac);
				fwrite(mul_fac_out, sizeof(float), buffer_indx, pOutputFile_mul_fac);
				buffer_indx = 0; // reset index
			//if (num_buffers_read % 10 == 0)	cout << "buffer number " << num_buffers_read << "\t" << add_fac_out[i] << "\t" << attn[i] << endl;
			}
		}

		num_buffers_read++;
	}
	fwrite(add_fac_out, sizeof(float), buffer_indx, pOutputFile_add_fac);
	fwrite(mul_fac_out, sizeof(float), buffer_indx, pOutputFile_mul_fac);

	// re-name original ADD FAC file and save new file under original name
	cout << "-> deleting existing *.add_fac and renaming new file to *.add_fac" << endl;

	if ( remove(infile_add_fac.c_str()) != 0 ){
		cerr << "!! could not remove previous file " << infile_add_fac << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully deleted previous file\n" << infile_add_fac << endl;
	}

	string outfile_add_fac_new = outfile_add_fac;
	outfile_add_fac_new.erase(outfile_add_fac.length() -4, 4);
	if (rename(outfile_add_fac.c_str() , outfile_add_fac_new.c_str()) != 0){
		cerr << "!! renaming file\n" << outfile_add_fac << "\nto\n" << outfile_add_fac_new << "\nfailed." << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully renamed file\n" << outfile_add_fac << "\nto\n" << outfile_add_fac_new << endl;
	}	 

	// re-name original MUL FAC file and save new file under original name
	cout << "-> deleting existing *.mul_fac and renaming new file to *.mul_fac" << endl;

	if ( remove(infile_mul_fac.c_str()) != 0 ){
		cerr << "!! could not remove previous file " << infile_mul_fac << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully deleted previous file\n" << infile_mul_fac << endl;
	}

	string outfile_mul_fac_new = outfile_mul_fac;
	outfile_mul_fac_new.erase(outfile_mul_fac.length() -4, 4);
	if (rename(outfile_mul_fac.c_str() , outfile_mul_fac_new.c_str()) != 0){
		cerr << "!! renaming file\n" << outfile_mul_fac << "\nto\n" << outfile_mul_fac_new << "\nfailed." << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully renamed file\n" << outfile_mul_fac << "\nto\n" << outfile_mul_fac_new << endl;
	}	 


	cout << "-> closing all files" << endl;
	delete[] pids_in;
	delete[] add_fac_out;
	delete[] add_fac;
	delete[] mul_fac;
	delete[] attn;
	fclose(pInputFile_lm);
	fclose(pInputFile_mul_fac);
	fclose(pInputFile_add_fac);
	fclose(pOutputFile_add_fac);
	fclose(pInputFile_attn);

	cout << "done" << endl;

}
