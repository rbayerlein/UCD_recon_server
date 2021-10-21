#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

#define BUFFER_SIZE 65536 // = elements per buffer


using namespace std;

	// see what of the following you can delete in the end! FIXME

int num_bins_sino = 549 * 420;
int num_bins_sino_block = 91 * 60;
int num_ax_block_ring = 112; // 14*8
int num_tx_block_ring = 120; // 5*24
int num_ax_crys_per_unit_w_gap = 85;
int num_ax_crys_per_block = 6;
int num_tx_crys_per_block = 7;
int num_lor_blkpair = 903;
int num_crystals_all = 564480; // 672*840
int num_plane_efficiencies = 451584; // 672*672

struct Lut {
	short nv,nu;
}; // block lut

struct ids {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

int main(int argc, char **argv) {

	if (argc != 5 ) {
		cout << "usage: " << argv[0] << " [lm_file_name] " << " [attn_fp_filename] " << " [user_name] " << " [nc_file]" << endl;
		cerr << "not enough input parameters" << endl;
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
	stringstream ss_add_fac, ss_mul_fac, ss_add_fac_out;
	ss_add_fac << fname_raw << ".add_fac";
	ss_add_fac_out << fname_raw << ".add_fac_out";	// for writing output
	ss_mul_fac << fname_raw << ".mul_fac";

	string infile_add_fac = ss_add_fac.str();
	string infile_mul_fac = ss_mul_fac.str();
	string outfile_add_fac = ss_add_fac_out.str();

	cout << "input and output add_fac and mul_fac files:\n" <<
	infile_lm << "\n" << infile_add_fac << "\n" << infile_mul_fac << "\n" << outfile_add_fac << endl;
	
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

	FILE *pInputFile_mul_fac =  fopen(infile_mul_fac.c_str(), "rb");
	if (pInputFile_mul_fac == NULL) { // check return value of file pointer
		cerr << infile_mul_fac << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pOutputFile_add_fac =  fopen(outfile_add_fac.c_str(), "wb");
	if (pOutputFile_add_fac == NULL) { // check return value of file pointer
		cerr << outfile_add_fac << " cannot be opened." << endl;
		exit(1);
	}	
	exit(1);

	cout << "-> getting normalization file from study directory: " << endl ;
	string ncfile_fullpath = argv[4];
	cout << ncfile_fullpath << endl;

	ifstream nc_read;
	nc_read.open(ncfile_fullpath.c_str(), ios::in | ios::binary);
	vector<float> nc_crys(num_crystals_all); 
	vector<float> nc_plane(num_plane_efficiencies);

	bool use_nc_default = false;
	if (!nc_read) {
		cout << "Could not open nc file, use default values (1) \n";
		use_nc_default = true;
		for (int nci = 0; nci < num_crystals_all; nci++) {
			nc_crys[nci] = 1.0;
		}
		for (int nci = 0; nci < num_plane_efficiencies; ++nci) {
			nc_plane[nci] = 1.0;
		}
	}

	if (!use_nc_default) {
		// get plane efficiencies
		nc_read.seekg( 2079453 * 4, nc_read.beg);
		for (int nc_c = 0;  nc_c < num_plane_efficiencies; ++nc_c) {
			nc_read.read(reinterpret_cast<char*>(&nc_plane[nc_c]), sizeof(float));
		}
		//get crystal efficiencies
		nc_read.seekg(3095517 * 4, nc_read.beg);
		for (int nc_c = 0; nc_c < num_crystals_all; nc_c++) {
			nc_read.read(reinterpret_cast<char*>(&nc_crys[nc_c]), sizeof(float));
		}
	}

	exit(1);

// ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo ==== main program start

	// read in events from lm file
	cout << "-> reading from list mode file...(may take a while)" << endl;

	// get buffers 
	ids * pids_in = new ids[BUFFER_SIZE];
	float *add_fac_out = new float[BUFFER_SIZE];
	float *add_fac = new float[BUFFER_SIZE];
	float *mul_fac = new float [BUFFER_SIZE];

	int num_buffers_read = 1;
	int buffer_indx = 0;
	int set_zero_ct = 0; // for debugging only
	
	while(!feof(pInputFile_lm)){
		//cout << "----> reading buffer number " << num_buffers_read << endl;

		int read_ct = fread(pids_in, sizeof(ids), BUFFER_SIZE, pInputFile_lm);
		fread(add_fac, sizeof(float), BUFFER_SIZE, pInputFile_add_fac);
		fread(mul_fac, sizeof(float), BUFFER_SIZE, pInputFile_mul_fac);
		//cout << "Events in this buffer (read_ct): "<< read_ct << endl;

		for (int i = 0; i < read_ct; ++i)
		{	
			// get crystal IDs
			short txCrysA = pids_in[i].txIDA;
			short txCrysB = pids_in[i].txIDB;
			short axCrysA = pids_in[i].axIDA - (pids_in[i].axIDA / num_ax_crys_per_unit_w_gap);
			short axCrysB = pids_in[i].axIDB - (pids_in[i].axIDB / num_ax_crys_per_unit_w_gap);

			// convert to block pids_in
			short txBiA = txCrysA / num_tx_crys_per_block;
			short axBiA = axCrysA / num_ax_crys_per_block;
			short txBiB = txCrysB / num_tx_crys_per_block;
			short axBiB = axCrysB / num_ax_crys_per_block;

			// get transaxial sinogram index
		//	int idx_tx_blk = blk_idx[txBiA][txBiB]; 
		//	int ind_blk_sino = idx_tx_blk + num_bins_sino_block * axBiA + num_bins_sino_block * num_ax_block_ring * axBiB; 

			// calculate correction factor
			float stemp;
		//	float stemp = (float)scatter_sino[ind_blk_sino];
		//	stemp = (float)stemp/num_lor_blkpair;

			stemp = stemp / ( nc_crys[axCrysA + 672*txCrysA] * nc_crys[axCrysB + 672*txCrysB]);
			
			stemp = stemp * mul_fac[i];			// apply dead time and decay correction to scatters (randoms already have it)
			float rtemp = add_fac[i];			// read in randoms from original add_fac file
			rtemp = rtemp + stemp;				// add scatters and randoms
			
			add_fac_out[i]=rtemp;
			buffer_indx++;

			if(buffer_indx == BUFFER_SIZE){
				fwrite(add_fac_out, sizeof(float), BUFFER_SIZE, pOutputFile_add_fac);
				buffer_indx = 0; // reset index
			}
		}

		num_buffers_read++;
	}
	fwrite(add_fac_out, sizeof(float), buffer_indx, pOutputFile_add_fac);
	cout << "set zero counter: "<< set_zero_ct << " (for debugging only)" << endl;

	// re-name original file and save new file under original name
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


	cout << "-> closing all files" << endl;
	delete[] pids_in;
	delete[] add_fac_out;
	fclose(pInputFile_lm);
	fclose(pInputFile_mul_fac);
	fclose(pInputFile_add_fac);
	fclose(pOutputFile_add_fac);

	cout << "done" << endl;

}
