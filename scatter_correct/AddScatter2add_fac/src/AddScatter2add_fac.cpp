#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <string>

#define BUFFER_SIZE 65536 // = elements per buffer


using namespace std;

int num_bins_sino = 549 * 420;
int num_bins_sino_block = 91 * 60;
int num_ax_block_ring = 112; // 14*8
int num_tx_block_ring = 120; // 5*24
int num_ax_crys_per_unit_w_gap = 85;
int num_ax_crys_per_block = 6;
int num_tx_crys_per_block = 7;
int num_lor_blkpair = 903;
//int num_lor_blkpair = 42*42;
int num_crystals_all_wo_gap = 564480;	//672*840
int num_plane_efficiencies_wo_gap = 451584; // 672*672
int num_plane_efficiencies = 461041; //679*679
int num_crystals_all = 570360; // 679*840

struct Lut {
	short nv,nu;
}; // block lut

struct ids {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

int main(int argc, char **argv) {

	if (argc != 8 ) {
		cout << "usage: " << argv[0] << " [lm_folder_name]  [scaled_sino]  [user_name]  [index_blockpairs_transaxial_2x91x60_int16]  [crys_eff]  [plane_eff]  [frame_number]" << endl;
		cerr << "not enough input parameters" << endl;
		exit(1);
	}

	cout << "Starting AddScatter2add_fac executable..." << endl;
	
	string infile_fullpath = argv[1];
	cout << "-> opening list-mode files from folder " << infile_fullpath << endl;

	int frame_number = atoi(argv[7]);
	stringstream ss_fn;
	ss_fn << "/lm_reorder_f" << frame_number << "_prompts";
	string lm_file_name_raw = ss_fn.str();

	stringstream ss_lm_file, ss_add_fac, ss_mul_fac, ss_add_fac_out, ss_attn_fac, ss_scat_fac;
	ss_lm_file << infile_fullpath << lm_file_name_raw << ".lm";
	ss_add_fac << infile_fullpath << lm_file_name_raw << ".add_fac";
	ss_add_fac_out << infile_fullpath << lm_file_name_raw << ".add_fac_out";	// for writing output
	ss_mul_fac << infile_fullpath << lm_file_name_raw << ".mul_fac.original";	// use original, since the other one contains only ones
	ss_attn_fac << infile_fullpath << lm_file_name_raw << ".attn_fac";
	ss_scat_fac << infile_fullpath << lm_file_name_raw << ".scat_fac";
	
	string infile_lm = ss_lm_file.str();
	string infile_add_fac = ss_add_fac.str();
	string infile_mul_fac = ss_mul_fac.str();
	string infile_attn_fac = ss_attn_fac.str();
	string outfile_add_fac = ss_add_fac_out.str();
	string outfile_scat_fac = ss_scat_fac.str();

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

	FILE *pInputFile_attn_fac =  fopen(infile_attn_fac.c_str(), "rb");
	if (pInputFile_attn_fac == NULL) { // check return value of file pointer
		cerr << infile_attn_fac << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pOutputFile_add_fac =  fopen(outfile_add_fac.c_str(), "wb");
	if (pOutputFile_add_fac == NULL) { // check return value of file pointer
		cerr << outfile_add_fac << " cannot be opened." << endl;
		exit(1);
	}	

	FILE *pOutputFile_scat_fac =  fopen(outfile_scat_fac.c_str(), "wb");
	if (pOutputFile_scat_fac == NULL) { // check return value of file pointer
		cerr << outfile_scat_fac << " cannot be opened." << endl;
		exit(1);
	}	

	cout << "-> opening scatter sino " << endl;
	string scatter_sino_path = argv[2];
	ifstream scatter_sino_read; 
	scatter_sino_read.open(scatter_sino_path.c_str(),  ios::in | ios::binary); 
	if (!scatter_sino_read) {
		cout << "could not open scatter sino" <<  endl;
		exit(1); 
	}
	vector<double> scatter_sino(num_bins_sino_block * num_ax_block_ring * num_ax_block_ring);
	for (int sci=0;  sci<(num_bins_sino_block * num_ax_block_ring * num_ax_block_ring); sci++) {
		scatter_sino_read.read(reinterpret_cast<char*>(&scatter_sino[sci]), sizeof(double));
	}
	scatter_sino_read.close();


	cout << "-> opening tof distribution " << endl;
	string user_name = argv[3];
	string tof_wt_path = "/home/";
	tof_wt_path.append(user_name);
	tof_wt_path.append("/code/explorer-master/read_lm/lut/tof_wt");
	ifstream tof_wt_read; 
	tof_wt_read.open(tof_wt_path.c_str(), ios::in | ios::binary);
	if (!tof_wt_read) {
		cerr <<  "could not open tof_wt file: " << tof_wt_path << endl; 
		exit(1);
	}
	vector<float> tof_wt(129);  
	for (int ti=0; ti<129; ti++) {
		tof_wt_read.read(reinterpret_cast<char*>(&tof_wt[ti]), sizeof(float)); 
	}
	tof_wt_read.close();  
	vector<double> tof_spectrum(129); 

	cout << "-> opening index_blockpairs_transaxial_2x91x60_int16 lut " << endl;
	string fname_lut = string(argv[4]);
	FILE* pFile_lut = fopen(fname_lut.c_str(), "rb");
	if (pFile_lut == NULL) {
		cout << fname_lut << " cannot be opened." << endl;
		exit(1);
	}
	Lut* plut = new Lut[num_bins_sino_block]; // block sino lut;  = 91 * 60;
	fread(plut, sizeof(Lut), num_bins_sino_block, pFile_lut);

	// set up block sino reverse lut
	vector< vector<int> > blk_idx(num_tx_block_ring, vector<int>(num_tx_block_ring, -1)); // initialize values to -1
	for (int i = 0; i < num_bins_sino_block; i++) {
		blk_idx[plut[i].nv][plut[i].nu] = i;
		blk_idx[plut[i].nu][plut[i].nv] = i;
	}

	//here read in crys eff and plane eff instead of nc file.  
	cout <<"-> getting crystal efficiencies from file: " << endl;
	string crys_fullpath = argv[5];
	cout << "\t" << crys_fullpath << endl;
	ifstream crys_read;
	crys_read.open(crys_fullpath.c_str(), ios::in | ios::binary);

	cout <<"-> getting plane efficiencies from file: " << endl;
	string plane_fullpath = argv[6];
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

// define coincidence time windows
	vector<float> coinc_window(8); // ns
	coinc_window[0] = 4500.0;
	coinc_window[1] = 4800.0;
	coinc_window[2] = 5400.0;
	coinc_window[3] = 6100.0;
	coinc_window[4] = 6900.0;
	coinc_window[5] = 6900.0;
	coinc_window[6] = 6900.0;
	coinc_window[7] = 6900.0;
	coinc_window[8] = 6900.0;

// ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo ==== main program start

	// read in events from lm file
	cout << "-> reading from list mode file...(may take a while)" << endl;

	// get buffers 
	ids * pids_in = new ids[BUFFER_SIZE];
	float *add_fac_out = new float[BUFFER_SIZE];
	float *add_fac = new float[BUFFER_SIZE];
	float *mul_fac = new float [BUFFER_SIZE];
	float *attn_fac = new float [BUFFER_SIZE];
	float *scat_fac = new float [BUFFER_SIZE];

	int num_buffers_read = 1;
	int buffer_indx = 0;

	int zero_ct=0;

	while(!feof(pInputFile_lm)){
		//cout << "----> reading buffer number " << num_buffers_read << endl;

		int read_ct = fread(pids_in, sizeof(ids), BUFFER_SIZE, pInputFile_lm);
		fread(add_fac, sizeof(float), BUFFER_SIZE, pInputFile_add_fac);
		fread(mul_fac, sizeof(float), BUFFER_SIZE, pInputFile_mul_fac);
		fread(attn_fac, sizeof(float), BUFFER_SIZE, pInputFile_attn_fac);
		//cout << "Events in this buffer (read_ct): "<< read_ct << endl;

		for (int i = 0; i < read_ct; ++i)
		{	
			// get crystal IDs
			short txCrysA = pids_in[i].txIDA;
			short txCrysB = pids_in[i].txIDB;
			short axCrysA = pids_in[i].axIDA; - (pids_in[i].axIDA / num_ax_crys_per_unit_w_gap);
			short axCrysB = pids_in[i].axIDB - (pids_in[i].axIDB / num_ax_crys_per_unit_w_gap);
			short axCrysA_w_gap = pids_in[i].axIDA;
			short axCrysB_w_gap = pids_in[i].axIDB;

			// unit
			int unitA = (int) axCrysA_w_gap / 85;
			int unitB = (int) axCrysB_w_gap / 85;
			float t_window = coinc_window[abs(unitA - unitB)];

			// convert to block pids_in
			short txBiA = txCrysA / num_tx_crys_per_block;
			short axBiA = axCrysA / num_ax_crys_per_block;
			short txBiB = txCrysB / num_tx_crys_per_block;
			short axBiB = axCrysB / num_ax_crys_per_block;

			// get transaxial sinogram index
			int idx_tx_blk = blk_idx[txBiA][txBiB]; 
			int ind_blk_sino = idx_tx_blk + num_bins_sino_block * axBiA + num_bins_sino_block * num_ax_block_ring * axBiB; 

			// calculate correction factor for sinogram block
			float stemp = (float)scatter_sino[ind_blk_sino];
/*
          if (axBiA == axBiB) {
				stemp = stemp / 2.0;
			}   else {
				stemp = stemp / 1.0;  
			}
*/
			stemp = (float)stemp/num_lor_blkpair;	// avg number of scatters per LOR

/*
			if(abs(pids_in[i].tof) < 64){
				stemp = stemp * tof_wt[pids_in[i].tof+64];
				//stemp = stemp / 128;
			}else{
				stemp = 0.0;
			}
*/
			stemp = stemp * (39.0625 / t_window);	// Average num of scatters per tof bin
			float pure_stemp = stemp;				// save the pure stemp, which does not contain added randoms, save further down in "scat_fac[i]"

			// add randoms (uncorrected for attenuation, etc!)
			float rtemp = add_fac[i];				// read in randoms from ORIGINAL add_fac file NOT containing normalization, attenuation and dead time!
			if (add_fac[i] == 0) zero_ct++;
	//		stemp = stemp + rtemp;					// add scatters and randoms 
			stemp = rtemp;							// remove this later FIXME !!

			//divide by mul_fac (apply dead time and decay correction)
			if(mul_fac[i] != 0){
				stemp /= mul_fac[i];
				pure_stemp /= mul_fac[i];
			}else{
				stemp = 0;
				pure_stemp = 0;
			}
			
			//divide stemp by crystal normalization i.e. multiply by UIH-defined correction factors
			pure_stemp *= (nc_crys[axCrysA_w_gap + 679*txCrysA] * nc_crys[axCrysB_w_gap + 679*txCrysB]);
			stemp *= (nc_crys[axCrysA_w_gap + 679*txCrysA] * nc_crys[axCrysB_w_gap + 679*txCrysB]);

			//divide stemp by plane normalization i.e. multiply by UIH-defined correction factors
			pure_stemp *= (nc_plane[axCrysA_w_gap + 679*axCrysB_w_gap]);
			stemp *= (nc_plane[axCrysA_w_gap + 679*axCrysB_w_gap]);

			//divide stemp by attn_fac
			if (attn_fac[i] !=0){ // catch dividing by zero. Should theoretically never happen as that would correspond to infnite attenuation
				stemp /= attn_fac[i];
				pure_stemp /= attn_fac[i];
			}else{
				stemp = 0;
				pure_stemp = 0;
			}

			scat_fac[i] = pure_stemp;	// write pure scatter correction factor to file
	//		add_fac_out[i] = stemp;		// save scatters and randoms to the output add_fac variable
			add_fac_out[i] = stemp + pure_stemp;	// revert this later FIXME !!
			

			buffer_indx++;

			if(buffer_indx == BUFFER_SIZE){
				fwrite(add_fac_out, sizeof(float), BUFFER_SIZE, pOutputFile_add_fac);
				fwrite(scat_fac , sizeof(float), BUFFER_SIZE, pOutputFile_scat_fac);
				buffer_indx = 0; // reset index

			}

		}

		num_buffers_read++;
	}

	cout << zero_ct << endl; 

	fwrite(add_fac_out, sizeof(float), buffer_indx, pOutputFile_add_fac);
	fwrite(scat_fac , sizeof(float), buffer_indx, pOutputFile_scat_fac);

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
	delete[] add_fac;
	delete[] attn_fac;
	delete[] scat_fac;
	delete[] mul_fac;
	fclose(pInputFile_lm);
	fclose(pInputFile_mul_fac);
	fclose(pInputFile_add_fac);
	fclose(pOutputFile_add_fac);
	fclose(pOutputFile_scat_fac);
	fclose(pInputFile_attn_fac);

	cout << "done" << endl;


}
