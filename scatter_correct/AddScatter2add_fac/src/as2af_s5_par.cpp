#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/time.h>
#include <sys/stat.h>
//#include <math.h>
#include <thread>
#include <future>
#include <unistd.h>
#include <chrono>
 
#define BUFFER_SIZE 65536 // = elements per buffer
#define N_TOF 27 //UIH defined [15,17,19,23,25,27], we use 27 to maintain compatibility with old format
#define N_TIME_BIN_PER_TOF 7

using namespace std;

const auto processor_count = std::thread::hardware_concurrency();

int num_bins_sino = 549 * 420;
int num_bins_sino_block = 91 * 60;
int num_ax_block_ring = 112; // 14*8
int num_tx_block_ring = 120; // 5*24
int num_ax_crys_per_unit_w_gap = 85;
int num_ax_crys_per_block = 6;
int num_tx_crys_per_block = 7;
//int num_lor_blkpair = 903;
int num_lor_blkpair = 42*42;
int num_crystals_all_wo_gap = 564480;	//672*840
int num_plane_efficiencies_wo_gap = 451584; // 672*672
int num_plane_efficiencies = 461041; //679*679
int num_crystals_all = 570360; // 679*840

vector<float> nc_crys(num_crystals_all); 
vector<float> nc_plane(num_plane_efficiencies);

vector<float> scatter_sino(num_bins_sino_block * num_ax_block_ring * num_ax_block_ring * N_TOF);

vector< vector<int> > blk_idx(num_tx_block_ring, vector<int>(num_tx_block_ring, -1)); // initialize values to -1

string infile_fullpath;
string infile_lm;
string infile_add_fac;
string infile_add_fac_orig;
string infile_mul_fac;
string infile_attn_fac;
string outfile_add_fac;
string outfile_scat_fac;

int frame_number;

struct Lut {
	short nv,nu;
}; // block lut

struct ids {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

struct filePos {
	long long fStartPos, fEndPos;
};

int AddScatter2add_fac(long long startPos, long long endPos, short threadID);

// ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo ==== main function

int main(int argc, char **argv) {

	

	if (argc != 7 ) {
		cout << "usage: " << argv[0] << " [lm_folder_name]  [scaled_sino]  [index_blockpairs_transaxial_2x91x60_int16]  [crys_eff]  [plane_eff]  [frame_number]" << endl;
		cerr << "not enough input parameters" << endl;
		exit(1);
	}

	cout << "Starting AddScatter2add_fac executable..." << endl;
	cout << "Available threads: " << processor_count << endl;


// Read In Listmode Files and save file names

	infile_fullpath = argv[1];
	cout << "-> opening list-mode files from folder " << infile_fullpath << endl;

	frame_number = atoi(argv[6]);
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
	
	infile_lm = ss_lm_file.str();
	infile_add_fac = ss_add_fac.str();
	infile_mul_fac = ss_mul_fac.str();
	infile_attn_fac = ss_attn_fac.str();
	outfile_add_fac = ss_add_fac_out.str();
	outfile_scat_fac = ss_scat_fac.str();

	// always use original add_fac file coming from the read_lm executable 
// check if original file exists
	cout << "-> checking if original add_fac file exists" << endl;
	infile_add_fac_orig = infile_add_fac;
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

	fclose(pInputFile_lm);
	fclose(pInputFile_mul_fac);
	fclose(pInputFile_add_fac);
	fclose(pOutputFile_add_fac);
	fclose(pOutputFile_scat_fac);

	cout << "-> Creating temp folder for parallel threads" << endl;
	string temp_out = infile_fullpath;
	temp_out.append("/temp");
	int status = mkdir(temp_out.c_str(), 0775);
	if (status != 0){
		cout << "----> temp folder was not created at this location\n" << temp_out << "\n\t(probably exists or directory is invalid.)" << endl;	// make this nicer! FIXME
	}
// Read in Look Up Tables, Scatter Sino and Detector Normalization Files

	cout << "-> opening scatter sino " << endl;
	string scatter_sino_path = argv[2];
	ifstream scatter_sino_read; 
	scatter_sino_read.open(scatter_sino_path.c_str(),  ios::in | ios::binary); 
	if (!scatter_sino_read) {
		cout << "could not open scatter sino" <<  endl;
		exit(1); 
	}
	
	for (int sci=0;  sci<(num_bins_sino_block * num_ax_block_ring * num_ax_block_ring * N_TOF); sci++) {
		scatter_sino_read.read(reinterpret_cast<char*>(&scatter_sino[sci]), sizeof(float));
	}
	scatter_sino_read.close();

	cout << "-> opening index_blockpairs_transaxial_2x91x60_int16 lut " << endl;
	string fname_lut = string(argv[3]);
	FILE* pFile_lut = fopen(fname_lut.c_str(), "rb");
	if (pFile_lut == NULL) {
		cout << fname_lut << " cannot be opened." << endl;
		exit(1);
	}
	Lut* plut = new Lut[num_bins_sino_block]; // block sino lut; num_bins_sino_block = 91 * 60;
	fread(plut, sizeof(Lut), num_bins_sino_block, pFile_lut);

// set up block sino reverse lut
	for (int i = 0; i < num_bins_sino_block; i++) {
		blk_idx[plut[i].nv][plut[i].nu] = i;
	//	blk_idx[plut[i].nu][plut[i].nv] = i;
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




// get number of events per thread
	ifstream infile;
	infile.open(infile_lm.c_str(), ios::in | ios::binary); // open list mode file
	infile.seekg(0, infile.end);
	long long file_size = infile.tellg(); //get size (events) of list mode file
	long long num_events = file_size / sizeof(ids);
	cout << "-> " << num_events << " total events" << endl;
	infile.close();

	filePos* fPos = new filePos[processor_count];
	long long events_per_thread = ceil((float)num_events/processor_count);
	cout << "-> Number of events per thread: " << events_per_thread << endl;
	for (int i = 0; i < processor_count-1; ++i)
	{
		fPos[i].fStartPos = i*events_per_thread;
		fPos[i].fEndPos = fPos[i].fStartPos + events_per_thread-1;
	}
	fPos[processor_count-1].fStartPos = fPos[processor_count-2].fEndPos+1;
	fPos[processor_count-1].fEndPos = num_events-1;



	future<int> fut[processor_count];
	thread t[processor_count];

	for (short i = 0; i < processor_count; ++i)
	{
		long long s = fPos[i].fStartPos;
		long long e = fPos[i].fEndPos;
		packaged_task<int(long long, long long, short)> task(&AddScatter2add_fac);// wrap the function
	    if ( task.valid() ) {
	    //	cout << "Task " << i << " is valid" << endl;
	    }else{
	    	cout << "Task " << i << " NOT valid!" << endl;
	    }
	    fut[i] = task.get_future();  						// get a future
	    t[i] = thread(move(task),s,e,i);// launch on a thread
	    usleep(10000);	// sleep for XXX microseconds
	}


// wait for the processes to finish
	sleep(1);
	bool process_finished[processor_count] = {false};
	short num_processes_done = 0;
	int interval = 4;
	int while_loop_counter = 1;
	bool not_done = true;
	while (not_done)
	{
		for (int i = 0; i < processor_count; ++i)
		{
			if(process_finished[i]) continue;
			future_status status = fut[i].wait_for(std::chrono::seconds(0));
	        if (status == future_status::deferred) {
	        //    cout << "status process " << i << ": deferred\n";
	        } else if (status == future_status::timeout) {
	        //    cout << "status process " << i << ": timeout\n";
	        } else if (status == future_status::ready) {
	        //    cout << "status process " << i << ": ready!\n";
	            num_processes_done++;
	            process_finished[i] = true;
	        }
		}// end of for
		//check status flag from all processes:
		cout << "----> # processes done: " << num_processes_done << " (time since start: " << interval*while_loop_counter << "s)" << endl;
		if(num_processes_done == (short)processor_count){
			not_done = false;
		}else{
			sleep(interval);
			while_loop_counter++;
		}
	}// end of while not_done

	for (int i = 0; i < processor_count; ++i)
	{
		t[i].join();
	}

	cout << "-> done all threads." << endl;

	usleep(10000);

// perform clean up

	// concatenate files
	cout << "-> concatenate files..." << endl;
	stringstream ss_concat_add_fac, ss_concat_scat_fac;
	ss_concat_add_fac << "cat";
	ss_concat_scat_fac << "cat";

	for (int i = 0; i < processor_count; ++i)
	{
		ss_concat_add_fac << " " << infile_fullpath << "/temp/lm_reorder_f" << frame_number << "_prompts.add_fac_out." << i ;
		ss_concat_scat_fac << " " << infile_fullpath << "/temp/lm_reorder_f" << frame_number << "_prompts.scat_fac." << i ;
	}
	ss_concat_add_fac << " > " << infile_fullpath << "/lm_reorder_f" << frame_number << "_prompts.add_fac_out";
	ss_concat_scat_fac << " > " << infile_fullpath << "/lm_reorder_f" << frame_number << "_prompts.scat_fac";
	string cmd_concat_1 = ss_concat_add_fac.str();
	string cmd_concat_2 = ss_concat_scat_fac.str();
	system(cmd_concat_1.c_str());
	system(cmd_concat_2.c_str());

	// delete temp folder
	cout << "-> deleting temp files..." ;
	stringstream ss_delete;
	ss_delete << "rm -r " << infile_fullpath << "/temp" << endl;
	string delete_temp = ss_delete.str();
	system(delete_temp.c_str());


	// re-name original file and save new file under original name
	cout << "\n-> deleting existing *.add_fac and renaming new file to *.add_fac" << endl;

	if ( remove(infile_add_fac.c_str()) != 0 ){
		cerr << "!! could not remove previous file " << infile_add_fac << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully deleted previous file\n\t" << infile_add_fac << endl;
	}

	string outfile_add_fac_new = outfile_add_fac;
	outfile_add_fac_new.erase(outfile_add_fac.length() -4, 4);
	if (rename(outfile_add_fac.c_str() , outfile_add_fac_new.c_str()) != 0){
		cerr << "!! renaming file\n" << outfile_add_fac << "\nto\n" << outfile_add_fac_new << "\nfailed." << endl;
		exit(1);
	}else{
		cerr << "----> sucessfully renamed file\n\t" << outfile_add_fac << "\nto\n\t" << outfile_add_fac_new << endl;
	}	 


	cout << "=====\nDONE." << endl;
}


// ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo ==== main function

int AddScatter2add_fac(long long startPos, long long endPos, short threadID){

	cout << "----> starting thread " << threadID << ", from evtNum " << startPos << " to " << endPos << endl;

	FILE *pInputFile_lm =  fopen(infile_lm.c_str(), "rb");
	FILE *pInputFile_add_fac =  fopen(infile_add_fac_orig.c_str(), "rb");
	FILE *pInputFile_mul_fac =  fopen(infile_mul_fac.c_str(), "rb");
	FILE *pInputFile_attn_fac =  fopen(infile_attn_fac.c_str(), "rb");
// outfile add_fac
	string outfile_add_fac_thread = infile_fullpath;
	outfile_add_fac_thread.append("/temp/lm_reorder_f");
	outfile_add_fac_thread.append(to_string(frame_number));
	outfile_add_fac_thread.append("_prompts.add_fac_out.");
	outfile_add_fac_thread.append(to_string(threadID));
	FILE *pOutputFile_add_fac =  fopen(outfile_add_fac_thread.c_str(), "wb");
// outfile scat_fac
	string outfile_scat_fac_thread = infile_fullpath;
	outfile_scat_fac_thread.append("/temp/lm_reorder_f");
	outfile_scat_fac_thread.append(to_string(frame_number));
	outfile_scat_fac_thread.append("_prompts.scat_fac.");
	outfile_scat_fac_thread.append(to_string(threadID));
	FILE *pOutputFile_scat_fac =  fopen(outfile_scat_fac_thread.c_str(), "wb");

// get buffers 
	ids * pids_in = new ids[BUFFER_SIZE];
	float *add_fac = new float[BUFFER_SIZE];
	float *mul_fac = new float [BUFFER_SIZE];
	float *attn_fac = new float [BUFFER_SIZE];
	float *scat_fac = new float [BUFFER_SIZE];
	float *add_fac_out = new float[BUFFER_SIZE];

// get start position
	fseek(pInputFile_lm, startPos*sizeof(ids), SEEK_SET); // SEEK_SET is the beginning of the file
	fseek(pInputFile_mul_fac, startPos*sizeof(float), SEEK_SET); // SEEK_SET is the beginning of the file
	fseek(pInputFile_add_fac, startPos*sizeof(float), SEEK_SET); // SEEK_SET is the beginning of the file
	fseek(pInputFile_attn_fac, startPos*sizeof(float), SEEK_SET); // SEEK_SET is the beginning of the file

	long long total_events_to_read = endPos - startPos +1;
	long long remaining_events = total_events_to_read;
	int current_BUFFER_SIZE;
	int buffer_indx=0;

	int lg26_ct=0;
	int sm0_ct=0;

// loop over events in file section
	bool still_reading = true;
	while(still_reading){
		if (remaining_events >= BUFFER_SIZE){
			current_BUFFER_SIZE = BUFFER_SIZE;
		}else {
			current_BUFFER_SIZE = remaining_events;
		}

		int read_ct = fread(pids_in, sizeof(ids), current_BUFFER_SIZE, pInputFile_lm);
		fread(add_fac, sizeof(float), current_BUFFER_SIZE, pInputFile_add_fac);
		fread(mul_fac, sizeof(float), current_BUFFER_SIZE, pInputFile_mul_fac);
		fread(attn_fac, sizeof(float), current_BUFFER_SIZE, pInputFile_attn_fac);

		for (int i = 0; i < read_ct; ++i)
		{	
			// get crystal IDs
			short txCrysA = pids_in[i].txIDA;
			short txCrysB = pids_in[i].txIDB;
			short axCrysA = pids_in[i].axIDA - (pids_in[i].axIDA / num_ax_crys_per_unit_w_gap);
			short axCrysB = pids_in[i].axIDB - (pids_in[i].axIDB / num_ax_crys_per_unit_w_gap);
			short axCrysA_w_gap = pids_in[i].axIDA;
			short axCrysB_w_gap = pids_in[i].axIDB;
			short TOF_AB = pids_in[i].tof;
			if(TOF_AB >= 0){
				TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
			}else{
				TOF_AB -=7;
				TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
			}
			if(TOF_AB+13 > 26) {TOF_AB = 13; lg26_ct++;}	// catch out of bound values
			if(TOF_AB+13 < 0) {TOF_AB = -13; sm0_ct++;}

			// unit
			int unitA = (int) axCrysA_w_gap / 85;
			int unitB = (int) axCrysB_w_gap / 85;

			// convert to block pids_in
			short txBiA = txCrysA / num_tx_crys_per_block;
			short axBiA = axCrysA / num_ax_crys_per_block;
			short txBiB = txCrysB / num_tx_crys_per_block;
			short axBiB = axCrysB / num_ax_crys_per_block;

			// get transaxial sinogram index
			int idx_tx_blk = blk_idx[txBiA][txBiB]; 		// transaxial sinogram index ,120×120
			int idx_tx_blk_reverse = blk_idx[txBiB][txBiA]; // transaxial sinogram index ,120×120,reverse direction
			int ind_blk_sino;

			if (idx_tx_blk !=-1){
				ind_blk_sino = idx_tx_blk 
				+ num_bins_sino_block * axBiB 
				+ num_bins_sino_block * num_ax_block_ring * axBiA
				+ num_bins_sino_block * num_ax_block_ring * num_ax_block_ring * (TOF_AB+13); 
			}else if(idx_tx_blk_reverse !=-1){
				TOF_AB=-TOF_AB;
				ind_blk_sino = idx_tx_blk 
				+ num_bins_sino_block * axBiA
				+ num_bins_sino_block * num_ax_block_ring * axBiB
				+ num_bins_sino_block * num_ax_block_ring * num_ax_block_ring * (TOF_AB+13); 
			}
			// calculate correction factor for sinogram block
			float stemp = (float)scatter_sino[ind_blk_sino] / N_TIME_BIN_PER_TOF;
			
			// divide by N_TIME_BIN_PER_TOF since that amount of tof bins is combined in the scatter sinogram

			stemp = (float)stemp/num_lor_blkpair;	// avg number of scatters per LOR

/*
          if (axBiA == axBiB) {
				stemp = stemp / 2.0;
			}   else {
				stemp = stemp / 1.0;  
			}

			if(abs(pids_in[i].tof) < 64){
				stemp = stemp * tof_wt[pids_in[i].tof+64];
				//stemp = stemp / 128;
			}else{
				stemp = 0.0;
			}
*/

//			stemp = stemp * (39.0625 / t_window);	// Average num of scatters per tof bin

			// add randoms (uncorrected for attenuation, etc!)
			float rtemp = add_fac[i];				// read in randoms from ORIGINAL add_fac file NOT containing normalization, attenuation and dead time!
			
			//divide by mul_fac (apply dead time and decay correction)
			if(mul_fac[i] != 0){
				rtemp /= mul_fac[i];
			//	stemp /= mul_fac[i];
			}else{
				rtemp = 0;
				stemp = 0;
			}
			
			//divide rtemp by crystal normalization i.e. multiply by UIH-defined correction factors
		//	stemp *= (nc_crys[axCrysA_w_gap + 679*txCrysA] * nc_crys[axCrysB_w_gap + 679*txCrysB]);
			rtemp *= (nc_crys[axCrysA_w_gap + 679*txCrysA] * nc_crys[axCrysB_w_gap + 679*txCrysB]);

			//divide rtemp by plane normalization i.e. multiply by UIH-defined correction factors
			stemp *= (nc_plane[axCrysA_w_gap + 679*axCrysB_w_gap]);
			rtemp *= (nc_plane[axCrysA_w_gap + 679*axCrysB_w_gap]);

			//divide rtemp by attn_fac
			if (attn_fac[i] !=0){ // catch dividing by zero. Should theoretically never happen as that would correspond to infnite attenuation
				rtemp /= attn_fac[i];
				stemp /= attn_fac[i];
			}else{
				rtemp = 0;
				stemp = 0;
			}

			scat_fac[i] = stemp;	// write pure scatter correction factor to file
	//		add_fac_out[i] = rtemp;		// save scatters and randoms to the output add_fac variable
			add_fac_out[i] = rtemp + stemp;	// revert this later FIXME !!
			

			buffer_indx++;

			if(buffer_indx == BUFFER_SIZE){
				fwrite(add_fac_out, sizeof(float), BUFFER_SIZE, pOutputFile_add_fac);
				fwrite(scat_fac , sizeof(float), BUFFER_SIZE, pOutputFile_scat_fac);
				buffer_indx = 0; // reset index
			}

		}// end of for

		remaining_events -= current_BUFFER_SIZE;
		if (remaining_events <= 0) still_reading = 0;
	} // end of while

	fwrite(add_fac_out, sizeof(float), buffer_indx, pOutputFile_add_fac);
	fwrite(scat_fac , sizeof(float), buffer_indx, pOutputFile_scat_fac);


	fclose(pInputFile_lm);
	fclose(pInputFile_mul_fac);
	fclose(pInputFile_add_fac);
	fclose(pInputFile_attn_fac);
	fclose(pOutputFile_add_fac);
	fclose(pOutputFile_scat_fac);

	return 0;
}
