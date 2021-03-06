/*! \file subsample.hh
 \brief Implementation of Class Subsample.cpp 
 Subsampling of the EXPLORER data to mimic shorter scanners or damages or axial gaps for TO scanners
 This class implements a quasi-LUT to check if a coincidence between a pair of crystals A and B is to be kept for reconstruction or not. 
 Decision is made based on the crystal efficiency array created by the matlab script respiratory_gui.m
*/

#include "Subsample.h"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Subsample::Subsample(std::string s, int t) /*!< Takes the full path to the .raw file as input. */
	: input_raw_fullpath(s)
	, firstTimeStamp(t)
{
	createLogFile = true;
	Initialize();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Subsample::~Subsample(){		// Destructor 

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Subsample::Initialize(){

	/// get file path to config file
	std::string str_find = "PET";
	std::string input_config_TEMP;
	size_t found = input_raw_fullpath.find(str_find);
	input_config_TEMP = input_raw_fullpath.substr(0, found);
	stringstream ss ;
	ss << input_config_TEMP << "UCD/multi_bed_recon.config" ;
	input_config = ss.str();

	if (createLogFile)
	{
		/// get system time to create timestamp for log file name
		std::time_t t = std::time(0);   	// get time now
	    std::tm* now = std::localtime(&t);

		/// set logfile out path
		stringstream ss_log;
		ss_log << input_config_TEMP << "UCD/subsample_" << now->tm_year+1900 << "-" << now->tm_mon + 1 
			<< "-" << now->tm_mday << "_" << now->tm_hour << ":" << now->tm_min << ".log" ;
		output_LOG = ss_log.str();
		o_LOG.open(output_LOG.c_str(), std::ios_base::app);

		o_LOG << input_raw_fullpath << endl;
	}

	cout << "====================================\nreading multi bed config file \n" << input_config << endl;

	ifstream config;
	config.open(input_config.c_str());
	if (!config) {
		cout << "could not open multi bed recon config file" << endl;
		exit(1);
	}
	std::string str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp(str_temp);
	ss_temp >> num_cycles;

	str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp2(str_temp);
	ss_temp2 >> start_ring;

	str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp3(str_temp);
	ss_temp3 >> rings_per_bed;

	str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp4(str_temp);
	ss_temp4 >> num_beds;	

	str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp5(str_temp);
	ss_temp5 >> bed_overlap;	


	str_temp =  "";
	getline(config, str_temp);
	stringstream ss_temp6(str_temp);
	ss_temp6 >>	time_per_bed;

	cout << "number of scans per bed position:\t" << num_cycles 
	<< "\nstart ring\t\t\t\t" << start_ring 
	<< "\nnumber of rings per bed position\t" << rings_per_bed
	<< "\nnumber of bed positions\t\t\t" << num_beds
	<< "\noverlap of bed positions\t\t" << bed_overlap 
	<< "\ntime per bed\t\t\t\t" << time_per_bed
	<< "\n===================================="<< endl;


	for (int cycle = 0; cycle < num_cycles; ++cycle)
	{
		for (int pos = 0; pos < num_beds; ++pos)
		{
		bed_time_start.push_back(cycle*time_per_bed*num_beds + pos*time_per_bed);
		bed_time_end.push_back(cycle*time_per_bed*num_beds + (pos+1)*time_per_bed);
		int this_ring_start = start_ring+pos*(rings_per_bed-bed_overlap);
		bed_ring_start.push_back(this_ring_start);
		bed_ring_end.push_back(this_ring_start+rings_per_bed-1);
		}
	}

	for (int i = 0; i < bed_time_start.size(); ++i)
	{
		cout << "start time " << i << ": " << bed_time_start.at(i) 
		<< ". end time: " << bed_time_end.at(i)
		<< ". start ring: " << bed_ring_start.at(i)
		<< ". end ring: " << bed_ring_end.at(i)
		<< endl;

	}

	cout << "first time stamp set: \t" << firstTimeStamp << "sec = " << (int)firstTimeStamp/3600 << "h:" << ((int)firstTimeStamp%3600)/60 << "min:" << ((int)firstTimeStamp%3600)%60 << "s" << endl;


// fill log file if desired
	if (createLogFile)
	{
		o_LOG << "number of scans per bed position:\t" << num_cycles 
		<< "\nstart ring\t\t\t\t" << start_ring 
		<< "\nnumber of rings per bed position\t" << rings_per_bed
		<< "\nnumber of bed positions\t\t\t" << num_beds
		<< "\noverlap of bed positions\t\t" << bed_overlap 
		<< "\ntime per bed\t\t\t\t" << time_per_bed
		<< "\n===================================="<< endl;

		for (int i = 0; i < bed_time_start.size(); ++i)
		{
			o_LOG << "start time " << i << ": " << bed_time_start.at(i) 
			<< ". end time: " << bed_time_end.at(i)
			<< ". start ring: " << bed_ring_start.at(i)
			<< ". end ring: " << bed_ring_end.at(i)
			<< endl;
		}

		o_LOG << "first time stamp set: \t" << firstTimeStamp << "sec = " << (int)firstTimeStamp/3600 << "h:" << ((int)firstTimeStamp%3600)/60 << "min:" << ((int)firstTimeStamp%3600)%60 << "s" << endl;
		o_LOG.close();
	}	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
	bool Subsample::KeepEvent(int axA, int axB, int transA, int transB, int timeStamp) // Takes the axial and transaxial coordinates of the two coincident crystals as input.
	{
		// adjust input indexing (0-indexing) to indexing used in this subsample code (matlab indexing starting at 1)
		axA++;
		axB++;
		transA++;
		transB++;
		
		int time_elapsed = timeStamp - firstTimeStamp;	// time elapsed since first time stamp. This essentially determines the bed position
		int pos = 0;
		while(pos < num_beds*num_cycles){	// go through all possible bed possitions and find the one with matching time frame
			if(time_elapsed>=bed_time_end.at(pos)) {
				pos++;
			}else{
				break;
			}
		}

		// if pos exceeds total number of bed positions, then time stamp of current event is larger than the maximum of the last bed position -> reject!
		if(pos >= num_beds*num_cycles) return false;	

		if(axA < bed_ring_start.at(pos) || axA > bed_ring_end.at(pos) || axB < bed_ring_start.at(pos) || axB > bed_ring_end.at(pos)) {
	//		cout << time_elapsed << "\t" << pos << "\t" << bed_time_start.at(pos) << "\t" << bed_time_end.at(pos) << endl;
	//		cout << axA << "\t" << axB << "\t" << bed_ring_start.at(pos) << "\t" << bed_ring_end.at(pos) << " -> reject!" << endl;
			return false;
		}else{
	//		cout << time_elapsed << "\t" << pos << "\t" << bed_time_start.at(pos) << "\t" << bed_time_end.at(pos) << endl;
	//		cout << axA << "\t" << axB << "\t" << bed_ring_start.at(pos) << "\t" << bed_ring_end.at(pos) << " -> keep!" << endl;			
			return true;
		}


	}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
