// Inserts coincidence data from one dataset (e.g. point source) into another dataset
// Eric Berg; July 2019
// Edwin Leung; April 2020

#define BUFFER_SIZE 65536

#include <iostream>
#include <fstream>
#include <math.h> // deprecated - consider using <cmath>
#include <random>
#include <chrono>
#include <algorithm> // all_of

using namespace std;

struct ids {
	short txID1, axID1, txID2, axID2, tof;
};

class Buffer {	// buffer for list mode data
public:
	ids *pBuf;
	unsigned long long index;
	unsigned long long read_count;

	Buffer();
	virtual ~Buffer();
};

Buffer::Buffer() {
	pBuf = new ids[BUFFER_SIZE];
	index = 0;
	read_count = 0;
}

Buffer::~Buffer() {
	delete[] pBuf;
}

class BufferF {	// buffer for correction _F_actors
public:
	float *pBuf;
	unsigned long long index;
	unsigned long long read_count;

	BufferF();
	virtual ~BufferF();
};

BufferF::BufferF() {
	pBuf = new float[BUFFER_SIZE];
	index = 0;
	read_count = 0;
}

BufferF::~BufferF() {
	delete[] pBuf;
}

int main(int argc, char **argv) {

	// check arguments
	if (argc != 10) {
		cout << "Usage: " << argv[0]
				<< " [fname_out] [fname1_in] ... [fname8_in]" << endl;
		cout << "All filenames must be specified with .lm extension." << endl;
		return 1;
	}

	// variables that are set by user in GUI
	string outname; // outname_p
	outname = argv[1]; // .lm
	string outname_s = outname.substr(0, outname.find_last_of("."))
			+ ".add_fac";
	string outname_s_orig = outname.substr(0, outname.find_last_of("."))
			+ ".add_fac.original";
	string outname_m = outname.substr(0, outname.find_last_of("."))
			+ ".mul_fac";	
	string outname_m_orig = outname.substr(0, outname.find_last_of("."))
			+ ".mul_fac.original";
	string outname_a = outname.substr(0, outname.find_last_of("."))
			+ ".attn_fac";

	vector<string> inname_p(8, ""); 	// 8 files
	vector<string> inname_s(8, ""); 	// 8 files
	vector<string> inname_s_o(8, "");	// 8 files
	vector<string> inname_m(8, ""); 	// 8 files
	vector<string> inname_a(8, ""); 	// 8 files
	vector<string> inname_m_o(8, ""); 	// 8 files

	for (auto i = 0; i < 8; i++) { // 8 files
		inname_p.at(i) = argv[i + 2]; // i+2 because files 1-8 are argv[2-9] (0-index)
		inname_s.at(i) = inname_p.at(i).substr(0,
				inname_p.at(i).find_last_of(".")) + ".add_fac";
		inname_s_o.at(i) = inname_p.at(i).substr(0,
				inname_p.at(i).find_last_of(".")) + ".add_fac.original";
		inname_m.at(i) = inname_p.at(i).substr(0,
				inname_p.at(i).find_last_of(".")) + ".mul_fac";
		inname_m_o.at(i) = inname_p.at(i).substr(0,
				inname_p.at(i).find_last_of(".")) + ".mul_fac.original";
		inname_a.at(i) = inname_p.at(i).substr(0,
				inname_p.at(i).find_last_of(".")) + ".attn_fac";
	}

	// open input files
	FILE **pInFile_p = new FILE*[8]; // 8 files
	FILE **pInFile_s = new FILE*[8]; // 8 files	
	FILE **pInFile_s_o = new FILE*[8]; // 8 files
	FILE **pInFile_m = new FILE*[8]; // 8 files
	FILE **pInFile_m_o = new FILE*[8]; // 8 files
	FILE **pInFile_a = new FILE*[8]; // 8 files
	Buffer *pInBuf_p = new Buffer[8]; // 8 files
	BufferF *pInBuf_s = new BufferF[8]; // 8 files
	BufferF *pInBuf_s_o = new BufferF[8]; // 8 files
	BufferF *pInBuf_m = new BufferF[8]; // 8 files
	BufferF *pInBuf_m_o = new BufferF[8]; // 8 files
	BufferF *pInBuf_a = new BufferF[8]; // 8 files

	for (auto i = 0; i < 8; i++) { // 8 files
		pInFile_p[i] = fopen(inname_p.at(i).c_str(), "rb"); // open list-mode file
		pInFile_s[i] = fopen(inname_s.at(i).c_str(), "rb"); // open .add_fac file
		pInFile_s_o[i] = fopen(inname_s_o.at(i).c_str(), "rb"); // open .add_fac original file
		pInFile_m[i] = fopen(inname_m.at(i).c_str(), "rb"); // open .mul_fac file
		pInFile_m_o[i] = fopen(inname_m_o.at(i).c_str(), "rb"); // open .mul_fac.original file
		pInFile_a[i] = fopen(inname_a.at(i).c_str(), "rb"); // open .attn_fac file
		if (pInFile_p[i] == NULL) {
			cout << "Cannot open " << inname_p.at(i) << endl;
			exit(1);
		}
		if (pInFile_s[i] == NULL) {
			cout << "Cannot open " << inname_s.at(i) << endl;
			exit(1);
		}
		if (pInFile_s_o[i] == NULL) {
			cout << "Cannot open " << inname_s_o.at(i) << endl;
			exit(1);
		}
		if (pInFile_m[i] == NULL) {
			cout << "Cannot open " << inname_m.at(i) << endl;
			exit(1);
		}
		if (pInFile_m_o[i] == NULL) {
			cout << "Cannot open " << inname_m_o.at(i) << endl;
			exit(1);
		}
		if (pInFile_a[i] == NULL) {
			cout << "Cannot open " << inname_a.at(i) << endl;
			exit(1);
		}
	}

	// open output files
	FILE *pOutFile = fopen(outname.c_str(), "wb");
	FILE *pOutFile_s = fopen(outname_s.c_str(), "wb");
	FILE *pOutFile_s_o = fopen(outname_s_orig.c_str(), "wb");
	FILE *pOutFile_m = fopen(outname_m.c_str(), "wb");
	FILE *pOutFile_m_o = fopen(outname_m_orig.c_str(), "wb");
	FILE *pOutFile_a = fopen(outname_a.c_str(), "wb");
	Buffer outBuf;
	BufferF outBuf_s, outBuf_s_o, outBuf_m, outBuf_a, outBuf_m_o;

	long long max_events_file = 0;

	vector<long long> file_size_p(8, 0); // 8 files
	vector<long long> num_events_p(8, 0); // 8 files

	// get number of events for each file
	for (auto i = 0; i < 8; i++) { // 8 files
		ifstream ifs;
		ifs.open(inname_p.at(i).c_str(), ios::in | ios::binary); // open list-mode file
		ifs.seekg(0, ifs.end);
		file_size_p.at(i) = ifs.tellg(); //get size of list mode file
		num_events_p.at(i) = file_size_p.at(i) / 10; // 10-byte events
		if (num_events_p.at(i) > max_events_file) {
			max_events_file = num_events_p.at(i);
		}
		ifs.seekg(0, ifs.beg); // rewind to the beginning of the file
	}

	// calculate probability for each file
	vector<double> p_prob(8, 0); // 8 files
	cout << "File sizes relative to largest file:" << endl;
	for (auto i = 0; i < 8; i++) { // 8 files
		p_prob.at(i) = (double) num_events_p.at(i) / (double) max_events_file;
		if ( p_prob.at(i) == 0) {
			p_prob.at(i) = 0.01;		// empty files still need to be processed, thus, probability cannot be zero! rbayerlein
		}
		cout << "*" << i << ".lm\t" << p_prob.at(i) << endl;
	}

	// **************		Main Run Program		***************************
	// ************************************************************************

	// random number generator
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // why unsigned int not unsigned long long?
	default_random_engine generator(seed);
	uniform_real_distribution<double> distribution(0.0, 1.0);

	vector<bool> run_fread(8, true); // 8 files
	bool run = true;

	// initial read
	for (auto i = 0; i < 8; i++) { // 8 files
		pInBuf_p[i].read_count = fread(pInBuf_p[i].pBuf, sizeof(ids),
		BUFFER_SIZE, pInFile_p[i]);
		pInBuf_s[i].read_count = fread(pInBuf_s[i].pBuf, sizeof(float),
		BUFFER_SIZE, pInFile_s[i]);
		pInBuf_m[i].read_count = fread(pInBuf_m[i].pBuf, sizeof(float),
		BUFFER_SIZE, pInFile_m[i]);
		pInBuf_m_o[i].read_count = fread(pInBuf_m_o[i].pBuf, sizeof(float),
		BUFFER_SIZE, pInFile_m_o[i]);
		pInBuf_a[i].read_count = fread(pInBuf_a[i].pBuf, sizeof(float),
		BUFFER_SIZE, pInFile_a[i]);
	}

	// main loop
	while (run) {

		double r_num = distribution(generator); // generate random number

		for (auto i = 0; i < 8; i++) { // 8 files
			// start probability check
			if (r_num < p_prob.at(i) && run_fread.at(i)) {
				if (pInBuf_p[i].index < pInBuf_p[i].read_count) {
					outBuf.pBuf[outBuf.index] =
							pInBuf_p[i].pBuf[pInBuf_p[i].index];
					pInBuf_p[i].index++;
					if (pInBuf_p[i].index == BUFFER_SIZE) {
						pInBuf_p[i].read_count = fread(pInBuf_p[i].pBuf,
								sizeof(ids),
								BUFFER_SIZE, pInFile_p[i]);
						pInBuf_p[i].index = 0;
					}
					outBuf.index++;
					if (outBuf.index == BUFFER_SIZE) {
						fwrite(outBuf.pBuf, sizeof(ids), BUFFER_SIZE, pOutFile);
						outBuf.index = 0;
					}

					// add_fac files
					outBuf_s.pBuf[outBuf_s.index] =
							pInBuf_s[i].pBuf[pInBuf_s[i].index];
					pInBuf_s[i].index++;
					if (pInBuf_s[i].index == BUFFER_SIZE) {
						pInBuf_s[i].read_count = fread(pInBuf_s[i].pBuf,
								sizeof(float),
								BUFFER_SIZE, pInFile_s[i]);
						pInBuf_s[i].index = 0;
					}
					outBuf_s.index++;
					if (outBuf_s.index == BUFFER_SIZE) {
						fwrite(outBuf_s.pBuf, sizeof(float), BUFFER_SIZE,
								pOutFile_s);
						outBuf_s.index = 0;
					}

					// add_fac original files
					outBuf_s_o.pBuf[outBuf_s_o.index] =
							pInBuf_s_o[i].pBuf[pInBuf_s_o[i].index];
					pInBuf_s_o[i].index++;
					if (pInBuf_s_o[i].index == BUFFER_SIZE) {
						pInBuf_s_o[i].read_count = fread(pInBuf_s_o[i].pBuf,
								sizeof(float),
								BUFFER_SIZE, pInFile_s_o[i]);
						pInBuf_s_o[i].index = 0;
					}
					outBuf_s_o.index++;
					if (outBuf_s_o.index == BUFFER_SIZE) {
						fwrite(outBuf_s_o.pBuf, sizeof(float), BUFFER_SIZE,
								pOutFile_s_o);
						outBuf_s_o.index = 0;
					}

					// mul fac files
					outBuf_m.pBuf[outBuf_m.index] =
							pInBuf_m[i].pBuf[pInBuf_m[i].index];
					pInBuf_m[i].index++;
					if (pInBuf_m[i].index == BUFFER_SIZE) {
						pInBuf_m[i].read_count = fread(pInBuf_m[i].pBuf,
								sizeof(float),
								BUFFER_SIZE, pInFile_m[i]);
						pInBuf_m[i].index = 0;
					}
					outBuf_m.index++;
					if (outBuf_m.index == BUFFER_SIZE) {
						fwrite(outBuf_m.pBuf, sizeof(float), BUFFER_SIZE,
								pOutFile_m);
						outBuf_m.index = 0;
					}

					// mul fac original files
					outBuf_m_o.pBuf[outBuf_m_o.index] =
							pInBuf_m_o[i].pBuf[pInBuf_m_o[i].index];
					pInBuf_m_o[i].index++;
					if (pInBuf_m_o[i].index == BUFFER_SIZE) {
						pInBuf_m_o[i].read_count = fread(pInBuf_m_o[i].pBuf,
								sizeof(float),
								BUFFER_SIZE, pInFile_m_o[i]);
						pInBuf_m_o[i].index = 0;
					}
					outBuf_m_o.index++;
					if (outBuf_m_o.index == BUFFER_SIZE) {
						fwrite(outBuf_m_o.pBuf, sizeof(float), BUFFER_SIZE,
								pOutFile_m_o);
						outBuf_m_o.index = 0;
					}

					// attn fac files
					outBuf_a.pBuf[outBuf_a.index] =
							pInBuf_a[i].pBuf[pInBuf_a[i].index];
					pInBuf_a[i].index++;
					if (pInBuf_a[i].index == BUFFER_SIZE) {
						pInBuf_a[i].read_count = fread(pInBuf_a[i].pBuf,
								sizeof(float),
								BUFFER_SIZE, pInFile_a[i]);
						pInBuf_a[i].index = 0;
					}
					outBuf_a.index++;
					if (outBuf_a.index == BUFFER_SIZE) {
						fwrite(outBuf_a.pBuf, sizeof(float), BUFFER_SIZE,
								pOutFile_a);
						outBuf_a.index = 0;
					}
				} else { // feof
					run_fread.at(i) = false;
					fclose(pInFile_p[i]);
					fclose(pInFile_s[i]);
					fclose(pInFile_s_o[i]);
					fclose(pInFile_m[i]);
					fclose(pInFile_m_o[i]);
					fclose(pInFile_a[i]);
				}
			} // end probability check
		} // end file for loop

		// check feof for all files
		if (all_of(run_fread.begin(), run_fread.end(), [](bool i) {
			return !i;
		})) {
			run = false;
		}

	} // end main loop

	// flush buffers
	fwrite(outBuf.pBuf, sizeof(ids), outBuf.index, pOutFile);
	fwrite(outBuf_s.pBuf, sizeof(float), outBuf_s.index, pOutFile_s);
	fwrite(outBuf_s_o.pBuf, sizeof(float), outBuf_s_o.index, pOutFile_s_o);
	fwrite(outBuf_m.pBuf, sizeof(float), outBuf_m.index, pOutFile_m);
	fwrite(outBuf_m_o.pBuf, sizeof(float), outBuf_m_o.index, pOutFile_m_o);
	fwrite(outBuf_a.pBuf, sizeof(float), outBuf_a.index, pOutFile_a);

	// clean up

	fclose(pOutFile);
	fclose(pOutFile_s);
	fclose(pOutFile_s_o);
	fclose(pOutFile_m);
	fclose(pOutFile_m_o);
	fclose(pOutFile_a);

	delete[] pInBuf_p;
	delete[] pInBuf_s;
	delete[] pInBuf_s_o;
	delete[] pInBuf_m;
	delete[] pInBuf_m_o;
	delete[] pInBuf_a;

	delete[] pInFile_p;
	delete[] pInFile_s;
	delete[] pInFile_s_o;
	delete[] pInFile_m;
	delete[] pInFile_m_o;
	delete[] pInFile_a;

	cout << "Listmode data compilation complete." << endl;

	return 0;

}
