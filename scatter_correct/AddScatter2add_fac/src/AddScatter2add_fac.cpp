#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

#define BUFFER_SIZE 65536 // = elements per buffer


using namespace std;

struct ids {
	short txID1, axID1, txID2, axID2, tof;
};

int main(int argc, char **argv) {

	cout << "Starting AddScatter2add_fac executable..." << endl;

	// body of program

	// read in events from lm file

	// read in corresponding correction value from sinogram

	// read in correction value from mul_fac file

	// apply correction and save to new add_fac file

	// re-name original file and save new file under original name

	cout << "done" << endl;

}