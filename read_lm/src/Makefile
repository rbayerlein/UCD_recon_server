bin = ../process_lm_sino_explorer_singlefile_release
inc = ../../include

CXX = g++
CXXFLAGS = -std=c++14

$(bin) : main.o $(inc)/coinc.o
	$(CXX) $(CXXFLAGS) -o $(bin) main.o $(inc)/coinc.o

main.o : main.cpp $(inc)/coinc.h

$(inc)/coinc.o : $(inc)/coinc.cpp $(inc)/coinc.h $(inc)/bankPairLUT.h

.PHONY : clean

clean :
	-rm -v $(bin) main.o $(inc)/coinc.o
