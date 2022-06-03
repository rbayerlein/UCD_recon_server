bin = process_lm_sino_explorer_singlefile_release
bin_local = process_lm_sino_explorer_local
inc = ../include

CXX = g++
CXXFLAGS = -std=c++14

all : $(bin) $(bin_local)
	
$(bin) : src/main.o $(inc)/coinc.o $(inc)/Subsample.o
	$(CXX) $(CXXFLAGS) -o $(bin) src/main.o $(inc)/coinc.o $(inc)/Subsample.o

$(bin_local) : src/main_local.o $(inc)/coinc.o $(inc)/Subsample.o
	$(CXX) $(CXXFLAGS) -o $(bin_local) src/main_local.o $(inc)/coinc.o $(inc)/Subsample.o

src/main.o : src/main.cpp $(inc)/coinc.h

src/main_local.o : src/main_local.cpp $(inc)/coinc.h

$(inc)/coinc.o : $(inc)/coinc.cpp $(inc)/coinc.h $(inc)/bankPairLUT.h

$(inc)/Subsample.o : $(inc)/Subsample.cpp $(inc)/Subsample.h

.PHONY : clean

clean :
	-rm -v $(bin) src/main.o src/main_local.o $(inc)/coinc.o $(inc)/Subsample.o
