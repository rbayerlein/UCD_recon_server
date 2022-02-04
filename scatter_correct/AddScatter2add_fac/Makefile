bin = bin/AddScatter2add_fac
bin_5d = bin/AddScatter2add_fac_sino5d
bin_5d_par = bin/as2af_s5_par

CXX = g++
CXXFLAGS = -std=c++14

all : AddScatter2add_fac.o AddScatter2add_fac_sino5d.o as2af_s5_par.o

AddScatter2add_fac.o : src/AddScatter2add_fac.o
	$(CXX) $(CXXFLAGS) -o $(bin) src/AddScatter2add_fac.cpp

AddScatter2add_fac_sino5d.o : src/AddScatter2add_fac_sino5d.o
	$(CXX) $(CXXFLAGS) -o $(bin_5d) src/AddScatter2add_fac_sino5d.cpp

as2af_s5_par.o : src/as2af_s5_par.o
	$(CXX) $(CXXFLAGS) -pthread -o $(bin_5d_par) src/as2af_s5_par.cpp

.PHONY : clean

clean :
	-rm -v $(bin) src/AddScatter2add_fac.o
	-rm -v $(bin_5d) src/AddScatter2add_fac_sino5d.o