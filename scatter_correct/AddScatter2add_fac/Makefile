bin = bin/AddScatter2add_fac

CXX = g++
CXXFLAGS = -std=c++14

$(bin) : src/AddScatter2add_fac.o
	$(CXX) $(CXXFLAGS) -o $(bin) src/AddScatter2add_fac.o

AddScatter2add_fac.o : src/AddScatter2add_fac.cpp

.PHONY : clean

clean :
	-rm -v $(bin) src/AddScatter2add_fac.o
